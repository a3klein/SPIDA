"""
GPU xy tiling for segmentation: plan -> execute -> stitch.

Backend-agnostic. Takes plain data plus an ``eval_tile`` callable.

This file has segmentation-specific policy, that utilities/tiling does not relate to.

There are two tiling layers, one here on the GPU level and one in masks.py at the CPU level post-segmentation.
* masks.py tiles ONE label array produced by ONE ``model.eval``.
    A cell crossing a tile edge keeps the same integer label in
    both tiles, so the fragments re-merge by a plain dissolve on (ID, z).
* here, each tile is an independent ``model.eval`` call, so labels carry
  no meaning across tiles. Merging must be spatial via IoU.

One outer tile therefore contains many inner tiles. Overlap here exists so a
boundary cell is fully contained in at least one tile and so the IoU match has
something to match on; it is also what drives the redundancy cost.

Tissue masking is applied three ways
  1. the grid covers the tissue bounding box rather than the whole mosaic;
  2. tiles with no tissue overlap never reach the GPU (``gate_tiles``);
  3. within a tile, out-of-mask pixels are zeroed by the reader so half-covered
     tiles cannot contribute cells from outside the tissue.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from spida.io.images import DEFAULT_STRIP_ROWS, read_stains
from ...utilities.tiling import (
    DEFAULT_TILE_MAX,
    DEFAULT_TILE_MIN,
    Tile,
    bbox_of_mask,
    build_grid,
    mask_integral_image,
    optimize_grid,
    rect_sums,
    redundancy,
    score_grid,
    tile_coverage,
)

logger = logging.getLogger(__package__)

__all__ = [
    # re-exported: the generic geometry, so a caller planning a grid needs one
    # import. IO names are deliberately absent -- see the import note above.
    "DEFAULT_STRIP_ROWS", "DEFAULT_TILE_MAX", "DEFAULT_TILE_MIN", "Tile",
    "bbox_of_mask", "build_grid", "mask_integral_image", "optimize_grid",
    "read_stains", "rect_sums", "redundancy", "score_grid", "tile_coverage",
    # defined here
    "DEFAULT_GPU_THREADS", "DEFAULT_STITCH_CONTAINMENT",
    "DEFAULT_STITCH_IOU", "PolygonStitcher", "UnionFind",
    "gate_tiles", "plan_grid", "polygon_iou", "polygon_overlap", "run_tiled",
    "stitch_polygons",
    "tile_polygons_to_mosaic",
]


def gate_tiles(tiles: list[Tile], mask: np.ndarray | None,
               min_coverage: float = 0.0) -> tuple[list[Tile], list[Tile]]:
    """Split tiles into (to segment, skipped) by tissue coverage.

    ``min_coverage=0.0`` keeps any tile touching tissue, which is the safe
    default: a tile with a single tissue pixel might still hold part of a cell
    whose centre is in a neighbour.
    """
    if mask is None:
        return list(tiles), []
    keep, skip = [], []
    for tile in tiles:
        (keep if tile_coverage(mask, tile) > min_coverage else skip).append(tile)
    logger.info("tissue gate: %d/%d tiles reach the GPU (%d skipped)",
                len(keep), len(tiles), len(skip))
    return keep, skip


def plan_grid(
    mosaic_shape: tuple[int, int],
    tile: int | tuple[int, int],
    overlap: int,
    mask: np.ndarray | None = None,
    min_coverage: float = 0.0,
) -> dict:
    """Build the full tiling plan and report its cost.

    Returns the bbox, every tile, the gated subsets and the derived metrics, so
    the caller can log a complete picture (and record it in ``meta.json``) before
    committing GPU time.
    """
    if mask is not None:
        bbox = bbox_of_mask(mask)
    else:
        bbox = (0, mosaic_shape[0], 0, mosaic_shape[1])

    tiles = build_grid(bbox, tile, overlap)
    keep, skip = gate_tiles(tiles, mask, min_coverage)

    y0, y1, x0, x1 = bbox
    bbox_area = (y1 - y0) * (x1 - x0)
    mosaic_area = mosaic_shape[0] * mosaic_shape[1]
    partial = 0
    if mask is not None:
        partial = sum(1 for t in keep if tile_coverage(mask, t) < 1.0)

    plan = {
        "bbox": bbox,
        "tiles": tiles,
        "to_segment": keep,
        "skipped": skip,
        "tile_size": tile,
        "overlap": overlap,
        "grid": (len({t.gi for t in tiles}), len({t.gj for t in tiles})),
        "n_tiles": len(tiles),
        "n_segment": len(keep),
        "n_skipped": len(skip),
        "n_partial": partial,
        "redundancy": redundancy(tiles, bbox),
        "bbox_frac_of_mosaic": bbox_area / mosaic_area,
        # pixels actually pushed through the GPU, relative to an unmasked grid
        "compute_frac": (len(keep) / len(tiles)) if tiles else 0.0,
    }
    # %s, not %d, for `tile`: it is an int OR a (tile_y, tile_x) tuple, and a
    # tuple against %d raises inside logging — which the logging module swallows,
    # so the line silently never appears.
    logger.info(
        "tiling plan: bbox=%s (%.1f%% of mosaic), grid=%dx%d, tile=%s overlap=%d, "
        "%d tiles -> %d to segment (%d partial, %d skipped), redundancy=%.2f",
        bbox, plan["bbox_frac_of_mosaic"] * 100, *plan["grid"], tile, overlap,
        plan["n_tiles"], plan["n_segment"], partial, plan["n_skipped"],
        plan["redundancy"],
    )
    return plan



# Stitching:
# Two fragments are declared the same cell as soon as 1 shared z-plane clears the IoU threshold.
# Polygon-based, not mask-based. The trade-off is that geometric IoU on simplified polygons isn't identical to voxel IoU.
DEFAULT_STITCH_IOU = 0.15

# Added a containment criteria to stitching
# Cells are stitched if one is more than 50% contained in the other (mirrors VPT rules)
DEFAULT_STITCH_CONTAINMENT = 0.5


class UnionFind:
    """Minimal union-find over hashable labels, with path compression."""

    def __init__(self) -> None:
        self.parent: dict[int, int] = {}

    def find(self, x: int) -> int:
        root = self.parent.setdefault(x, x)
        while root != self.parent[root]:
            root = self.parent[root]
        while self.parent[x] != root:            # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def polygon_overlap(a, b) -> tuple[float, float]:
    """``(iou, containment)`` for two geometries; ``(0.0, 0.0)`` if disjoint."""

    if not a.intersects(b):
        return 0.0, 0.0
    inter = a.intersection(b).area
    if inter <= 0:
        return 0.0, 0.0
    area_a, area_b = a.area, b.area
    union = area_a + area_b - inter
    smaller = min(area_a, area_b)
    return (inter / union if union > 0 else 0.0,
            inter / smaller if smaller > 0 else 0.0)


def polygon_iou(a, b) -> float:
    """Intersection-over-union only."""
    return polygon_overlap(a, b)[0]


def _overlap_rect(a: "Tile", b: "Tile"):
    """Shared rectangle of two tiles as ``(y0, y1, x0, x1)``, or None."""
    y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
    x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
    if y1 <= y0 or x1 <= x0:
        return None
    return y0, y1, x0, x1


class PolygonStitcher:
    """Incremental cross-tile stitcher: matches pairs as soon as both tiles land.

    Built for the constraint that **nothing but GPU time should extend the
    run**. A batch stitch after the last tile would add its whole cost to the
    wall clock; here a neighbour pair is matched the moment the second of the two
    tiles has been polygonised, so almost all matching happens while the GPU is
    still working on later tiles. Only the final relabel + fragment collapse is
    left for the end, and both are O(polygons) with no geometry predicates.

    Two things make the per-pair work cheap:

    * **Per-tile indices, never a global frame.** The previous batch version
      filtered a concatenated frame with ``geometry.intersects(region)`` for every
      pair, i.e. O(all polygons) shapely calls per pair. Here each tile keeps its
      own ``z -> (geoms, gids, bounds, STRtree)`` index, so a pair only ever
      touches its own two tiles.
    * **Bounds arithmetic before geometry.** Restricting to the overlap rect is a
      vectorised numpy comparison on the bounds array.

    Thread-safe: ``add_tile`` may be called concurrently from the post pool. The
    lock covers only bookkeeping (id allocation, the ready set, union-find
    writes); the expensive geometry runs outside it.
    """

    def __init__(self, tiles: list["Tile"], *,
                 iou_threshold: float = DEFAULT_STITCH_IOU,
                 containment_threshold: float = DEFAULT_STITCH_CONTAINMENT,
                 merge: str = "union") -> None:
        import threading

        if merge not in ("largest", "union"):
            raise ValueError(f"merge must be 'largest' or 'union', got {merge!r}")
        self.tiles = list(tiles)
        self.iou_threshold = iou_threshold
        # >1.0 makes disables containment test.
        self.containment_threshold = (containment_threshold
                                      if 0.0 < containment_threshold <= 1.0
                                      else 2.0)
        self.merge = merge
        self._by_key = {(t.gi, t.gj): t for t in self.tiles}
        self._frames: dict[tuple[int, int], object] = {}
        self._index: dict[tuple[int, int], dict] = {}
        self._uf = UnionFind()
        self._next_id = 1
        self._lock = threading.Lock()
        self.n_containment_merges = 0
        self.n_merges = 0
        self.n_pairs = 0

    # -- ingestion -------------------------------------------------------
    def add_tile(self, key: tuple[int, int], gdf) -> None:
        """Register one tile's polygons and match against any ready neighbours."""
        if gdf is None or len(gdf) == 0:
            with self._lock:
                self._frames[key] = None
                self._index[key] = {}
            return

        with self._lock: # id allocation only
            base = self._next_id
            self._next_id += int(gdf["ID"].max()) + 1

        g = gdf.copy()
        g["_gid"] = g["ID"].to_numpy().astype(np.int64) + base
        index = self._build_index(g) # outside the lock: pure CPU

        with self._lock:
            self._frames[key] = g
            self._index[key] = index
            for gid in g["_gid"].to_numpy():
                self._uf.find(int(gid))
            # only when the latter tile lands so unordered pairs fire just once.
            ready = [
                nb for nb in self._neighbours(key)
                if nb in self._index and self._index[nb]
            ]
            snapshot = {nb: self._index[nb] for nb in ready}

        for nb in ready:                      # geometry outside the lock
            pairs = self._match(key, index, nb, snapshot[nb])
            if pairs:
                with self._lock:
                    for a, b in pairs:
                        self._uf.union(a, b)
                    self.n_merges += len(pairs)
            with self._lock:
                self.n_pairs += 1

    def _neighbours(self, key):
        gi, gj = key
        for dgi, dgj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            nb = (gi + dgi, gj + dgj)
            if nb in self._by_key:
                yield nb

    @staticmethod
    def _build_index(g) -> dict:
        """``z -> (geoms, gids, bounds, STRtree)`` for one tile."""
        from shapely.strtree import STRtree

        out = {}
        bounds_all = g.bounds.to_numpy()
        z_all = g["z"].to_numpy()
        gid_all = g["_gid"].to_numpy()
        geom_all = g.geometry.to_numpy()
        for z in np.unique(z_all):
            sel = z_all == z
            geoms = list(geom_all[sel])
            out[int(z)] = (geoms, gid_all[sel], bounds_all[sel], STRtree(geoms))
        return out

    def _match(self, key_a, idx_a, key_b, idx_b) -> list[tuple[int, int]]:
        """Candidate unions between two neighbouring tiles, any-z."""
        rect = _overlap_rect(self._by_key[key_a], self._by_key[key_b])
        if rect is None:
            return []
        y0, y1, x0, x1 = rect
        pairs: list[tuple[int, int]] = []

        for z, (geoms_a, gids_a, bounds_a, _tree_a) in idx_a.items():
            entry_b = idx_b.get(z)
            if entry_b is None:
                continue
            geoms_b, gids_b, _bounds_b, tree_b = entry_b
            # bounds arithmetic first; polygon coords are (x, y)
            keep = ((bounds_a[:, 0] < x1) & (bounds_a[:, 2] > x0) &
                    (bounds_a[:, 1] < y1) & (bounds_a[:, 3] > y0))
            if not keep.any():
                continue
            for i in np.flatnonzero(keep):
                geom = geoms_a[i]
                for j in tree_b.query(geom):
                    iou, containment = polygon_overlap(geom, geoms_b[j])
                    if iou >= self.iou_threshold:
                        pairs.append((int(gids_a[i]), int(gids_b[j])))
                    elif containment >= self.containment_threshold:
                        # A sliver: one side is a thin lip of the other, cut off
                        # at the tile edge. Low IoU by construction, so the test
                        # above cannot see it. See DEFAULT_STITCH_CONTAINMENT.
                        self.n_containment_merges += 1
                        pairs.append((int(gids_a[i]), int(gids_b[j])))
        return pairs

    # -- finalisation ----------------------------------------------------
    def finish(self):
        """Resolve union-find, collapse fragments, and renumber densely."""
        import geopandas as gpd
        import pandas as pd

        frames = [f for f in self._frames.values() if f is not None and len(f)]
        if not frames:
            logger.warning("PolygonStitcher: no polygons to stitch")
            return gpd.GeoDataFrame(
                {"ID": np.array([], dtype=np.int64),
                 "z": np.array([], dtype=np.int64),
                 "Geometry": gpd.GeoSeries([], dtype="geometry")},
                geometry="Geometry", crs=None,
            )
        combined = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True),
                                    geometry="Geometry", crs=None)
        combined["_final"] = [self._uf.find(int(g)) for g in combined["_gid"]]
        logger.info(
            "stitch: %d fragments from %d tiles, %d pairs compared, %d merges "
            "-> %d cells", len(combined), len(frames), self.n_pairs,
            self.n_merges, combined["_final"].nunique(),
        )
        # How the merges split across the two criteria is tuning detail, not
        # something an operator needs on every run.
        logger.debug("stitch: %d of %d merges came from containment rather than "
                     "IoU (threshold %.2f)", self.n_containment_merges,
                     self.n_merges, self.containment_threshold)

        out = _collapse_fragments(combined, self.merge)
        codes, _ = out["_final"].factorize(sort=True)
        out = out.assign(ID=codes.astype(np.int64) + 1)
        out = out[["ID", "z", "Geometry"]].reset_index(drop=True)
        return gpd.GeoDataFrame(out, geometry="Geometry", crs=None)


def stitch_polygons(
    tile_polys: dict,
    tiles: list["Tile"],
    *,
    iou_threshold: float = DEFAULT_STITCH_IOU,
    containment_threshold: float = DEFAULT_STITCH_CONTAINMENT,
    merge: str = "union",
):
    """Batch convenience wrapper over :class:`PolygonStitcher`.

    Equivalent to feeding every tile in and calling ``finish``. The pipeline uses
    the incremental class directly so matching overlaps GPU work; this exists for
    callers that already have all tiles in hand.
    """
    stitcher = PolygonStitcher(
        tiles, iou_threshold=iou_threshold,
        containment_threshold=containment_threshold, merge=merge)
    for key in sorted(tile_polys):
        stitcher.add_tile(key, tile_polys[key])
    return stitcher.finish()


def _collapse_fragments(combined, merge: str):
    """One geometry per (final id, z)."""
    import geopandas as gpd
    import pandas as pd

    dup = combined.duplicated(subset=["_final", "z"], keep=False)
    singles = combined[~dup]
    multi = combined[dup]
    if len(multi) == 0:
        return combined

    if merge == "largest":
        idx = multi.geometry.area.groupby(
            [multi["_final"], multi["z"]]).idxmax().to_numpy()
        kept = multi.loc[idx]
    else:
        kept = multi.dissolve(by=["_final", "z"]).reset_index()

    return gpd.GeoDataFrame(
        pd.concat([singles, kept], ignore_index=True),
        geometry="Geometry", crs=None,
    )


# ---------------------------------------------------------------------------
# 3-stage pipeline
#
# The GPU thread runs NOTHING but the model call. Every other cost — reading,
# downscaling, tissue-zeroing, polygonising, coordinate transforms — happens on
# worker threads either side of it. It is the invariant a multi-GPU-per-node pool 
# needs, and it means CPU-side work like sub-tiling can be added later for free.
#
# Threads rather than processes: the heavy stages (tifffile IO, numpy block
# means, shapely) all release the GIL, and threads avoid pickling large arrays.
# ---------------------------------------------------------------------------

_SENTINEL = object()
DEFAULT_GPU_THREADS = 2

def _default_workers() -> int:
    """CPUs actually allocated to this job, not the machine's total."""
    import os

    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:                       # pragma: no cover - non-Linux
        return max(1, os.cpu_count() or 1)


def tile_polygons_to_mosaic(
    gdf,
    tile: "Tile",
    scale: int
):
    """
    Scale a tile's polygons back to full resolution and place them in the mosaic.
    
    Parameters
    ----------
    gdf : GeoDataFrame
        Polygons in tile coordinates.
    tile : Tile
        The tile's bounding box in mosaic coordinates.
    scale : int
        The downscale factor applied to the tile when reading it.
    Returns
    -------
    GeoDataFrame
        Polygons in mosaic coordinates.
    
    Notes
    -----
    Polygon coordinates from ``masks.py`` are ``(x, y)`` (it swaps axes before
    contouring), so the translation is x-first.
    """
    if gdf is None or len(gdf) == 0:
        return gdf
    out = gdf.copy()
    out.geometry = out.geometry.affine_transform(
        [scale, 0, 0, scale, tile.x0, tile.y0]
    )
    return out


def run_tiled(
    plan: dict,
    stain_planes: dict[str, list[Path]],
    eval_tile,
    *,
    scale: int = 4,
    mask: np.ndarray | None = None,
    poly_tile_size: int = 1000,
    poly_tile_overlap: int = 200,
    poly_tolerance: float = 0.5,
    n_prep: int | None = None,
    n_post: int | None = None,
    n_gpu_threads: int = DEFAULT_GPU_THREADS,
    queue_depth: int = 2,
    iou_threshold: float = DEFAULT_STITCH_IOU,
    containment_threshold: float = DEFAULT_STITCH_CONTAINMENT,
    merge: str = "union",
    strip_rows: int = DEFAULT_STRIP_ROWS,
):
    """Segment a tiled mosaic and return one stitched polygon GeoDataFrame.

    Parameters:
    plan (dict): output of :func:`plan_grid` — supplies ``to_segment``.
    stain_planes (dict): stain name -> z-ordered plane paths, as built by the
        backend's ``_resolve_stains``. Channel order follows the dict's
        insertion order (nuclear first, cytoplasmic second) -- the same dict the
        whole-image path hands to ``read_stains``, which is what keeps the two
        paths' channel conventions identical by construction.
    eval_tile (callable): ``(C, Y, X, Z) -> masks``. Injected so this module
        stays cellpose-free; the backend passes a closure over its config.
    scale (int): xy downsample applied during the read, and the factor polygons
        are scaled back up by.
    mask (np.ndarray | None): full-resolution tissue mask; zeroes out-of-tissue
        pixels before the block mean.
    n_gpu_threads (int): concurrent ``eval_tile`` callers sharing ONE model. See
        :data:`DEFAULT_GPU_THREADS` for why this is >1.
    queue_depth (int): bounded queues cap how many tiles are in flight, which is
        what keeps peak RAM proportional to a few tiles rather than the grid.

    Returns a GeoDataFrame of ``ID``, ``z``, ``Geometry`` in full-resolution
    mosaic pixel coordinates — the same contract as the single-call path.
    """
    import queue
    import threading
    import time

    from .backends.masks import masks_to_geodataframe

    tiles = list(plan["to_segment"])
    if not tiles:
        raise ValueError("plan has no tiles to segment")

    n_prep = n_prep or max(2, min(4, _default_workers()))
    n_post = n_post or max(2, min(4, _default_workers()))
    eval_q: queue.Queue = queue.Queue(maxsize=queue_depth)
    post_q: queue.Queue = queue.Queue(maxsize=queue_depth)
    tile_q: queue.Queue = queue.Queue()
    for t in tiles:
        tile_q.put(t)

    # Stitching with neighbor tiles.
    stitcher = PolygonStitcher(
        tiles, iou_threshold=iou_threshold,
        containment_threshold=containment_threshold, merge=merge)
    errors: list[BaseException] = []
    stop = threading.Event()

    def _record(exc: BaseException) -> None:
        errors.append(exc)
        stop.set()

    def _put(q: queue.Queue, item) -> bool:
        # Block while full, but give up once stop is set: the consumers may be gone.
        while not stop.is_set():
            try:
                q.put(item, timeout=0.5)
                return True
            except queue.Full:
                continue
        return False

    def prep_worker() -> None:
        try:
            while not stop.is_set():
                try:
                    tile = tile_q.get_nowait()
                except queue.Empty:
                    return
                img = read_stains(                    # (C, Y, X, Z)
                    stain_planes,
                    window=(tile.y0, tile.y1, tile.x0, tile.x1),
                    downscale=scale, strip_rows=strip_rows, mask=mask,
                )
                _put(eval_q, (tile, img))
        except BaseException as exc:                  # noqa: BLE001
            _record(exc)

    def post_worker() -> None:
        try:
            while not stop.is_set():
                try:
                    item = post_q.get(timeout=0.5)
                except queue.Empty:
                    continue
                if item is _SENTINEL:
                    post_q.task_done()
                    return
                tile, masks = item
                try:
                    # n_workers=1 is REQUIRED here. this runs on a worker thread, and 
                    # mp.Pool's fork() from a multi-threaded process can deadlock the child.
                    gdf = masks_to_geodataframe(
                        masks, tolerance=poly_tolerance,
                        tile_size=poly_tile_size, overlap=poly_tile_overlap,
                        n_workers=1,
                    )
                    gdf = tile_polygons_to_mosaic(gdf, tile, scale)
                    stitcher.add_tile((tile.gi, tile.gj), gdf)
                finally:
                    post_q.task_done()
        except BaseException as exc:                  # noqa: BLE001
            _record(exc)

    total = len(tiles)
    progress = {"done": 0}
    progress_lock = threading.Lock()

    # Where the GPU stage's wall-clock actually goes.
    # If wait_s is large, prep is not keeping up and GPU idle time is high.
    # If block_s is large, post is not keeping up and GPU idle time is high.
    # If eval_s dominates, the GPU thread is busy inside model.eval the whole time, GPU idle time is a part of the method itself.
    stage_s = {"wait_s": 0.0, "eval_s": 0.0, "block_s": 0.0}

    def gpu_worker() -> None:
        """The ONLY stage touching the GPU. Runs nothing but the model call."""
        try:
            while not stop.is_set():
                t_wait = time.perf_counter()
                try:
                    item = eval_q.get(timeout=0.5)
                except queue.Empty:
                    with progress_lock:
                        stage_s["wait_s"] += time.perf_counter() - t_wait
                    continue
                waited = time.perf_counter() - t_wait
                if item is _SENTINEL:
                    return
                tile, img = item
                t_eval = time.perf_counter()
                masks = eval_tile(img)
                evaled = time.perf_counter() - t_eval
                t_block = time.perf_counter()
                with progress_lock:
                    progress["done"] += 1
                    done = progress["done"]
                    stage_s["wait_s"] += waited
                    stage_s["eval_s"] += evaled
                    snapshot = dict(stage_s)
                # Log every tile. These runs can die: spot preemption,
                # OOM on an undersized node, wall-clock timeout. Polygons are
                # only written after `finish()`, so this log is the sole record
                # of progress and of the eval/wait/post split.
                busy = snapshot["eval_s"] + snapshot["wait_s"] + snapshot["block_s"]
                logger.info(
                    "tiled: %d/%d tiles segmented (tile %d,%d) | this tile "
                    "wait=%.1fs eval=%.1fs | cumulative eval=%.0f%% wait=%.0f%% "
                    "post-block=%.0f%%",
                    done, total, tile.gi, tile.gj, waited, evaled,
                    100 * snapshot["eval_s"] / busy if busy else 0,
                    100 * snapshot["wait_s"] / busy if busy else 0,
                    100 * snapshot["block_s"] / busy if busy else 0,
                )
                _put(post_q, (tile, masks))
                with progress_lock:
                    stage_s["block_s"] += time.perf_counter() - t_block
        except BaseException as exc:                  # noqa: BLE001
            _record(exc)

    preps = [threading.Thread(target=prep_worker, name=f"prep{i}", daemon=True)
             for i in range(n_prep)]
    gpus = [threading.Thread(target=gpu_worker, name=f"gpu{i}", daemon=True)
            for i in range(max(1, n_gpu_threads))]
    posts = [threading.Thread(target=post_worker, name=f"post{i}", daemon=True)
             for i in range(n_post)]
    for th in (*preps, *gpus, *posts):
        th.start()

    try:
        # Shut the stages down in pipeline order so nothing is dropped.
        # Sentinels are only sent once the upstream stage has actually finished.
        for th in preps:
            th.join()
        for _ in gpus:
            _put(eval_q, _SENTINEL)  # unbounded put hung here when a gpu thread failed
        for th in gpus:
            th.join()
    except BaseException as exc:                      # noqa: BLE001
        _record(exc)
    finally:
        for _ in posts:
            _put(post_q, _SENTINEL)
        for th in posts:
            th.join(timeout=120)
        stop.set()
        for th in preps:
            th.join(timeout=10)
        for th in gpus:
            th.join(timeout=10)

    done = progress["done"]

    if errors:
        raise RuntimeError(
            f"tiled segmentation failed after {done}/{total} tiles"
        ) from errors[0]
    if done != total:
        raise RuntimeError(
            f"tiled segmentation ended early: {done}/{total} tiles processed"
        )

    # Everything except the final relabel + fragment collapse already happened
    # while the GPU was busy.
    busy = stage_s["eval_s"] + stage_s["wait_s"] + stage_s["block_s"]
    logger.info(
        "tiled: GPU-stage time breakdown over %d tiles x %d thread(s): "
        "in-eval %.0fs (%.0f%%), waiting for prep %.0fs (%.0f%%), "
        "blocked on post %.0fs (%.0f%%)",
        done, max(1, n_gpu_threads),
        stage_s["eval_s"], 100 * stage_s["eval_s"] / busy if busy else 0,
        stage_s["wait_s"], 100 * stage_s["wait_s"] / busy if busy else 0,
        stage_s["block_s"], 100 * stage_s["block_s"] / busy if busy else 0,
    )
    logger.debug("tiled: finalising stitch (%d pairs already matched in-flight)",
                 stitcher.n_pairs)
    return stitcher.finish()
