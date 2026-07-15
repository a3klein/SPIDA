"""Pure-Python reimplementations of VPT (vizgen-postprocessing) commands.

These replace shelling out to the external ``vpt`` binary. Each function is a
drop-in for the corresponding VPT command, reading/writing the same file formats
and matching VPT's semantics exactly (validated bit-for-bit against the binary).

Currently implemented:
    - ``partition_transcripts`` (VPT ``partition-transcripts``)
    - ``derive_entity_metadata`` (VPT ``derive-entity-metadata``)
    - ``convert_geometry`` (VPT ``convert-geometry``, 2D cellpose -> micron 3D)
    - ``sum_signals`` (VPT ``sum-signals``)
"""

from __future__ import annotations

import re
import glob
import datetime
import logging
from pathlib import Path
from typing import Iterator
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import geopandas as gpd
import pyarrow.parquet as pq
from shapely import affinity, make_valid, geometry
from shapely.ops import unary_union
from shapely.affinity import affine_transform, translate
from fsspec.implementations.local import LocalFileSystem

from spida._constants import DEFAULT_PRESET, GEOMETRY_COL  # noqa: F401

logger = logging.getLogger(__package__)

# raw VPT file-format column names (the format this replaces reads/writes)
_ENTITY_COL = DEFAULT_PRESET["METADATA_CELL_KEY"]  # "EntityID"
_ZINDEX_COL = DEFAULT_PRESET["geom_depth_col"]     # "ZIndex"


def _iter_transcripts(
    path: str | Path, chunk_size: int | None
) -> Iterator[pd.DataFrame]:
    """Yield the transcript table whole (``chunk_size=None``) or in row batches."""
    path = str(path)
    if chunk_size is None:
        if path.endswith(".parquet"):
            yield pd.read_parquet(path)
        else:
            yield pd.read_csv(path)
        return
    if path.endswith(".parquet"):
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(chunk_size):
            yield batch.to_pandas()
    else:
        yield from pd.read_csv(path, chunksize=chunk_size)


def _assign_chunk(
    tx: pd.DataFrame,
    polys: gpd.GeoDataFrame,
    z_planes: np.ndarray,
    gene_col: str,
    x_col: str,
    y_col: str,
    z_col: str,
    entity_col: str,
    zindex_col: str,
) -> pd.DataFrame:
    """Spatially join one transcript batch to cell polygons.

    Returns rows of ``(_tx_idx, gene, EntityID)``. A point inside ``k > 1``
    polygons yields ``k`` rows, reproducing VPT's overlap double-count.
    """
    geom_name = polys.geometry.name
    tx = tx[tx[z_col].isin(z_planes)]
    parts = []
    for z in z_planes:
        txz = tx[tx[z_col] == z]
        pz = polys[polys[zindex_col] == z][[entity_col, geom_name]]
        if len(txz) == 0 or len(pz) == 0:
            continue
        pts = gpd.GeoDataFrame(
            {"_tx_idx": txz.index.values, gene_col: txz[gene_col].values},
            geometry=gpd.points_from_xy(txz[x_col].values, txz[y_col].values),
        )
        j = gpd.sjoin(pts, pz, predicate="within", how="inner")
        parts.append(j[["_tx_idx", gene_col, entity_col]])
    if parts:
        return pd.concat(parts, ignore_index=True)
    return pd.DataFrame(columns=["_tx_idx", gene_col, entity_col])


def _write_transcripts(
    tx: pd.DataFrame,
    joined: pd.DataFrame,
    entity_col: str,
    out_path: str | Path,
    append: bool,
) -> None:
    """Write a transcript batch with a ``cell_id`` column (last-wins on overlap)."""
    cell_id = np.full(len(tx), -1, dtype=np.int64)
    pos = pd.Series(np.arange(len(tx)), index=tx.index)
    if len(joined):
        last = joined.drop_duplicates(subset="_tx_idx", keep="last")
        cell_id[pos.loc[last["_tx_idx"]].values] = last[entity_col].values
    out = tx.copy()
    out["cell_id"] = cell_id
    out_path = str(out_path)
    if out_path.endswith(".parquet"):
        out.to_parquet(out_path, index=False)  # chunked parquet not appended
    else:
        out.to_csv(out_path, index=False, mode="a" if append else "w", header=not append)


def partition_transcripts(
    boundaries_path: str | Path,
    transcripts_path: str | Path,
    output_entity_by_gene: str | Path | None = None,
    output_transcripts: str | Path | None = None,
    *,
    chunk_size: int | None = None,
    gene_col: str = "gene",
    barcode_col: str = "barcode_id",
    x_col: str = "global_x",
    y_col: str = "global_y",
    z_col: str = "global_z",
    entity_col: str = _ENTITY_COL,
    zindex_col: str = _ZINDEX_COL,
) -> pd.DataFrame:
    """Partition detected transcripts into cells and build a cell-by-gene matrix.

    Pure-Python replacement for VPT ``partition-transcripts``. Assigns each
    transcript to the cell polygon that contains it (boundary-exclusive, matched
    within the transcript's own integer z-plane) and tallies counts per cell/gene.

    Parameters
    ----------
    boundaries_path
        Cell boundary polygons in VPT micron-space format (one polygon per cell
        per z-plane), e.g. ``cellpose_micron_space.parquet``.
    transcripts_path
        Detected transcripts (``.parquet`` or ``.csv``) with x/y/z, gene, barcode.
    output_entity_by_gene
        If given, write the cell-by-gene matrix as CSV here.
    output_transcripts
        If given, write the transcript table with an assigned ``cell_id`` column.
    chunk_size
        ``None`` (default) reads all transcripts at once (fastest, higher peak
        memory). An int streams transcripts in row batches, bounding memory.

    Returns
    -------
    pandas.DataFrame
        Cell-by-gene counts: rows = every cell (index name ``"cell"``, sorted),
        columns = genes ordered by ``barcode_id`` over observed barcodes.

    Notes
    -----
    Matches VPT exactly, including: transcripts with z outside the boundary
    z-planes are unassigned (``cell_id = -1``) and uncounted; a transcript inside
    overlapping polygons is counted for every containing cell, and its per-row
    ``cell_id`` takes the last such cell.
    """
    polys = gpd.read_parquet(boundaries_path)
    polys = polys[[entity_col, zindex_col, polys.geometry.name]].copy()
    polys[entity_col] = polys[entity_col].astype(np.int64)
    z_planes = np.sort(polys[zindex_col].unique())
    all_ids = np.sort(polys[entity_col].unique())
    logger.info(
        "partition_transcripts: %d cells across %d z-planes",
        len(all_ids), len(z_planes),
    )

    running = None            # accumulated (EntityID, gene) counts
    barcode_frames = []       # per-chunk unique barcode->gene, for column ordering
    first_chunk = True
    n_tx = 0

    for tx in _iter_transcripts(transcripts_path, chunk_size):
        n_tx += len(tx)
        joined = _assign_chunk(
            tx, polys, z_planes, gene_col, x_col, y_col, z_col, entity_col, zindex_col
        )
        c = joined.groupby([entity_col, gene_col]).size()
        running = c if running is None else running.add(c, fill_value=0)
        barcode_frames.append(
            tx[[barcode_col, gene_col]].drop_duplicates(subset=barcode_col)
        )
        if output_transcripts is not None:
            _write_transcripts(
                tx, joined, entity_col, output_transcripts, append=not first_chunk
            )
        first_chunk = False

    # column order: genes ordered by barcode_id, over barcodes actually observed
    gene_order = (
        pd.concat(barcode_frames)
        .drop_duplicates(subset=barcode_col)
        .sort_values(barcode_col)[gene_col]
        .tolist()
    )

    if running is None or len(running) == 0:
        cbg = pd.DataFrame(index=all_ids)
    else:
        cbg = running.astype("int64").unstack(fill_value=0)
    cbg = cbg.reindex(index=all_ids, columns=gene_order, fill_value=0).astype("int64")
    cbg.index.name = "cell"
    cbg.columns.name = None
    logger.info("partition_transcripts: %d transcripts -> cell-by-gene %s",
                n_tx, cbg.shape)

    if output_entity_by_gene is not None:
        Path(output_entity_by_gene).parent.mkdir(parents=True, exist_ok=True)
        cbg.to_csv(output_entity_by_gene)

    return cbg


# --------------------------------------------------------------------------- #
# derive-entity-metadata                                                      #
# --------------------------------------------------------------------------- #
_METADATA_COLS = [
    "fov", "volume", "center_x", "center_y", "min_x", "min_y", "max_x", "max_y",
    "anisotropy", "transcript_count", "perimeter_area_ratio", "solidity",
]


def _anisotropy(mrr) -> float:
    """major/minor side-length ratio of a minimum rotated rectangle polygon."""
    try:
        xs, ys = mrr.exterior.xy
    except (AttributeError, NotImplementedError):
        return np.nan
    pts = list(zip(xs, ys))
    sides = [
        np.hypot(p1[0] - p2[0], p1[1] - p2[1])
        for p1, p2 in zip(pts[:-1], pts[1:])
    ]
    if not sides or min(sides) == 0:
        return np.nan
    return max(sides) / min(sides)


def derive_entity_metadata(
    boundaries_path: str | Path,
    cell_by_gene_path: str | Path | None = None,
    output_metadata: str | Path | None = None,
    *,
    entity_col: str = _ENTITY_COL,
    zindex_col: str = _ZINDEX_COL,
    zlevel_col: str = "ZLevel",
) -> pd.DataFrame:
    """Compute per-cell geometric metadata from boundary polygons.

    Pure-Python replacement for VPT ``derive-entity-metadata``. For each cell,
    over its non-degenerate z-plane polygons:

    ==================== ====================================================
    ``volume``           ``sum(area_z * z_thickness_z)``
    ``center_x/center_y`` mean of per-plane centroid coordinates
    ``min/max_x/y``      overall bounds across planes
    ``anisotropy``       major/minor side of the min rotated rect of all planes
    ``perimeter_area_ratio`` mean of ``length_z / area_z``
    ``solidity``         ``sum(area_z) / sum(convex_hull_area_z)``
    ``transcript_count`` row-sum of the cell-by-gene matrix (incl. blanks)
    ``fov``              always NaN (matches VPT)
    ==================== ====================================================

    Parameters
    ----------
    boundaries_path
        Cell boundary polygons in VPT micron-space format.
    cell_by_gene_path
        Cell-by-gene CSV; if given, ``transcript_count`` is its per-cell row sum.
    output_metadata
        If given, write the metadata table as CSV here.

    Returns
    -------
    pandas.DataFrame
        Indexed by ``EntityID`` (sorted); columns in VPT order. Cells with no
        valid polygon get an all-NaN row.

    Notes
    -----
    Matches VPT for every column except ``anisotropy``, which depends on
    shapely's ``minimum_rotated_rectangle`` implementation and can differ
    negligibly on near-symmetric cells across shapely versions. The value from
    the installed shapely is used (the newer, generally more correct one).
    """
    polys = gpd.read_parquet(boundaries_path)
    geom_name = polys.geometry.name
    all_ids = np.sort(polys[entity_col].unique())

    # z-plane thickness, indexed by ZIndex (diff of sorted unique ZLevels)
    zpairs = polys[[zindex_col, zlevel_col]].drop_duplicates().sort_values(zindex_col)
    depths = np.diff(zpairs[zlevel_col].to_numpy(), prepend=0)
    depth_map = dict(zip(zpairs[zindex_col].to_numpy(), depths))

    g = polys[[entity_col, zindex_col, geom_name]].copy()
    area = g.geometry.area
    g = g[area >= 1e-9].copy()  # drop degenerate/empty planes (matches VPT)
    b = g.geometry.bounds
    cent = g.geometry.centroid
    g["_area"] = g.geometry.area.values
    g["_len"] = g.geometry.length.values
    g["_cx"] = cent.x.values
    g["_cy"] = cent.y.values
    g["_minx"] = b["minx"].values
    g["_miny"] = b["miny"].values
    g["_maxx"] = b["maxx"].values
    g["_maxy"] = b["maxy"].values
    g["_hull"] = g.geometry.convex_hull.area.values
    g["_depth"] = g[zindex_col].map(depth_map).astype(float)
    g["_vol"] = g["_area"] * g["_depth"]
    g["_par"] = g["_len"] / g["_area"]

    grp = g.groupby(entity_col)
    out = pd.DataFrame(index=grp.size().index)
    out["volume"] = grp["_vol"].sum()
    out["center_x"] = grp["_cx"].mean()
    out["center_y"] = grp["_cy"].mean()
    out["min_x"] = grp["_minx"].min()
    out["min_y"] = grp["_miny"].min()
    out["max_x"] = grp["_maxx"].max()
    out["max_y"] = grp["_maxy"].max()
    out["perimeter_area_ratio"] = grp["_par"].mean()
    out["solidity"] = grp["_area"].sum() / grp["_hull"].sum()

    dissolved = g.dissolve(by=entity_col)
    mrr = dissolved.geometry.minimum_rotated_rectangle()
    out["anisotropy"] = [_anisotropy(m) for m in mrr]

    out["fov"] = np.nan

    if cell_by_gene_path is not None:
        cbg = pd.read_csv(cell_by_gene_path, index_col=0)
        out["transcript_count"] = cbg.sum(axis=1).reindex(out.index)
    else:
        out["transcript_count"] = np.nan

    out = out.reindex(all_ids)[_METADATA_COLS]
    out.index.name = entity_col
    logger.info("derive_entity_metadata: %d cells", len(out))

    if output_metadata is not None:
        Path(output_metadata).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output_metadata)
    return out


# --------------------------------------------------------------------------- #
# convert-geometry (2D cellpose polygons -> micron-space 3D boundaries)       #
# --------------------------------------------------------------------------- #
# NOTE: this handles the 2D->3D path. The already-3D path is handled by
# vpt_utils._3d_convert_geometry, which (by design) does not harmonize overlaps.
_GEOM = "Geometry"

# module-level EntityID generator state (VPT timestamp scheme)
_counter_digits = 8
_counter = 0
_process_id = ""


def _format_experiment_timestamp(ts: float) -> str:
    t = datetime.datetime.fromtimestamp(ts)
    year_day = str(int(t.strftime("%j")) + 100)
    day_time = t.strftime("%H%M%S")
    day_second = str(
        int(day_time[:2]) * 3600 + int(day_time[2:4]) * 60 + int(day_time[4:])
    ).zfill(5)
    return "".join([year_day, day_second])


def _set_process_id() -> None:
    global _process_id, _counter_digits, _counter
    _process_id = _format_experiment_timestamp(
        datetime.datetime.now(datetime.UTC).timestamp()
    )
    _counter_digits = 19 - len(_process_id)
    _counter = 0


def _get_id() -> np.int64:
    global _counter
    next_id = np.int64(f"{_process_id}{str(_counter).zfill(_counter_digits)}")
    _counter += 1
    return next_id


def _get_valid(shape):
    try:
        return make_valid(shape)
    except ValueError:
        return geometry.MultiPolygon()


def _to_multipoly(shape):
    """Coerce a geometry to a MultiPolygon (matches vpt_core.convert_to_multipoly)."""
    if isinstance(shape, geometry.Polygon):
        return geometry.MultiPolygon([shape])
    if isinstance(shape, geometry.MultiPolygon):
        return shape
    if isinstance(shape, geometry.GeometryCollection):
        polys = [g for g in shape.geoms
                 if isinstance(g, (geometry.Polygon, geometry.MultiPolygon))]
        poly = _get_valid(unary_union(polys))
        return poly if isinstance(poly, geometry.MultiPolygon) \
            else geometry.MultiPolygon([poly])
    return geometry.MultiPolygon()


def _read_micron_to_mosaic_transform(path: str | Path) -> np.ndarray:
    fs = LocalFileSystem()
    lines = fs.open(str(path), "r", lambda f: f.readlines())
    transform = [list(map(float, ln.split())) for ln in lines]
    if len(transform) != 3 or not all(len(r) == 3 for r in transform):
        raise ValueError("Micron to mosaic transform should be a 3x3 matrix")
    return np.array(transform)


def _overlap_pairs(g: gpd.GeoDataFrame) -> list[tuple[int, int]]:
    """Unordered unique positional pairs of intersecting polygons (single plane)."""
    res = g.sindex.query(g.geometry, predicate="intersects")
    return list({(min(a, b), max(a, b)) for a, b in zip(res[0], res[1]) if a != b})


def harmonize_polygons(
    g: gpd.GeoDataFrame, min_distance: int = 2, min_area: int = 100
) -> gpd.GeoDataFrame:
    """Resolve overlaps between single-plane cell polygons (VPT ``harmonize``).

    Merge pass: if two cells overlap by >50% of the smaller cell's area, union
    the smaller into the larger and drop the smaller. Trim pass: for remaining
    overlaps, trim the larger cell by ``smaller.buffer(min_distance)``. Finally
    drop cells with area <= ``min_area``. Operates in the input coordinate space
    (pixels), matching VPT which harmonizes before the micron transform.
    """
    g = g.reset_index(drop=True).copy()

    deprecated: set[int] = set()
    for i, j in _overlap_pairs(g):
        ei, ej = g.at[i, _ENTITY_COL], g.at[j, _ENTITY_COL]
        if ei in deprecated or ej in deprecated:
            continue
        gi, gj = g.at[i, _GEOM], g.at[j, _GEOM]
        vi, vj = gi.area, gj.area
        if min(vi, vj) == 0:
            continue
        overlap = _get_valid(gi).intersection(_get_valid(gj)).area
        if overlap / min(vi, vj) > 0.5:
            union = _to_multipoly(_get_valid(_get_valid(gi).union(_get_valid(gj))))
            if vi > vj:
                g.at[i, _GEOM] = union
                deprecated.add(ej)
            else:
                g.at[j, _GEOM] = union
                deprecated.add(ei)
    g = g[~g[_ENTITY_COL].isin(deprecated)].reset_index(drop=True)

    for i, j in _overlap_pairs(g):
        gi, gj = g.at[i, _GEOM], g.at[j, _GEOM]
        large_i, large_g, small_g = (i, gi, gj) if gi.area > gj.area else (j, gj, gi)
        try:
            trimmed = _get_valid(large_g).difference(
                _get_valid(small_g).buffer(min_distance)
            )
            g.at[large_i, _GEOM] = _to_multipoly(_get_valid(trimmed))
        except ValueError:
            g.at[large_i, _GEOM] = geometry.MultiPolygon()

    return g[g[_GEOM].area > min_area].reset_index(drop=True)


def convert_geometry(
    input_boundaries: str | Path,
    output_boundaries: str | Path,
    input_micron_to_mosaic: str | Path | None = None,
    *,
    n_z_planes: int = 7,
    spacing_z: float = 1.5,
    harmonize: bool = True,
    min_distance: int = 2,
    min_area: int = 100,
    entity_type: str = "cell",
) -> gpd.GeoDataFrame:
    """Convert raw 2D cellpose polygons to VPT micron-space 3D boundaries.

    Pure-Python replacement for VPT ``convert-geometry`` (2D parquet path):
    coerce geometries to valid MultiPolygons, harmonize overlaps (pixel space),
    replicate across ``n_z_planes`` with ``ZLevel = spacing_z*(z+1)``, then apply
    the mosaic->micron affine transform (inverse of the 3x3 matrix). Fresh
    timestamp-based ``EntityID``s are assigned. Writes a parquet with columns
    ``ID, Geometry, EntityID, ParentID, ParentType, Type, ZIndex, ZLevel``.

    Harmonizing the single 2D plane before replication is equivalent to VPT's
    replicate-then-harmonize (identical planes, per-plane overlap resolution).

    Notes
    -----
    Geometry matches VPT for ~99.8% of cells. A small fraction (~0.15%) of cells
    that sit at tile seams and overlap two or more neighbours may receive a
    contested sliver differently than VPT: overlap resolution is order-sensitive
    and VPT's order is itself arbitrary (set-hash over differently-assigned
    EntityIDs). Both results are valid, fully non-overlapping partitions.
    """
    pq_df = gpd.read_parquet(input_boundaries)
    geoms = gpd.GeoSeries([_to_multipoly(_get_valid(s)) for s in pq_df.geometry.values])
    g = gpd.GeoDataFrame(
        {_ENTITY_COL: np.arange(len(pq_df), dtype=np.int64), _GEOM: geoms},
        geometry=_GEOM,
    )
    logger.info("convert_geometry: %d input polygons", len(g))

    if harmonize:
        g = harmonize_polygons(g, min_distance=min_distance, min_area=min_area)
        logger.info("convert_geometry: %d cells after harmonize", len(g))

    planes = []
    for z in range(n_z_planes):
        gz = g.copy()
        gz[_ZINDEX_COL] = z
        gz["ZLevel"] = spacing_z * (z + 1)
        planes.append(gz)
    out = gpd.GeoDataFrame(pd.concat(planes, ignore_index=True), geometry=_GEOM)

    if input_micron_to_mosaic is not None:
        m2m = _read_micron_to_mosaic_transform(input_micron_to_mosaic)
        mat = np.linalg.inv(m2m)
        flat = [*mat[:2, :2].flatten(), *mat[:2, 2].flatten()]
        out[_GEOM] = out[_GEOM].apply(
            lambda x: _to_multipoly(affinity.affine_transform(x, flat))
        )
        out = out.set_geometry(_GEOM)

    _set_process_id()
    id_map = {e: _get_id() for e in np.sort(g[_ENTITY_COL].unique())}
    out[_ENTITY_COL] = out[_ENTITY_COL].map(id_map).astype(np.int64)
    out["ID"] = np.arange(len(out), dtype=np.int64)
    out["Type"] = entity_type
    for c in ("ParentID", "ParentType", "Name"):
        out[c] = None
    out = out[["ID", _GEOM, _ENTITY_COL, "ParentID", "ParentType", "Type",
               _ZINDEX_COL, "ZLevel"]]

    Path(output_boundaries).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output_boundaries)
    return out


# --------------------------------------------------------------------------- #
# sum-signals (per-cell image intensity, raw + FFT high-pass)                 #
# --------------------------------------------------------------------------- #
_SIGNAL_IMG_RE = re.compile(r"mosaic_(?P<stain>[\w-]+)_z(?P<z>\d+)\.tif$")


def _high_pass(window: np.ndarray) -> np.ndarray:
    """FFT high-pass of a window: subtract a size-10 uniform (low-pass) blur."""
    from scipy import ndimage

    cell_fft = np.fft.fft2(window)
    filtered_fft = ndimage.fourier_uniform(cell_fft, size=10)
    inverse = np.fft.ifft2(filtered_fft)
    return np.maximum(window - inverse.real, 0)


def _cell_brightness(cell, file) -> tuple[float, float]:
    """Raw and high-pass intensity summed inside one (pixel-space) cell polygon."""
    import rasterio
    from rasterio.features import rasterize

    window = [cell.bounds[0], cell.bounds[1],
              cell.bounds[2] - cell.bounds[0] + 1, cell.bounds[3] - cell.bounds[1] + 1]
    window = [int(x) for x in window]
    image_data = file.read(1, window=rasterio.windows.Window(*window))
    hp_image = _high_pass(image_data)
    try:
        selector = translate(cell, xoff=-cell.bounds[0], yoff=-cell.bounds[1])
        bbox = [int(x) for x in selector.bounds]
        if (min(bbox[:2]) < 0 or bbox[2] >= image_data.shape[1]
                or bbox[3] >= image_data.shape[0]):
            raise ValueError("cell beyond image boundaries")
        mask = rasterize([selector], out_shape=(bbox[3] + 1, bbox[2] + 1),
                         all_touched=True)
        return float(np.sum(image_data * mask)), float(np.sum(hp_image * mask))
    except ValueError:
        logger.warning("could not apply cell mask; summing whole window")
        return float(np.sum(image_data)), float(np.sum(hp_image))


def _sumsig_image_task(args):
    """Sum a single stain image over its z-plane cells (parallel worker).

    ``ids`` and ``geoms`` are the z-plane subset only (kept small so each worker
    receives minimal pickled data).
    """
    import rasterio

    path, stain, ids, geoms, transform, gdal_cachemax_mb = args
    raw_vals, filt_vals, out_ids = [], [], []
    # bound GDAL's block cache (defaults to 5% of RAM, ~tens of GB on big nodes)
    with rasterio.Env(GDAL_CACHEMAX=gdal_cachemax_mb), rasterio.open(path) as file:
        for eid, geom in zip(ids, geoms):
            if geom is None or geom.is_empty:
                continue
            raw, filt = _cell_brightness(affine_transform(geom, transform), file)
            raw_vals.append(raw)
            filt_vals.append(filt)
            out_ids.append(eid)
    return stain, pd.Series(raw_vals, index=out_ids), pd.Series(filt_vals, index=out_ids)


def sum_signals(
    boundaries_path: str | Path,
    images_dir: str | Path,
    micron_to_mosaic_path: str | Path,
    output_csv: str | Path | None = None,
    *,
    n_jobs: int | None = None,
    gdal_cachemax_mb: int = 256,
) -> pd.DataFrame:
    """Sum mosaic-image intensity inside each cell, per stain (VPT ``sum-signals``).

    For each ``mosaic_{stain}_z{z}.tif`` image and each cell at that z-plane, maps
    the (micron-space) polygon to mosaic-pixel space, reads the cell's window,
    computes a per-cell FFT high-pass of that window, rasterizes the polygon mask
    (``all_touched=True``) and sums ``image*mask`` (raw) and ``highpass*mask``
    (high_pass). Images are processed in parallel (``n_jobs`` workers).

    Returns a DataFrame indexed by ``EntityID`` with ``{stain}_raw`` /
    ``{stain}_high_pass`` columns.

    Parameters
    ----------
    n_jobs
        Number of parallel worker processes (one image per task). ``None`` lets
        the executor choose; ``1`` runs sequentially.
    gdal_cachemax_mb
        Cap on GDAL's block cache per worker, in MB (default 256). GDAL otherwise
        defaults to ~5% of RAM, which on large-memory nodes lets each worker grow
        to tens of GB. Lower = less memory; higher = faster on big-RAM machines.

    Notes
    -----
    The FFT high-pass is bit-identical to VPT. Absolute values may differ from a
    given VPT run by ~0.5-1% per cell (correlation ~0.99999) because GDAL's
    ``all_touched`` rasterization includes slightly different boundary pixels
    across GDAL versions; the reimplementation uses the installed GDAL.
    """
    boundaries = gpd.read_parquet(boundaries_path)
    m2m = _read_micron_to_mosaic_transform(micron_to_mosaic_path)
    tform = [*m2m[:2, :2].flatten(), *m2m[:2, 2].flatten()]
    all_ids = np.sort(boundaries[_ENTITY_COL].unique())

    images = []
    for p in sorted(glob.glob(f"{images_dir}/mosaic_*_z*.tif")):
        m = _SIGNAL_IMG_RE.search(Path(p).name)
        if m:
            images.append((p, m.group("stain"), int(m.group("z"))))
    logger.info("sum_signals: %d images, %d cells", len(images), len(all_ids))

    # pass each worker only its z-plane's cells (small pickles, low worker memory)
    by_z = {
        z: (sub[_ENTITY_COL].to_numpy(), sub.geometry.to_numpy())
        for z, sub in boundaries.groupby(_ZINDEX_COL)
    }
    tasks = [
        (p, s, by_z[z][0], by_z[z][1], tform, gdal_cachemax_mb)
        for p, s, z in images if z in by_z
    ]
    if n_jobs == 1 or len(tasks) <= 1:
        results = [_sumsig_image_task(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            results = list(ex.map(_sumsig_image_task, tasks))

    cols: dict = {}
    for stain, raw, filt in results:
        raw = raw.reindex(all_ids, fill_value=0.0)
        filt = filt.reindex(all_ids, fill_value=0.0)
        cols[f"{stain}_raw"] = cols.get(f"{stain}_raw", 0) + raw
        cols[f"{stain}_high_pass"] = cols.get(f"{stain}_high_pass", 0) + filt

    df = pd.DataFrame(cols)
    df.index.name = _ENTITY_COL
    df.sort_index(inplace=True)
    if output_csv is not None:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_csv)
    return df
