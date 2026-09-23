"""Image tile-filtering strategies for the deconvolution / segmentation tiling.

Selects which tiles are passed to the downstream algorithm and drops the rest. Strategies
are kept SEPARATE (never combined):

- ``otsu``        : keep a tile if its max intensity exceeds ``threshold_otsu(maxes)/2``
                    (the historical SPIDA decon default; intensity-based).
- ``tissue_mask`` : keep a tile if it overlaps any tissue pixel of a precomputed tissue
                    mask (transcript-based; added in Stage 3).
- ``none``        : keep all tiles.

Each strategy returns ``(kept_tiles, kept_tile_info, tile_maxes, thr)`` and, matching the
original inline behavior, records ``max_intensity`` (and ``thresholded`` for otsu) on every
entry of ``tile_info`` so the ``--plot_thr`` grid visualization keeps working.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import skimage as ski

logger = logging.getLogger(__package__)


def _record_tile_maxes(tiles, tile_info):
    """Set ``max_intensity`` on each tile_info entry and return the list of maxes."""
    for _t, _ti in zip(tiles, tile_info):
        _ti["max_intensity"] = np.max(_t)
    return [ti["max_intensity"] for ti in tile_info]


def otsu_tile_filter(tiles, tile_info):
    """Keep tiles with max intensity above ``threshold_otsu(all maxes) / 2``."""
    tile_maxes = _record_tile_maxes(tiles, tile_info)
    thr = ski.filters.threshold_otsu(np.asarray(tile_maxes)) / 2
    for _ti in tile_info:
        _ti["thresholded"] = _ti["max_intensity"] > thr
    ret = [(t, ti) for t, ti in zip(tiles, tile_info) if ti["thresholded"]]
    kept_tiles = [t for t, _ in ret]
    kept_info = [ti for _, ti in ret]
    logger.info("Otsu tile filter: thr=%s, kept %d/%d tiles", thr, len(kept_tiles), len(tiles))
    return kept_tiles, kept_info, tile_maxes, thr


def no_tile_filter(tiles, tile_info):
    """Keep all tiles (records maxes for plotting; thr=0)."""
    tile_maxes = _record_tile_maxes(tiles, tile_info)
    return list(tiles), list(tile_info), tile_maxes, 0


def _load_mask(mask):
    """Return a 2D mask array from a path (memmap; imread fallback) or pass an ndarray through."""
    if isinstance(mask, (str, Path)):
        import tifffile
        try:
            return tifffile.memmap(mask)
        except (ValueError, MemoryError):
            return tifffile.imread(mask)
    return mask


def _tile_rowcol(tile_info_entry) -> tuple[int, int, int, int]:
    """
    Returns ``(r0, c0, r1, c1)`` from a 2D ``(row, col)`` or 3D ``(z, row, col)`` position.

    ``tile_image_with_overlap`` branches on ``len(shape)``: the 2D path emits
    ``(row, col)`` positions, the 3D path emits ``(d, row, col)``. ``decon_image``
    passes the full ``(z, H, W)`` volume, so 3D deconvolution always produces
    3-element positions. 
    """
    pos, end = tile_info_entry["position"], tile_info_entry["end_position"]
    if len(pos) == 3:
        return pos[1], pos[2], end[1], end[2]
    if len(pos) == 2:
        return pos[0], pos[1], end[0], end[1]
    raise ValueError(
        f"tile position must have 2 (row, col) or 3 (z, row, col) elements, "
        f"got {len(pos)}: {pos!r}"
    )


def tissue_mask_tile_filter(tiles, tile_info, mask_tif):
    """Keep a tile if it overlaps any tissue pixel of the (pixel-aligned) tissue mask.

    ``mask_tif`` is a path to the boolean ``tissue_mask.tif`` (or an ndarray). Records
    ``max_intensity`` (for the plot_thr histogram) and ``thresholded`` (= overlaps tissue)
    on each tile_info entry; thr is not intensity-based here, so it is 0.
    """
    mask = _load_mask(mask_tif)
    tile_maxes = _record_tile_maxes(tiles, tile_info)  # for the plot_thr histogram/grid
    for _ti in tile_info:
        r0, c0, r1, c1 = _tile_rowcol(_ti)
        _ti["thresholded"] = bool(np.any(mask[r0:r1, c0:c1]))
    ret = [(t, ti) for t, ti in zip(tiles, tile_info) if ti["thresholded"]]
    kept_tiles = [t for t, _ in ret]
    kept_info = [ti for _, ti in ret]
    logger.info("Tissue-mask tile filter: kept %d/%d tiles", len(kept_tiles), len(tiles))
    return kept_tiles, kept_info, tile_maxes, 0


def select_tiles(tiles, tile_info, method: str = "otsu", **kwargs):
    """Dispatch to a tile-filter strategy.

    Returns ``(kept_tiles, kept_tile_info, tile_maxes, thr)``. ``method`` is one of
    ``{"otsu", "tissue_mask", "none"}`` (kept separate, never combined). ``tissue_mask``
    requires ``mask_tif=`` (path or ndarray). Extra kwargs pass to the chosen strategy.
    """
    if method == "otsu":
        return otsu_tile_filter(tiles, tile_info)
    if method in ("none", None):
        return no_tile_filter(tiles, tile_info)
    if method == "tissue_mask":
        mask_tif = kwargs.get("mask_tif")
        if mask_tif is None:
            raise ValueError("method='tissue_mask' requires mask_tif=<path or ndarray>")
        return tissue_mask_tile_filter(tiles, tile_info, mask_tif)
    raise ValueError(f"Unknown tile filter method: {method!r}")
