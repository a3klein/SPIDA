"""Convert segmentation label masks (numpy arrays) to boundary polygons.

Shared by the mask-based backends (cellpose, mesmer). The raw output of these
segmenters is a big integer label array; this module vectorizes each labelled
region into a polygon (tiled + parallelized for large mosaics) and returns a
GeoDataFrame with columns ``ID`` (label), ``Geometry``, ``z`` (plane).

No cellpose/deepcell dependency — importable from any env.
"""

from __future__ import annotations

import logging
import itertools
import multiprocessing as mp
from functools import partial

import numpy as np
import geopandas as gpd
import skimage as ski
from shapely.geometry import Polygon
from tqdm import tqdm

from spida.utilities.tiling import tile_image_with_overlap, merge_tile_polygons

logger = logging.getLogger(__package__)

# default worker count; callers may override (e.g. tests set masks.max_cpu = 1)
max_cpu = mp.cpu_count()


def convert_mask(mask_id, masks, tolerance: float = 0.5):
    """Convert a single labelled region to polygon(s), one per occupied z-plane.

    Parameters
    ----------
    mask_id : int
        The label value to extract.
    masks : np.ndarray
        A (Z, Y, X) label array (2D input is expanded to a single z-plane).
    tolerance : float
        Max geometry displacement for simplification (higher -> fewer vertices).
    """
    polygons = []
    cells = []
    layers = []
    mask = masks == mask_id
    for m in np.where(mask.any(axis=(1, 2)))[0]:
        # pad so contours of edge-touching regions close
        padded_mask = np.pad(mask[m], 1, mode="constant")
        contours = ski.measure.find_contours(padded_mask, 0.5)
        contours = [contour - 1 for contour in contours]
        if len(contours[0]) > 3:
            polygons.append(
                Polygon(contours[0]).simplify(
                    tolerance=tolerance, preserve_topology=True
                )
            )
            cells.append(mask_id)
            layers.append(m)
    return polygons, cells, layers


def _masks_to_polygons(masks: np.ndarray, tolerance: float = 0.5) -> gpd.GeoDataFrame:
    """Convert every labelled region in a (tile of a) mask array to polygons."""
    if masks.ndim == 2:
        masks = np.expand_dims(masks, axis=0)
    masks = masks.swapaxes(-1, -2)
    mask_ids = np.unique(masks)[1:]  # exclude background (0)

    polygons = []
    for mask_id in mask_ids:
        polygons_batch, cells_batch, layers_batch = convert_mask(
            mask_id, masks, tolerance=tolerance
        )
        polygons.append((polygons_batch, cells_batch, layers_batch))

    all_polygons = list(
        itertools.chain.from_iterable(p for p, _, _ in polygons)
    )
    all_cells = list(
        itertools.chain.from_iterable(c for _, c, _ in polygons)
    )
    all_layers = list(
        itertools.chain.from_iterable(z for _, _, z in polygons)
    )

    return gpd.GeoDataFrame(
        {"ID": all_cells, "Geometry": all_polygons, "z": all_layers},
        geometry="Geometry",
        crs=None,
    )


def masks_to_geodataframe(
    masks: np.ndarray,
    tolerance: float = 0.5,
    tile_size: int | tuple = 1000,
    overlap: int | tuple = 200,
) -> gpd.GeoDataFrame:
    """Convert a full label mask (2D or 3D) to a polygon GeoDataFrame.

    Tiles the mask (with overlap), vectorizes each tile in parallel, and merges
    tile-local polygons back into global coordinates.
    """
    tiles, tile_info = tile_image_with_overlap(
        masks, tile_size=tile_size, overlap=overlap
    )
    logger.info(f"Generated {len(tiles)} tiles.")

    parallel_func = partial(_masks_to_polygons, tolerance=tolerance)
    with mp.Pool(max_cpu) as pool:
        geo_list = list(tqdm(pool.map(parallel_func, tiles), total=len(tiles)))

    merged_polygons = merge_tile_polygons(tile_info, geo_list, geom_col="Geometry")

    # Dissolve each cell's per-tile fragments into one polygon per (cell, z-plane).
    # Tiling is an internal implementation detail (for parallelism/memory), so the
    # returned polygons must not leak tile fragments. Dissolving by ["ID", "z"]
    # handles both cases: in 2D `z` is constant (== dissolve by ID); in 3D it keeps
    # z-planes separate (dissolving by ID alone would collapse a cell's planes).
    if not merged_polygons.empty and "ID" in merged_polygons.columns:
        merged_polygons = merged_polygons.dissolve(by=["ID", "z"]).reset_index()

    logger.info("Finished mask -> polygon conversion")
    return merged_polygons
