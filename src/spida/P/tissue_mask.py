"""KDE-domain tissue mask for transcript QC (successor to the hex ``transcript_qc``).

Computes a smooth tissue footprint from transcript density and persists it three ways:
  1. a micron-space pass-region polygon (``ShapesModel``) consumed by the post-hoc
     segmentation filter (``apply_qc_filter_to_segmentation``);
  2. a native-grid boolean ``Labels2DModel`` in the SpatialData object (for overlay/QC);
  3. a pixel-aligned boolean ``.tif`` written into the raw images directory, so the
     deconvolution / segmentation steps can consume it without opening the zarr or the
     micron->mosaic transform.

Pipeline (per region), mirroring the validated dev prototype:
  bin all transcripts onto a fine micron grid -> density at a FIXED reference scale
  (bandwidth-independent) >= ``RHO_MIN`` -> binary_closing -> drop connected components
  < ``A_MIN_UM2`` -> erode boundary by ``ERODE_SIGMA`` reference sigmas.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import tifffile
from rasterio.features import rasterize, shapes as rio_shapes
from rasterio.transform import Affine as RioAffine
from scipy.ndimage import binary_closing, binary_erosion, gaussian_filter, label
from shapely.affinity import affine_transform

logger = logging.getLogger(__package__)

# --- Calibrated defaults (validated in spida_dev/v1/tissue_mask) ---
GB = 5.0            # micron grid spacing (<= REF_UM/2 for a faithful density estimate)
RHO_MIN = 5.0       # min transcript density (tx/100um^2) to count as tissue (inclusive)
REF_UM = 25.0       # FIXED reference scale for the density estimate (bandwidth-independent)
A_MIN_UM2 = 50000.0  # drop connected components smaller than this (~0.05 mm^2)
ERODE_SIGMA = 1.0   # boundary erosion in reference-sigma units; erosion_um = erode_sigma * REF_UM

QC_MASK_SHAPES_KEY = "qc_mask"
QC_MASK_LABELS_KEY = "tissue_mask"
QC_MASK_TIF_NAME = "tissue_mask.tif"


# --------------------------------------------------------------------------- compute

def _ensure_xy(points, x_col: str = "x", y_col: str = "y") -> tuple[np.ndarray, np.ndarray]:
    """Extract x/y arrays from a (dask) DataFrame of transcript points."""
    if hasattr(points, "compute"):
        points = points[[x_col, y_col]].compute()
    return np.asarray(points[x_col].values), np.asarray(points[y_col].values)


def bin_transcripts(x: np.ndarray, y: np.ndarray, gb: float = GB):
    """Bin points onto a micron grid; return (tot counts [ny,nx], x0, y0)."""
    x0, y0 = float(np.min(x)), float(np.min(y))
    bx = np.floor((x - x0) / gb).astype(np.int64)
    by = np.floor((y - y0) / gb).astype(np.int64)
    nx, ny = int(bx.max()) + 1, int(by.max()) + 1
    tot = np.zeros((ny, nx), dtype=np.float64)
    np.add.at(tot, (by, bx), 1.0)
    return tot, x0, y0


def tissue_mask(
    tot: np.ndarray,
    gb: float = GB,
    rho_min: float = RHO_MIN,
    ref_um: float = REF_UM,
    a_min_um2: float = A_MIN_UM2,
    erode_sigma: float = ERODE_SIGMA,
) -> dict:
    """KDE-domain tissue mask; returns dict of stages for inspection.

    ``erode_sigma`` is in multiplicative units of the reference sigma: the boundary is
    eroded by ``erode_sigma * ref_um`` microns, so it tracks the density-estimate edge
    creep regardless of ``ref_um``. Keys: ``dens100`` (density field, tx/100um^2),
    ``raw`` (floor mask), ``big`` (after speck drop), ``final`` (after erosion).
    """
    Tref = gaussian_filter(tot, sigma=max(0.5, ref_um / gb), mode="constant")
    dens100 = 100.0 * Tref / gb**2
    raw = dens100 >= rho_min
    closed = binary_closing(raw, iterations=1)

    lab, n = label(closed)
    if n == 0:
        return {"dens100": dens100, "raw": raw, "big": closed, "final": closed}
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep = np.where(sizes >= a_min_um2 / gb**2)[0]
    big = np.isin(lab, keep)

    if erode_sigma and erode_sigma > 0:
        final = binary_erosion(big, iterations=max(1, round(erode_sigma * ref_um / gb)))
    else:
        final = big
    return {"dens100": dens100, "raw": raw, "big": big, "final": final}


def mask_to_polygons(mask: np.ndarray, x0: float, y0: float, gb: float = GB) -> gpd.GeoDataFrame:
    """Polygonize the boolean micron-grid mask -> pass-region polygons (micron coords).

    Returns a GeoDataFrame with a ``filtered`` column set to ``False`` (pass), matching the
    convention the post-hoc filter expects (keep geometries where ``filtered == qc_pass_value``).
    """
    # rasterio Affine maps (col,row)->(x,y): x = x0 + gb*col, y = y0 + gb*row
    transform = RioAffine(gb, 0.0, x0, 0.0, gb, y0)
    geoms = [
        shapely.geometry.shape(geom)
        for geom, val in rio_shapes(mask.astype(np.uint8), mask=mask, transform=transform)
        if val == 1
    ]
    return gpd.GeoDataFrame({"filtered": [False] * len(geoms)}, geometry=geoms, crs=None)


def read_micron_to_mosaic(path: str | Path) -> np.ndarray:
    """Read the 3x3 micron->mosaic-pixel affine (whitespace-delimited, MERSCOPE/SPIDA format)."""
    with open(path) as fh:
        rows = [list(map(float, ln.split())) for ln in fh if ln.strip()]
    return np.asarray(rows, dtype=float)


def rasterize_pixel_mask(gdf_micron: gpd.GeoDataFrame, m2m: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Transform micron pass-polygons to mosaic-pixel space and rasterize to ``shape`` (rows, cols)."""
    sx, _, tx = m2m[0]
    _, sy, ty = m2m[1]
    params = [sx, 0.0, 0.0, sy, tx, ty]  # shapely: (a,b,d,e,xoff,yoff)
    geoms_px = [affine_transform(g, params) for g in gdf_micron.geometry]
    if not geoms_px:
        return np.zeros(shape, dtype=bool)
    arr = rasterize(((g, 1) for g in geoms_px), out_shape=shape, fill=0, dtype="uint8")
    return arr.astype(bool)


def _labels_transformation(x0: float, y0: float, gb: float):
    """SpatialData transformation mapping the native grid -> global (micron): x=x0+gb*col."""
    from spatialdata.transformations import Scale, Sequence, Translation
    return Sequence([Scale([gb, gb], axes=("x", "y")), Translation([x0, y0], axes=("x", "y"))])


def _mosaic_shape(images_dir: Path, channel: str = "DAPI") -> tuple[int, int]:
    """(rows, cols) of a raw mosaic image (excluding decon outputs), header-only."""
    cands = sorted(
        p for p in images_dir.glob(f"mosaic_{channel}_z*.tif")
        if ".decon." not in p.name and "_z_stack_" not in p.name
    )
    if not cands:
        cands = sorted(p for p in images_dir.glob("mosaic_*_z*.tif")
                       if ".decon." not in p.name and "_z_stack_" not in p.name)
    with tifffile.TiffFile(cands[0]) as f:
        return tuple(f.series[0].shape)


# --------------------------------------------------------------------------- persist

def _overwrite_element(sdata, key: str, element) -> None:
    """Set + write a SpatialData element, replacing any on-disk copy (transcript_qc pattern)."""
    try:
        if hasattr(sdata, "elements_paths_on_disk"):
            paths = sdata.elements_paths_on_disk()
            if any(p.endswith(key) for p in paths) and hasattr(sdata, "delete_element_from_disk"):
                sdata.delete_element_from_disk(key)
    except Exception:
        pass
    sdata[key] = element
    if hasattr(sdata, "write_element"):
        sdata.write_element(key)


def run_tissue_mask(
    sdata,
    points_key: str,
    images_dir: str | Path,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    gb: float = GB,
    rho_min: float = RHO_MIN,
    ref_um: float = REF_UM,
    a_min_um2: float = A_MIN_UM2,
    erode_sigma: float = ERODE_SIGMA,
    shapes_key: str = QC_MASK_SHAPES_KEY,
    labels_key: str = QC_MASK_LABELS_KEY,
    write_raster: bool = True,
    m2m_path: str | Path | None = None,
    persist: bool = True,
):
    """Compute the tissue mask and (optionally) persist polygon + labels + images-dir tif.

    Returns (stages, grid_polygons, x0, y0). ``images_dir`` holds the raw mosaics and the
    ``micron_to_mosaic_pixel_transform.csv``; the boolean ``tissue_mask.tif`` is written there.
    """
    images_dir = Path(images_dir)
    logger.info("run_tissue_mask started (points_key=%s)", points_key)
    x, y = _ensure_xy(sdata[points_key], x_col=x_col, y_col=y_col)
    logger.info("Read %d transcripts", len(x))

    tot, x0, y0 = bin_transcripts(x, y, gb=gb)
    stages = tissue_mask(tot, gb=gb, rho_min=rho_min, ref_um=ref_um,
                         a_min_um2=a_min_um2, erode_sigma=erode_sigma)
    final = stages["final"]
    logger.info("Tissue mask: grid=%s, tissue=%d px (%.1f%%)",
                final.shape, int(final.sum()), 100 * final.mean())

    grid = mask_to_polygons(final, x0, y0, gb=gb)
    logger.info("Polygonized pass region into %d polygon(s)", len(grid))

    if not persist:
        return stages, grid, x0, y0

    from spatialdata.models import Labels2DModel, ShapesModel

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        _overwrite_element(sdata, shapes_key, ShapesModel.parse(grid))
        labels = Labels2DModel.parse(
            final.astype(np.uint8), dims=("y", "x"),
            transformations={"global": _labels_transformation(x0, y0, gb)},
        )
        _overwrite_element(sdata, labels_key, labels)
    logger.info("Persisted shapes '%s' and labels '%s'", shapes_key, labels_key)

    if write_raster:
        if m2m_path is None:
            m2m_path = images_dir / "micron_to_mosaic_pixel_transform.csv"
        m2m = read_micron_to_mosaic(m2m_path)
        shape = _mosaic_shape(images_dir)
        px_mask = rasterize_pixel_mask(grid, m2m, shape)
        out_tif = images_dir / QC_MASK_TIF_NAME
        tifffile.imwrite(out_tif, (px_mask.astype(np.uint8) * 255), compression="deflate")
        logger.info("Wrote pixel mask %s (%s, %.1f MB, tissue=%.1f%%)",
                    out_tif, px_mask.shape, out_tif.stat().st_size / 1e6, 100 * px_mask.mean())

    logger.info("run_tissue_mask completed")
    return stages, grid, x0, y0
