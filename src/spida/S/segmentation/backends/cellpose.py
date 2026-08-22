"""
Cellpose segmentation backend.

Emits raw output - ``polygons.parquet`` with ``ID`` (global label),
``z`` and ``Geometry`` in full-resolution pixel space. Emits ``meta.json``, 
a recording of the runconfigurations.

All tunables are described by :class:`CellposeConfig`; see
``segment-region cellpose --list-params``.
"""

import os
import json
import time
from dotenv import load_dotenv  # type: ignore

import numpy as np
import geopandas as gpd
from pathlib import Path
import warnings
import logging

# Must come before importing cellpose.
# `cellpose.models` pulls in torch, which loads the system libstdc++
# importing cv2 after that raises `version CXXABI_1.3.15 not found`.
import cv2  

from ...io.mosaic import find_stain_planes, load_tissue_mask, probe_mosaic_shape, read_stains
from .masks import masks_to_geodataframe
from .config import ConfigError
from .cellpose_config import (
    CELLPOSE_ENV_VAR,
    SPIDA_ENV_VAR,
    CellposeConfig,
    resolve_model,
)

# DO NOT MOVE THIS BLOCK INTO A FUNCTION.
#
# cellpose resolves its model cache (``models.MODEL_DIR``) from
# CELLPOSE_LOCAL_MODELS_PATH at IMPORT time, so bridging SPIDA's own key across
# has to happen before the first `import cellpose.models`.
load_dotenv()
if os.getenv(SPIDA_ENV_VAR) and not os.getenv(CELLPOSE_ENV_VAR):
    os.environ[CELLPOSE_ENV_VAR] = os.environ[SPIDA_ENV_VAR]


logger = logging.getLogger(__package__)

MODE_2D = "2d"
MODE_STITCHED = "3d_stitched"
MODE_TRUE_3D = "3d_true"

META_FILE = "meta.json"

def _require_gpu():
    from cellpose import core
    if not core.use_gpu():
        raise ImportError("No GPU access")


def _resolve_stains(
    cfg: CellposeConfig,
    root_dir: str
) -> dict[str, list[Path]]:
    """
    Parameters:
        cfg (CellposeConfig): Configuration for cellpose segmentation.
        root_dir (str): Directory holding the region's mosaic stain TIFFs.

    Returns:
        dict[str, list[Path]]: ``{stain: z-ordered planes}`` for this run, nucleus FIRST.
    """
    images = Path(root_dir)
    stains = {
        cfg.nuc_stain_name: find_stain_planes(
            images, cfg.nuc_stain_name, cfg.image_ext)
    }
    if cfg.cyto_stain_name:
        stains[cfg.cyto_stain_name] = find_stain_planes(
            images, cfg.cyto_stain_name, cfg.image_ext)
    for stain, planes in stains.items():
        logger.info("%s: %d plane(s) - %s", stain, len(planes),
                    [f.name for f in planes])
    return stains


def _load_and_check_mask(
    tissue_mask_path,
    mosaic_shape: tuple[int, int]
) -> np.ndarray | None:
    """
    Load the tissue mask and assert it is on the mosaic's pixel grid.

    Parameters:
        tissue_mask_path (str | Path | None): Path to the tissue mask file.
        mosaic_shape (tuple[int, int]): Shape of the mosaic image.
    
    Returns:
        np.ndarray | None: The loaded tissue mask, or None if no path was provided.
    """
    if not tissue_mask_path:
        return None
    mask = load_tissue_mask(tissue_mask_path)
    if mask.shape != mosaic_shape:
        raise ConfigError(
            f"tissue mask {tissue_mask_path} has shape {mask.shape} but the "
            f"mosaic is {mosaic_shape}; they must be pixel-aligned."
        )
    logger.info("tissue mask: %.1f%% of the mosaic is tissue",
                100 * float(mask.mean()))
    return mask


# Load in cellpose model
def _load_model(cfg: CellposeConfig):
    """Load the cellpose model named by ``cfg.model``."""
    from cellpose import models

    _require_gpu()
    resolution = resolve_model(cfg.model, cfg.model_dir)
    resolution.log()
    model = models.CellposeModel(
        pretrained_model=resolution.pretrained_model,
        gpu=True,
        use_bfloat16=cfg.use_bfloat16,
    )
    return model, resolution

# Different options for collapzing multiple z-slices into a single plane for a 2D segmentation. 
def reduce_z(img: np.ndarray, z_reduce: str) -> np.ndarray:
    """Collapse a ``(C, Y, X, Z)`` stack to a single ``(C, Y, X)`` plane."""
    if img.ndim != 4:
        return img
    if z_reduce == "mip":
        logger.info("z_reduce='mip': maximum-intensity projection over %d planes",
                    img.shape[-1])
        return np.max(img, axis=-1)
    mid = img.shape[-1] // 2
    logger.info("z_reduce='middle': taking z-plane %d of %d", mid, img.shape[-1])
    return img[..., mid]


# GPU-free for unit-testing purposes. 
def build_eval_call(
    cfg: CellposeConfig,
    img: np.ndarray,
    anisotropy: float | None = None
):
    """
    Map cfg.seg_mode + input shape onto (image, model.eval kwargs).

    Parameters:
    cfg (CellposeConfig): Configuration for cellpose segmentation.
    img (np.ndarray): Input image array, either (C, Y, X) for single-plane or (C, Y, X, Z) for a stack.
    anisotropy (float | None): Anisotropy factor for 3D segmentation, if applicable.
    """
    common = dict(
        batch_size=cfg.batch_size,
        flow_threshold=cfg.flow_threshold,
        cellprob_threshold=cfg.cellprob_threshold,
        diameter=cfg.diameter,
        min_size=cfg.cp_min_size,
        normalize={"tile_norm_blocksize": cfg.tile_norm_blocksize},
    )

    if cfg.seg_mode == MODE_TRUE_3D:
        if img.ndim != 4:
            raise ConfigError(
                f"seg_mode='{MODE_TRUE_3D}' needs a multi-plane z-stack "
                f"(C, Y, X, Z); got shape {img.shape}. Use seg_mode='{MODE_2D}' "
                f"for single-plane input."
            )
        img = np.transpose(img, (3, 1, 2, 0)) # (C,Y,X,Z) -> (Z,Y,X,C)
        kwargs = dict(common, do_3D=True, z_axis=0, channel_axis=-1)
        if anisotropy is not None: # anisotropy is a physical property of the image. (XY vs. Z spacing)
            kwargs["anisotropy"] = anisotropy
        if cfg.flow3d_smooth_z > 0: # smoothing in 3D is relevant when anisotropy is applied to control for rings/holes.
            kwargs["flow3D_smooth"] = [float(cfg.flow3d_smooth_z), 0.0, 0.0]
        return img, kwargs

    if cfg.seg_mode == MODE_STITCHED:
        if img.ndim == 4:
            img = np.transpose(img, (3, 1, 2, 0)) # (C,Y,X,Z) -> (Z,Y,X,C)
            return img, dict(common, do_3D=False, z_axis=0, channel_axis=-1,
                             stitch_threshold=cfg.stitch_threshold)
        # Degenerate: a single plane has nothing to stitch across.
        logger.warning(
            "seg_mode='%s' but the input has a single z-plane; running plain 2D "
            "(stitch_threshold ignored).", MODE_STITCHED,
        )

    # 2D: explicitly requested, or a stitched request with single-plane input.
    if img.ndim == 3:
        img = np.moveaxis(img, 0, -1) # (C,Y,X) -> (Y,X,C)
        channel_axis = -1
    elif img.ndim == 2:
        channel_axis = None
    else:
        raise ConfigError(
            f"seg_mode='{cfg.seg_mode}' needs a single z-plane; got shape "
            f"{img.shape}. Set z_reduce to collapse the stack first."
        )
    return img, dict(common, do_3D=False, channel_axis=channel_axis,
                     stitch_threshold=0)


def _log_eval_call(
    what: str,
    cfg: CellposeConfig,
    shape,
    eval_kwargs: dict
) -> None:
    """
    Log one `model.eval` invocation's resolved kwargs.
    """
    logger.info("%s (seg_mode=%s, model=%s, image=%s):",
                what, cfg.seg_mode, cfg.model, shape)
    for key, value in eval_kwargs.items():
        logger.info("  %s: %s", key, value)


# Wraps the internal cellpose calls: loads model + config and runs the eval call.
def _cellpose_wrapper(
    img,
    cfg: CellposeConfig,
    anisotropy: float | None = None
):
    """
    Load the model and run one ``model.eval``. Returns ``(masks, flows, styles, resolution)``.

    Parameters:
    img (np.ndarray): Input image array, either (C, Y, X) for single-plane or (C, Y, X, Z) for a stack.
    cfg (CellposeConfig): Configuration for cellpose segmentation.
    anisotropy (float | None): Anisotropy factor for 3D segmentation, if applicable.
    """
    model, resolution = _load_model(cfg)
    img, eval_kwargs = build_eval_call(cfg, img, anisotropy)

    _log_eval_call("Running cellpose eval", cfg, img.shape, eval_kwargs)
    masks, flows, styles = model.eval(img, **eval_kwargs)
    return masks, flows, styles, resolution


# For optional image preprocessing just to cellpose
def _apply_clahe(img: np.ndarray, cfg: CellposeConfig) -> np.ndarray:
    """
    CLAHE contrast enhancement (single-plane input only, as before).

    Parameters:
    img (np.ndarray): Input image array.
    cfg (CellposeConfig): Configuration for cellpose segmentation.
    """
    logger.info("Applying CLAHE with clipLimit=%s and tileGridSize=%s",
                cfg.clahe_clip_limit, cfg.clahe_tile_grid_size)
    clahe = cv2.createCLAHE(
        clipLimit=cfg.clahe_clip_limit,
        tileGridSize=(cfg.clahe_tile_grid_size, cfg.clahe_tile_grid_size),
    )
    img = img.astype(np.uint16)
    for channel in range(img.shape[0]):
        img[channel] = clahe.apply(img[channel])
    return (img / 256).astype("uint8")


# for meta
def _spida_version() -> str | None:
    try:
        from importlib.metadata import version
        return version("spida")
    except Exception:                                  # pragma: no cover
        return None


# For the case where a tile returns no detected polygons. 
def _empty_polygons() -> gpd.GeoDataFrame:
    """A well-formed empty result."""
    return gpd.GeoDataFrame(
        {
            "ID": np.array([], dtype=np.int64),
            "z": np.array([], dtype=np.int64),
            "Geometry": gpd.GeoSeries([], dtype="geometry"),
        },
        geometry="Geometry",
        crs=None,
    )

# This is for serializing the tiling plan to JSON.
def _plan_meta(plan: dict) -> dict:
    """
    Parameters:
    plan (dict): The tiling plan dictionary.

    Returns:
    dict: A JSON-safe summary of the tiling plan.
    """
    summary = {k: v for k, v in plan.items()
               if k not in ("tiles", "to_segment", "skipped")}
    summary["segment_tiles"] = [
        {"gi": t.gi, "gj": t.gj, "y0": t.y0, "x0": t.x0,
         "height": t.height, "width": t.width}
        for t in plan["to_segment"]
    ]
    return summary


# The single-call path: one ``model.eval`` over the whole input, no tiling.
def _segment_whole(
    cfg: CellposeConfig,
    root_dir: str,
    region: str,
    tissue_mask_path,
    timing: dict
) -> tuple[gpd.GeoDataFrame, dict]:
    """
    Single ``model.eval`` over an entire image.

    Parameters:
    cfg (CellposeConfig): Configuration for cellpose segmentation.
    root_dir (str): Directory holding the region's mosaic stain TIFFs.
    region (str): Region name — used for logging and provenance only.
    tissue_mask_path (str | Path | None): Full-resolution boolean tissue mask.
        Out-of-tissue pixels are zeroed during the read, so cells cannot be
        called outside the tissue. Previously only the tiled path could do this
        and the untiled path warned the mask was ignored; both now share one
        reader, so the mask costs nothing to honour here.
    timing (dict): Dictionary to record timing information.

    Returns:
        tuple[gpd.GeoDataFrame, dict]: A tuple containing the GeoDataFrame of 
        detected polygons and a dictionary with segmentation information.
    """
    t0 = time.perf_counter()
    stains = _resolve_stains(cfg, root_dir)
    mosaic_shape = probe_mosaic_shape(next(iter(stains.values()))[0])
    mask = _load_and_check_mask(tissue_mask_path, mosaic_shape)

    img = read_stains(stains, downscale=cfg.scale, mask=mask)   # (C, Y, X, Z)
    del mask

    # Collapse a single-plane stain to (C, Y, X) for the 2D path.
    if img.shape[-1] == 1:
        img = img[..., 0]
    timing["load_s"] = time.perf_counter() - t0
    loaded_shape = tuple(int(n) for n in img.shape)
    logger.info("Loaded image shape (downscaled by %d): %s", cfg.scale, loaded_shape)

    # 2D mode needs a single plane; collapse the stack if we were given one.
    if cfg.seg_mode == MODE_2D:
        img = reduce_z(img, cfg.z_reduce)
        logger.info("Image shape after z-reduction: %s", img.shape)

    if cfg.apply_clahe:
        if img.ndim == 3:
            img = _apply_clahe(img, cfg)
        else:
            logger.warning(
                "apply_clahe=True is only supported for single-plane input "
                "(got shape %s); skipping CLAHE.", img.shape,
            )

    # read in anisotropy from path if auto.
    anisotropy = cfg.resolved_anisotropy(root_dir)

    logger.info("STARTING SEGMENTATION (seg_mode=%s)", cfg.seg_mode)
    t0 = time.perf_counter()
    masks, _flows, _styles, resolution = _cellpose_wrapper(img, cfg, anisotropy) # call wrapper - actual segmentation.
    timing["segment_s"] = time.perf_counter() - t0
    logger.info("SEGMENTATION COMPLETED in %.1fs", timing["segment_s"])

    n_labels = int(len(np.unique(masks)) - 1)
    logger.info("Masks detected: %d; mask shape: %s", n_labels, masks.shape)

    # Inner (polygon-extraction) tiling. Not the same as GPU tiling.
    t0 = time.perf_counter()
    gdf = masks_to_geodataframe(
        masks,
        tolerance=cfg.poly_tolerance,
        tile_size=cfg.poly_tile_size,
        overlap=cfg.poly_tile_overlap,
    )
    if gdf.empty:
        logger.warning(
            "cellpose produced NO cells for region %s (mask labels: %d). Writing an "
            "empty polygons.parquet; check the input images and thresholds.",
            region, n_labels,
        )
        gdf = _empty_polygons()
    gdf.geometry = gdf.geometry.affine_transform([cfg.scale, 0, 0, cfg.scale, 0, 0])
    timing["polygons_s"] = time.perf_counter() - t0

    info = {
        # Key off the mask, not the input image: seg_mode decides whether z survives.
        "is_3d": masks.ndim == 3,
        "resolution": resolution,
        "anisotropy": anisotropy,
        "image": {"loaded_shape": loaded_shape, "eval_ndim": int(img.ndim)},
        "masks": {"shape": [int(n) for n in masks.shape], "n_labels": n_labels},
    }
    return gdf, info

# Multiple calls to model.eval with tiling across GPUs.
def _segment_tiled(
    cfg: CellposeConfig,
    root_dir: str,
    region: str,
    tissue_mask_path,
    timing: dict
) -> tuple[gpd.GeoDataFrame, dict]:
    """
    Builds tiles and calls one independent ``model.eval`` per tile, stitched by IoU.

    Parameters: 
    cfg (CellposeConfig): Configuration for cellpose segmentation.
    root_dir (str): Directory holding the region's mosaic stain TIFFs.
    region (str): Region name — used for logging and provenance only.
    tissue_mask_path (str | Path | None): Full-resolution boolean tissue mask.
    timing (dict): Dictionary to record timing information.

    Returns: 
        tuple[gpd.GeoDataFrame, dict]: A tuple containing the GeoDataFrame of 
        detected polygons and a dictionary with segmentation information.
    """
    from ..tiling import (
        DEFAULT_TILE_MIN,
        optimize_grid,
        plan_grid,
        run_tiled,
    )

    stain_planes = _resolve_stains(cfg, root_dir)
    mosaic_shape = probe_mosaic_shape(next(iter(stain_planes.values()))[0])
    n_planes = len(next(iter(stain_planes.values())))
    logger.info("tiled: mosaic %s, %d channel(s) x %d z-plane(s)",
                mosaic_shape, len(stain_planes), n_planes)

    t0 = time.perf_counter()
    mask = _load_and_check_mask(tissue_mask_path, mosaic_shape)
    timing["mask_s"] = time.perf_counter() - t0

    # Grid search block.
    tile = cfg.resolved_gpu_tile()
    if tile == "auto":
        if mask is None:
            # Fatal by choice
            raise ConfigError(
                "gpu_tile_size='auto' optimises the grid against the tissue mask, "
                "but no mask is available for this region. Run `tissue-mask` "
                "first, pass --tissue_mask_path, or give an explicit "
                "--gpu_tile_size."
            )
        # The search floor must respect min_gpu_tile()
        best, _scored = optimize_grid(
            mask, overlap=cfg.gpu_tile_overlap,
            tile_min=max(DEFAULT_TILE_MIN, cfg.min_gpu_tile()),
        )
        tile = (best["tile_y"], best["tile_x"])
        logger.info("gpu_tile_size='auto' -> %dx%d", *tile)

    plan = plan_grid(mosaic_shape, tile, cfg.gpu_tile_overlap, mask=mask)
    # build_grid clamps the tile to the tissue bbox, so a genuinely narrow
    # region can land under the floor through no fault of the config.
    if min(plan["tile_size"] if isinstance(plan["tile_size"], tuple)
           else (plan["tile_size"],)) < cfg.min_gpu_tile():
        logger.warning(
            "tiling: the tissue bbox %s clamped the tile to %s, below the "
            "%d px floor (poly_tile_size x scale). Segmentation is still "
            "correct.",
            plan["bbox"], plan["tile_size"], cfg.min_gpu_tile(),
        )

    # 1 model for the whole run, shared by every GPU thread.
    model, resolution = _load_model(cfg)
    anisotropy = cfg.resolved_anisotropy(root_dir)
    seen: dict = {}

    # The tiled counterpart of `_cellpose_wrapper`. Enables sharing of model across threads.
    def eval_tile(img: np.ndarray) -> np.ndarray:
        """``(C, Y, X, Z)`` -> label array. Called concurrently; must stay pure."""
        if cfg.seg_mode == MODE_2D:
            img = reduce_z(img, cfg.z_reduce)
        arr, eval_kwargs = build_eval_call(cfg, img, anisotropy)
        if "logged" not in seen:
            seen["logged"] = True
            _log_eval_call("tiled eval", cfg, arr.shape, eval_kwargs)
            seen["eval_shape"] = tuple(int(n) for n in arr.shape)
        masks, _flows, _styles = model.eval(arr, **eval_kwargs)
        seen.setdefault("mask_ndim", int(np.ndim(masks)))
        return masks

    logger.info("STARTING TILED SEGMENTATION (%d tiles, seg_mode=%s, %d GPU thread(s))",
                plan["n_segment"], cfg.seg_mode, cfg.n_gpu_threads)
    t0 = time.perf_counter()
    gdf = run_tiled( # Actually calling the segmentation.
        plan, stain_planes, eval_tile,
        scale=cfg.scale,
        mask=mask,
        poly_tile_size=cfg.poly_tile_size,
        poly_tile_overlap=cfg.poly_tile_overlap,
        poly_tolerance=cfg.poly_tolerance,
        n_gpu_threads=cfg.n_gpu_threads,
        n_prep=cfg.n_prep_workers,
        n_post=cfg.n_post_workers,
        iou_threshold=cfg.gpu_tile_iou,
        containment_threshold=cfg.gpu_tile_containment,
        merge=cfg.gpu_merge,
    )
    timing["segment_s"] = time.perf_counter() - t0
    logger.info("TILED SEGMENTATION COMPLETED in %.1fs", timing["segment_s"])

    # The rest of this function is about setting the appropriate metadata. 
    n_labels = int(gdf["ID"].nunique()) if len(gdf) else 0
    if gdf.empty:
        logger.warning(
            "cellpose produced NO cells for region %s across %d tiles. Writing an "
            "empty polygons.parquet; check the input images and thresholds.",
            region, plan["n_segment"],
        )
        gdf = _empty_polygons()

    info = {
        "is_3d": seen.get("mask_ndim", 2) == 3,
        "resolution": resolution,
        "anisotropy": anisotropy,
        "image": {
            "mosaic_shape": [int(n) for n in mosaic_shape],
            "n_z_planes": n_planes,
            "channels": list(stain_planes),
            "tile_eval_shape": seen.get("eval_shape"),
        },
        "masks": {"shape": None, "n_labels": n_labels},
        "tiling": _plan_meta(plan) | {
            "tissue_mask_path": str(tissue_mask_path) if tissue_mask_path else None,
            "iou_threshold": cfg.gpu_tile_iou,
            "containment_threshold": cfg.gpu_tile_containment,
            "merge": cfg.gpu_merge,
            "n_gpu_threads": cfg.n_gpu_threads,
            "n_prep_workers": cfg.n_prep_workers,
            "n_post_workers": cfg.n_post_workers,
        },
    }
    return gdf, info


def run_cellpose(root_dir: str, output_dir: str, region: str, *,
                 tissue_mask_path: str | Path | None = None, **kwargs):
    """Run cellpose segmentation for one region.

    Two segmentation paths, selected by ``gpu_tile_size``: the default ``None``
    runs a single ``model.eval`` over the whole input, while any other value
    tiles the mosaic, segments each tile independently through one shared model
    and stitches the polygons back together by IoU. Both return polygons in
    full-resolution pixel space, so everything downstream of the model eval is identical.

    Parameters:
    root_dir (str): Directory holding the region's mosaic stain TIFFs.
    output_dir (str): Fully-resolved output directory for this region+segmentation
        (the caller resolves the layout; nothing is appended here).
    region (str): Region name — used for logging and provenance only.
    tissue_mask_path (str | Path | None): Full-resolution boolean tissue mask.
        Used three ways on the tiled path — it sets the grid extent, drops tiles
        with no tissue, and zeroes out-of-tissue pixels before downscaling —
        and is resolved by the caller (``segment_region``), not here. Ignored
        with a warning on the untiled path, which segments its whole input.
    **kwargs: Any :class:`CellposeConfig` field. Unknown names raise, and
        deprecated ones (``do_3D``, ``project_3d_to_2d``, ``tile_size``,
        ``overlap``, ``model_name``, ...) are translated with a warning. Run
        ``segment-region cellpose --list-params`` for the full surface.

    Returns:
        bool: True if the produced masks are a 3D stack, False if a single plane.
    """
    cfg = CellposeConfig.from_kwargs(**kwargs)

    logger.info("run_cellpose: region=%s", region)
    logger.info("  root_dir   = %s", root_dir)
    logger.info("  output_dir = %s", output_dir)
    logger.info("  tissue_mask_path = %s", tissue_mask_path)
    for name, value in cfg.to_meta().items():
        logger.info("  %s = %s", name, value)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    timing: dict[str, float] = {}
    tiled = cfg.resolved_gpu_tile() is not None
    if tiled:
        gdf, info = _segment_tiled(cfg, root_dir, region, tissue_mask_path, timing)
    else:
        gdf, info = _segment_whole(cfg, root_dir, region, tissue_mask_path, timing)
    is_3d = info["is_3d"]

    t0 = time.perf_counter()
    n_raw = len(gdf)
    logger.info("Filtering out polygons with area <= %d (full-resolution px^2)",
                cfg.min_size)
    gdf = gdf[gdf.geometry.area > cfg.min_size]
    n_post_area = len(gdf)
    logger.info("  %d -> %d polygons", n_raw, n_post_area)

    if is_3d:
        logger.info("Filtering out cells spanning fewer than %d z-planes", cfg.min_z)
        per_cell = gdf.groupby("ID")["z"].count()
        keep = per_cell[per_cell >= cfg.min_z].index
        gdf = gdf[gdf["ID"].isin(keep)].copy()
        logger.info("  %d -> %d polygons", n_post_area, len(gdf))
    n_final = len(gdf)
    timing["filter_s"] = time.perf_counter() - t0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf.to_parquet(out / "polygons.parquet", index=False)

    resolution = info["resolution"]
    meta = {
        "method": "cellpose",
        "region": region,
        "path": "tiled" if tiled else "single",
        "config": cfg.to_meta(),
        "model": {
            "name": resolution.name,
            "pretrained_model": resolution.pretrained_model,
            "model_dir": str(resolution.model_dir),
            "exists_locally": resolution.exists_locally,
            "downloaded": resolution.will_download,
            "size_bytes": resolution.size_bytes,
        },
        "anisotropy_resolved": info["anisotropy"],
        "tissue_mask_path": str(tissue_mask_path) if tissue_mask_path else None,
        "image": info["image"],
        "masks": info["masks"],
        "polygon_counts": {
            "post_dissolve": n_raw,
            "post_area": n_post_area,
            "final": n_final,
        },
        "timing_s": {k: round(v, 3) for k, v in timing.items()},
        "versions": {
            "cellpose": getattr(__import__("cellpose"), "version", None),
            "spida": _spida_version(),
        },
    }
    if "tiling" in info:
        meta["tiling"] = info["tiling"]
    (out / META_FILE).write_text(json.dumps(meta, indent=2, default=str))
    logger.info("Wrote %s and %s", out / "polygons.parquet", out / META_FILE)

    return is_3d
