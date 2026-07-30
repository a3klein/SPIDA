"""Cellpose segmentation backend.

Emits RAW output only — ``polygons.parquet`` with ``ID`` (the global cellpose
label), ``z`` and ``Geometry`` in full-resolution pixel space, plus a
``meta.json`` recording exactly how the run was configured. The
segmentation-schema fields (``ZIndex``/``ZLevel``/``EntityID``) are derived
downstream by ``ingest_polygons``, so those conventions live in one place.

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

from natsort import natsorted
import skimage as ski
import cv2  # NOTE: must precede the cellpose import — see below

from .masks import masks_to_geodataframe
from .config import ConfigError
from .cellpose_config import (
    CELLPOSE_ENV_VAR,
    SPIDA_ENV_VAR,
    CellposeConfig,
    resolve_model,
)

# .env is read BEFORE cellpose is imported on purpose: cellpose resolves its
# model cache (``models.MODEL_DIR``) from CELLPOSE_LOCAL_MODELS_PATH at *import*
# time, so bridging SPIDA's own key across has to happen first to have any
# effect. The CLI does the same via ``config.export_env_only_config``; this
# covers programmatic callers that bypass the CLI.
load_dotenv()
if os.getenv(SPIDA_ENV_VAR) and not os.getenv(CELLPOSE_ENV_VAR):
    os.environ[CELLPOSE_ENV_VAR] = os.environ[SPIDA_ENV_VAR]

# The `import cv2` above must stay ABOVE this import: `cellpose.models` pulls in
# torch, which loads the system libstdc++; if that happens first, cv2 then fails
# with `version CXXABI_1.3.15 not found`.
from cellpose import models, core, io  # type: ignore  # noqa: E402

# io.logger_setup() # TODO: Maybe put this in a temporary file since this creates race conditions between multiple runs.

logger = logging.getLogger(__package__)

MODE_2D = "2d"
MODE_STITCHED = "3d_stitched"
MODE_TRUE_3D = "3d_true"

META_FILE = "meta.json"


def _require_gpu():
    if not core.use_gpu():
        raise ImportError("No GPU access")


def _load_image(
    image_path: Path,
    image_ext: str = ".tif",
    nuc_stain_name: str = "DAPI",
    cyto_stain_name: str | None = None,
    downscale: int | None = None,
):
    """
    Load an image from the specified path.

    Parameters:
    image_path (str): Path to the image file.
    image_ext (str): Extension of the image file (default is '.tif').
    nuc_stain_name (str): Name of the nuclear stain (default is 'DAPI').
    cyto_stain_name (str): Name of the cytoplasmic stain (default is 'PolyT').

    Returns:
    np.ndarray: Loaded image as a numpy array.
    """
    files = natsorted(
        [
            f
            for f in image_path.glob("*" + image_ext)
            if "_masks" not in f.name and "_flows" not in f.name
        ]
    )
    logger.info(f"Found {len(files)} images in {image_path}")

    nuc_files = []
    for f in files:
        if nuc_stain_name in f.name:
            nuc_files.append(f)

    if len(nuc_files) == 0:
        raise ValueError(
            f"No nuclear stain file found with name '{nuc_stain_name}' in {image_path}"
        )
    elif len(nuc_files) == 1:
        nuc_file = nuc_files[0]
        nuc_img = io.imread(nuc_file)
        if downscale and downscale > 1:
            factors = (1,) * (nuc_img.ndim - 2) + (downscale, downscale)
            nuc_img = ski.transform.downscale_local_mean(nuc_img, factors)
    else:
        # Stack 3D image
        nuc_files = natsorted(nuc_files)
        stack_3d = []
        for _fn in nuc_files:
            logger.info(f"  {_fn}")
            _img = io.imread(_fn)
            if downscale and downscale > 1:
                _img = ski.transform.downscale_local_mean(
                    _img, (downscale, downscale)
                )
            stack_3d.append(_img)
        nuc_img = np.stack(stack_3d, axis=-1)

    logger.info(f"Loaded images: Nuclear - {nuc_files}, Image Shape: {nuc_img.shape}")

    if cyto_stain_name is not None:
        cyto_files = []
        for f in files:
            if cyto_stain_name in f.name:
                cyto_files.append(f)

        if len(cyto_files) == 0:
            raise ValueError(
                f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}"
            )
        elif len(cyto_files) == 1:
            cyto_file = cyto_files[0]
            cyto_img = io.imread(cyto_file)
            if downscale and downscale > 1:
                factors = (1,) * (cyto_img.ndim - 2) + (downscale, downscale)
                cyto_img = ski.transform.downscale_local_mean(cyto_img, factors)
        else:
            # Stack 3D image
            cyto_files = natsorted(cyto_files)
            stack_3d = []
            for _fn in cyto_files:
                logger.info(f"  {_fn}")
                _img = io.imread(_fn)
                if downscale and downscale > 1:
                    _img = ski.transform.downscale_local_mean(
                        _img, (downscale, downscale)
                    )
                stack_3d.append(_img)
            cyto_img = np.stack(stack_3d, axis=-1)

        logger.info(f"Loaded images: Cytoplasmic - {cyto_files}, Image Shape: {cyto_img.shape}")
        return np.stack([nuc_img, cyto_img], axis=0)
    else:
        logger.info(
            f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}, returning nuclear image only."
        )
        return np.expand_dims(nuc_img, axis=0)


def _load_model(cfg: CellposeConfig):
    """Load the cellpose model named by ``cfg.model``.

    Name resolution and the existence check live in
    :func:`~.cellpose_config.resolve_model`, which never lets cellpose silently
    substitute a different model. Returns ``(model, resolution)`` so the caller
    can record exactly which weights were used.
    """
    _require_gpu()
    resolution = resolve_model(cfg.model, cfg.model_dir)
    resolution.log()
    model = models.CellposeModel(
        pretrained_model=resolution.pretrained_model,
        gpu=True,
        use_bfloat16=cfg.use_bfloat16,
    )
    return model, resolution


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


def build_eval_call(cfg: CellposeConfig, img: np.ndarray,
                    anisotropy: float | None = None):
    """Map ``cfg.seg_mode`` + input shape onto ``(image, model.eval kwargs)``.

    Pure and GPU-free so the mode semantics are unit-testable. ``img`` arrives
    channels-first — ``(C, Y, X)`` for one plane or ``(C, Y, X, Z)`` for a stack
    — and is transposed to cellpose's channels-last expectation.
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
        img = np.transpose(img, (3, 1, 2, 0))          # (C,Y,X,Z) -> (Z,Y,X,C)
        kwargs = dict(common, do_3D=True, z_axis=0, channel_axis=-1)
        if anisotropy is not None:
            kwargs["anisotropy"] = anisotropy
        return img, kwargs

    if cfg.seg_mode == MODE_STITCHED:
        if img.ndim == 4:
            img = np.transpose(img, (3, 1, 2, 0))      # (C,Y,X,Z) -> (Z,Y,X,C)
            return img, dict(common, do_3D=False, z_axis=0, channel_axis=-1,
                             stitch_threshold=cfg.stitch_threshold)
        # Degenerate: a single plane has nothing to stitch across. This is the
        # standard case for --force_2d_segmentation, where IMAGE_EXT selects one
        # z-slice file, so it must stay a warning rather than an error.
        logger.warning(
            "seg_mode='%s' but the input has a single z-plane; running plain 2D "
            "(stitch_threshold ignored).", MODE_STITCHED,
        )

    # 2D: explicitly requested, or a stitched request with single-plane input.
    if img.ndim == 3:
        img = np.moveaxis(img, 0, -1)                  # (C,Y,X) -> (Y,X,C)
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


def _cellpose_wrapper(img, cfg: CellposeConfig, anisotropy: float | None = None):
    """Load the model and run one ``model.eval``. Returns
    ``(masks, flows, styles, resolution)``."""
    model, resolution = _load_model(cfg)
    img, eval_kwargs = build_eval_call(cfg, img, anisotropy)

    logger.info("Running cellpose eval (seg_mode=%s, model=%s, image=%s):",
                cfg.seg_mode, cfg.model, img.shape)
    for key, value in eval_kwargs.items():
        logger.info("  %s: %s", key, value)

    masks, flows, styles = model.eval(img, **eval_kwargs)
    return masks, flows, styles, resolution


def _apply_clahe(img: np.ndarray, cfg: CellposeConfig) -> np.ndarray:
    """CLAHE contrast enhancement (single-plane input only, as before)."""
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


def _spida_version() -> str | None:
    try:
        from importlib.metadata import version
        return version("spida")
    except Exception:                                  # pragma: no cover
        return None


def run_cellpose(root_dir: str, output_dir: str, region: str, **kwargs):
    """Run cellpose segmentation for one region.

    Parameters:
    root_dir (str): Directory holding the region's mosaic stain TIFFs.
    output_dir (str): Fully-resolved output directory for this region+segmentation
        (the caller resolves the layout; nothing is appended here).
    region (str): Region name — used for logging and provenance only.
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
    for name, value in cfg.to_meta().items():
        logger.info("  %s = %s", name, value)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    timing: dict[str, float] = {}
    t0 = time.perf_counter()
    img = _load_image(
        Path(root_dir),
        image_ext=cfg.image_ext,
        nuc_stain_name=cfg.nuc_stain_name,
        cyto_stain_name=cfg.cyto_stain_name,
        downscale=cfg.scale,
    )
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

    anisotropy = cfg.resolved_anisotropy(root_dir)

    logger.info("STARTING SEGMENTATION (seg_mode=%s)", cfg.seg_mode)
    t0 = time.perf_counter()
    masks, _flows, _styles, resolution = _cellpose_wrapper(img, cfg, anisotropy)
    timing["segment_s"] = time.perf_counter() - t0
    logger.info("SEGMENTATION COMPLETED in %.1fs", timing["segment_s"])

    n_labels = int(len(np.unique(masks)) - 1)
    logger.info("Masks detected: %d; mask shape: %s", n_labels, masks.shape)

    # Key off the MASK, not the input image: seg_mode decides whether z survives.
    is_3d = masks.ndim == 3
    logger.info("Mask dimensionality: %dD; treating as %s for polygon conversion.",
                masks.ndim, "3D" if is_3d else "2D")

    # Inner (polygon-extraction) tiling. NOT the same as GPU tiling: this splits a
    # single label array whose IDs are already globally unique, so a cell crossing
    # a tile edge keeps its ID and its fragments re-merge by (ID, z) in the
    # dissolve inside masks_to_geodataframe.
    t0 = time.perf_counter()
    gdf = masks_to_geodataframe(
        masks,
        tolerance=cfg.poly_tolerance,
        tile_size=cfg.poly_tile_size,
        overlap=cfg.poly_tile_overlap,
    )
    n_raw = len(gdf)

    if gdf.empty:
        # A blank mask yields a frame carrying only tiling internals — no ID/z —
        # which would crash the min_z groupby and produce a parquet that
        # ingest_polygons cannot read. Normalise to a well-formed empty frame.
        logger.warning(
            "cellpose produced NO cells for region %s (mask labels: %d). Writing an "
            "empty polygons.parquet; check the input images and thresholds.",
            region, n_labels,
        )
        gdf = gpd.GeoDataFrame(
            {
                "ID": np.array([], dtype=np.int64),
                "z": np.array([], dtype=np.int64),
                "Geometry": gpd.GeoSeries([], dtype="geometry"),
            },
            geometry="Geometry",
            crs=None,
        )

    gdf.geometry = gdf.geometry.affine_transform([cfg.scale, 0, 0, cfg.scale, 0, 0])

    logger.info("Filtering out polygons with area <= %d (full-resolution px^2)",
                cfg.min_size)
    gdf = gdf[gdf.geometry.area > cfg.min_size]
    n_post_area = len(gdf)
    logger.info("  %d -> %d polygons", n_raw, n_post_area)

    if is_3d:
        logger.info("Filtering out cells spanning fewer than %d z-planes", cfg.min_z)
        # Group on ID (the globally-unique cellpose label), which is constant
        # across a cell's z-planes and across inner tiles.
        per_cell = gdf.groupby("ID")["z"].count()
        keep = per_cell[per_cell >= cfg.min_z].index
        gdf = gdf[gdf["ID"].isin(keep)].copy()
        logger.info("  %d -> %d polygons", n_post_area, len(gdf))
    n_final = len(gdf)
    timing["polygons_s"] = time.perf_counter() - t0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf.to_parquet(out / "polygons.parquet", index=False)

    meta = {
        "method": "cellpose",
        "region": region,
        "config": cfg.to_meta(),
        "model": {
            "name": resolution.name,
            "pretrained_model": resolution.pretrained_model,
            "model_dir": str(resolution.model_dir),
            "exists_locally": resolution.exists_locally,
            "downloaded": resolution.will_download,
            "size_bytes": resolution.size_bytes,
        },
        "anisotropy_resolved": anisotropy,
        "image": {"loaded_shape": loaded_shape, "eval_ndim": int(img.ndim)},
        "masks": {"shape": [int(n) for n in masks.shape], "n_labels": n_labels},
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
    (out / META_FILE).write_text(json.dumps(meta, indent=2, default=str))
    logger.info("Wrote %s and %s", out / "polygons.parquet", out / META_FILE)

    return is_3d
