"""Cellpose run configuration + model resolution.

* :class:`CellposeConfig` — the full parameter surface of the cellpose backend.
* :func:`resolve_model` — turn a model *name* into something safe to hand to
  ``CellposeModel(pretrained_model=...)``, never silently substituting a
  different model.

Model resolution policy
-----------------------
Cellpose's own fallback is dangerous for us: given a path that does not exist
it rewrites the request to ``cpsam`` and merely logs a warning, so asking for
``cpsam_v2`` against a missing file yields *cpsam results labelled cpsam_v2*.
:func:`resolve_model` closes that hole by deciding up front:

1. name not known to SPIDA            -> raise
2. weights present in the model dir   -> pass the absolute path (exact load)
3. weights absent, cellpose can fetch -> pass the bare name, log the download
4. weights absent, cellpose cannot    -> raise, with the upgrade hint
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Literal

import numpy as np

from spida.io.formats.merscope import TRANSFORM_FILE
from .config import BackendConfig, ConfigError

logger = logging.getLogger(__package__)

#: Environment variable cellpose itself honours for its model cache. cellpose reads it
#: once, when ``cellpose.models`` is first imported (``models.MODEL_DIR``, not expanded).
CELLPOSE_ENV_VAR = "CELLPOSE_LOCAL_MODELS_PATH"

#: Models SPIDA knows about, with the notes that should inform picking one.
KNOWN_MODELS: dict[str, str] = {
    "cpsam": "Cellpose-SAM (ViT-L).",
    "cpsam_v2": "Improved SAM.",
    "cpdino": "DINOv3 ViT-L backbone. Needs the `dinov3` package.",
    "cpdino-vitb": "DINOv3 ViT-B backbone. ~2.2x faster and lowest VRAM. Needs `dinov3`.",
}

#: Models whose weights cannot be loaded without the DINOv3 package.
REQUIRES_DINOV3: frozenset[str] = frozenset({"cpdino", "cpdino-vitb"})


def resolve_model_dir(explicit: str | Path | None = None) -> Path:
    """Resolve the directory holding cellpose weights.

    Precedence: ``explicit`` (``--model_dir``, or ``CELLPOSE_MODEL_PATH`` from the SPIDA config,
    filled in by ``segment-region``) > ``CELLPOSE_LOCAL_MODELS_PATH`` (cellpose's own) >
    cellpose's default ``~/.cellpose/models``. Always ``~``-expanded.
    """
    for candidate in (explicit, os.getenv(CELLPOSE_ENV_VAR)):
        if candidate:
            return Path(candidate).expanduser()
    return Path.home() / ".cellpose" / "models"


def point_cellpose_at(directory: str | Path) -> None:
    """Make cellpose download and cache weights in ``directory`` (the one SPIDA reads).

    Sets ``CELLPOSE_LOCAL_MODELS_PATH`` for an import that has not happened yet, and
    updates ``cellpose.models.MODEL_DIR`` if cellpose was already imported in this process
    (its value is fixed at import time, and ``cache_model_path`` reads it on each call).
    """
    import sys

    directory = Path(directory).expanduser()
    os.environ[CELLPOSE_ENV_VAR] = str(directory)
    loaded = sys.modules.get("cellpose.models")
    if loaded is not None:
        loaded.MODEL_DIR = directory


def _cellpose_known_models() -> set[str]:
    """Names the installed cellpose can resolve and download on its own."""
    try:
        from cellpose import models as _cp_models  # noqa: PLC0415
    except Exception:                                # pragma: no cover
        return set()
    return set(getattr(_cp_models, "MODEL_NAMES", []) or [])


@dataclass(frozen=True)
class ModelResolution:
    """Outcome of :func:`resolve_model`."""

    name: str
    #: Value to pass as ``pretrained_model``: an absolute path when the weights
    #: are on disk, otherwise the bare name so cellpose fetches them.
    pretrained_model: str
    model_dir: Path
    exists_locally: bool
    will_download: bool
    size_bytes: int | None = None

    def log(self) -> None:
        if self.exists_locally:
            size_gb = (self.size_bytes or 0) / 1e9
            logger.info("cellpose model '%s': loading %s (%.2f GB)",
                        self.name, self.pretrained_model, size_gb)
        else:
            logger.warning(
                "cellpose model '%s' NOT FOUND in %s — cellpose will DOWNLOAD it "
                "there now. Point CELLPOSE_MODEL_PATH (config) or --model_dir at a "
                "shared directory to avoid re-downloading on every job.",
                self.name, self.model_dir,
            )


def resolve_model(model: str, model_dir: str | Path | None = None) -> ModelResolution:
    """Resolve a model *name* to a safe ``pretrained_model`` value.

    Never falls back to a different model: either the requested weights are
    used, or an exception is raised. See the module docstring for the policy.
    """
    if model not in KNOWN_MODELS:
        raise ConfigError(
            f"unknown cellpose model {model!r}. Known models: "
            f"{', '.join(sorted(KNOWN_MODELS))}."
        )

    directory = resolve_model_dir(model_dir)
    path = directory / model
    if path.is_file():
        return ModelResolution(
            name=model, pretrained_model=str(path), model_dir=directory,
            exists_locally=True, will_download=False,
            size_bytes=path.stat().st_size,
        )

    # Not on disk. Only let cellpose fetch it if this cellpose version actually
    # knows the name — otherwise it would quietly substitute cpsam.
    known = _cellpose_known_models()
    if known and model not in known:
        raise ConfigError(
            f"cellpose model {model!r} is not present in {directory} and the "
            f"installed cellpose does not know how to fetch it "
            f"(it knows: {', '.join(sorted(known)) or 'nothing'}). "
            f"Either place the weights in {directory} or upgrade the cellpose "
            f"pixi env — cpsam_v2/cpdino/cpdino-vitb need cellpose >= 4.2.1.1."
        )
    return ModelResolution(
        name=model, pretrained_model=model, model_dir=directory,
        exists_locally=False, will_download=True,
    )


#: Conventional filename of the micron->mosaic transform inside a region's images dir.
M2M_FILENAME = TRANSFORM_FILE


def read_native_xy_um(path: str | Path) -> float:
    """Micron size of one full-resolution mosaic pixel.

    Reads the 3x3 micron->pixel transform, whose element ``[0][0]`` is pixels per
    micron in x. ``path`` may be either the CSV file itself or a directory
    containing it (the region's images dir), so callers can supply an explicit
    location or rely on the conventional one.
    """
    path = Path(path).expanduser()
    if path.is_dir():
        path = path / M2M_FILENAME
    if not path.is_file():
        raise ConfigError(
            f"cannot derive anisotropy: {path} not found. Pass "
            f"--micron_to_mosaic_path=<file> to point at it explicitly, or give "
            f"an explicit --anisotropy=<float> to skip the derivation."
        )
    rows = [list(map(float, line.split())) for line in
            path.read_text().strip().splitlines() if line.strip()]
    matrix = np.array(rows)
    if matrix.shape != (3, 3):
        raise ConfigError(f"{path}: expected a 3x3 transform, got {matrix.shape}.")
    px_per_um = float(matrix[0][0])
    if px_per_um <= 0:
        raise ConfigError(f"{path}: non-positive pixels-per-micron ({px_per_um}).")
    return 1.0 / px_per_um


#: ``<tile_y>x<tile_x>`` — the non-square form of ``gpu_tile_size``.
_GPU_TILE_RE = re.compile(r"^\s*(\d+)\s*[x×]\s*(\d+)\s*$", re.IGNORECASE)

def parse_gpu_tile_size(value: int | str | None):
    """Normalise ``gpu_tile_size`` to ``None`` | ``"auto"`` | ``(tile_y, tile_x)``.

    Accepts an int (square tile), ``"<y>x<x>"`` (non-square) or ``"auto"``
    (optimise against the tissue mask). ``None`` means the untiled single-call
    path — today's behaviour.
    """
    if value is None or value == "auto":
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return (value, value)
    text = str(value).strip()
    if text.isdigit():
        # Programmatic callers bypass BackendConfig's str->int coercion.
        return (int(text), int(text))
    match = _GPU_TILE_RE.match(text)
    if match is None:
        raise ConfigError(
            f"gpu_tile_size={value!r} is not understood. Use an int for a square "
            f"tile (e.g. 9000), '<y>x<x>' for a non-square one (e.g. "
            f"'12000x9500'), or 'auto' to optimise the grid against the tissue "
            f"mask. Omit it entirely for the untiled single-call path."
        )
    return (int(match.group(1)), int(match.group(2)))


def derive_anisotropy(micron_per_z: float, native_xy_um: float, scale: int) -> float:
    """z:xy voxel-spacing ratio *as cellpose sees it*, i.e. after xy downsampling.

    ``anisotropy = z_step / (native_xy * scale)``. This is dictated by the data
    geometry, not a free parameter, and is coupled to ``scale``: change one and
    the other must change with it.
    """
    if scale < 1:
        raise ConfigError(f"scale must be >= 1, got {scale}.")
    return micron_per_z / (native_xy_um * scale)


@dataclass(frozen=True)
class CellposeConfig(BackendConfig):
    """Every parameter the cellpose backend accepts.

    Defaults reproduce current SPIDA production behaviour, with one deliberate
    exception: ``use_bfloat16`` now defaults to ``True`` (measured ~1.85x
    speedup at mean per-cell IoU 0.985 vs fp32).
    """

    METHOD: ClassVar[str] = "cellpose"

    # --- mode -------------------------------------------------------------
    seg_mode: Literal["2d", "3d_stitched", "3d_true"] = "3d_stitched"
    z_reduce: Literal["middle", "mip"] = "middle"

    # --- model ------------------------------------------------------------
    model: str = "cpsam"
    model_dir: str | None = None
    use_bfloat16: bool = True

    # --- geometry ---------------------------------------------------------
    scale: int = 4
    anisotropy: float | Literal["auto"] = "auto"
    micron_to_mosaic_path: str | None = None
    micron_per_z: float = 1.5
    n_z_planes: int = 7

    # --- cellpose eval ----------------------------------------------------
    flow_threshold: float = 0.0
    cellprob_threshold: float = -2.0
    tile_norm_blocksize: int = 0
    diameter: float | None = None
    cp_min_size: int = 100
    batch_size: int = 32
    stitch_threshold: float = 0.2
    flow3d_smooth_z: float = 1.0
    max_size_fraction: float = 1.0

    # --- inputs -----------------------------------------------------------
    nuc_stain_name: str = "DAPI"
    cyto_stain_name: str | None = None
    image_ext: str = ".tif"

    # --- optional preprocessing -------------------------------------------
    apply_clahe: bool = False
    clahe_clip_limit: float = 20.0
    clahe_tile_grid_size: int = 8

    # --- GPU tiling (whole-mosaic segmentation) ---------------------------
    gpu_tile_size: int | str | None = None
    gpu_tile_overlap: int = 500
    gpu_tile_iou: float = 0.15
    gpu_tile_containment: float = 0.5
    gpu_merge: Literal["union", "largest"] = "union"
    n_gpu_threads: int = 2
    n_prep_workers: int | None = None
    n_post_workers: int | None = None

    # --- polygon post-processing ------------------------------------------
    min_size: int = 100
    min_z: int = 3
    poly_tile_size: int = 1000
    poly_tile_overlap: int = 200
    poly_tolerance: float = 0.5
    n_poly_workers: int | None = None

    ALIASES: ClassVar[dict[str, str]] = {}  # deprecated names archived 2026-09-28

    CHOICES: ClassVar[dict[str, tuple]] = {
        "seg_mode": ("2d", "3d_stitched", "3d_true"),
        "z_reduce": ("middle", "mip"),
        "model": tuple(KNOWN_MODELS),
        "gpu_merge": ("union", "largest"),
    }

    HELP: ClassVar[dict[str, str]] = {
        "seg_mode": "2d = one z-slice, emitted as a single-plane mask; "
                    "3d_stitched = 2D per slice linked by stitch_threshold; "
                    "3d_true = joint 3D flows using anisotropy.",
        "z_reduce": "For seg_mode=2d only: take the middle z-slice, or a "
                    "maximum-intensity projection over z.",
        "model": "Cellpose model NAME (not a path). "
                 + " | ".join(f"{k}: {v}" for k, v in KNOWN_MODELS.items()),
        "model_dir": "Directory holding the weights (downloads land here too). Defaults to "
                     "CELLPOSE_MODEL_PATH in the config, then "
                     f"${CELLPOSE_ENV_VAR}, then ~/.cellpose/models.",
        "use_bfloat16": "Half-precision inference. ~1.85x faster at mean "
                        "per-cell IoU 0.985 vs fp32.",
        "scale": "xy downsample factor applied before inference. Coupled to anisotropy.",
        "anisotropy": "seg_mode=3d_true only. 'auto' derives z_step / "
                      "(native_xy * scale) from the region's "
                      "micron_to_mosaic transform, which is the physically correct value.",
        "micron_to_mosaic_path": (
            f"Explicit path to the micron-to-mosaic transform used to derive anisotropy. "
            f"Default None auto-discovers <images_dir>/{M2M_FILENAME}; set this only when the transform "
            f"lives elsewhere. Accepts a file or a containing directory."
        ),
        "micron_per_z": "Micron thickness of one z-step.",
        "n_z_planes": "Number of z-planes in the input stack.",
        "flow_threshold": "Cellpose's intenral flow threshold.",
        "cellprob_threshold": "Cellpose's internal cell probability threshold.",
        "tile_norm_blocksize": "Cellpose normalisation block size (0 = global).",
        "diameter": "Expected cell diameter in pixels; Cellpose internal parameter.",
        "cp_min_size": "Cellpose's internal minimum mask size, in downsampled pixels. Distinct from --min_size.",
        "batch_size": "Internal patches cellpose pushes through the network per forward pass.",
        "stitch_threshold": "seg_mode=3d_stitched only: IoU for linking masks across adjacent z-slices.",
        "flow3d_smooth_z": "seg_mode=3d_true only: Gaussian sigma applied to the 3D flow field along Z before mask reconstruction"
            "(cellpose `flow3D_smooth`). XY smoothing is deliberately not exposed.",
        "max_size_fraction": "Cellpose drops masks larger than this fraction of the image; its own default "
            "0.4 silently deletes large masks, so SPIDA uses 1.0 (keep all).",
        "cyto_stain_name": "Second channel (e.g. PolyT). Cheap to add (~1% GPU) "
                           "but changes the target from nucleus to whole cell.",
        "image_ext": "File extension of the input images (i.e. `.tif` or `.decon.tif`).",
        "apply_clahe": "Contrast-limited adaptive histogram equalisation before "
                       "inference.",
        "gpu_tile_size": (
            "Outer (GPU) tile edge in full-resolution pixels. None (default) "
            "runs the whole region in one model.eval call. "
            "An int gives a square tile; '<y>x<x>' a non-square one; "
            "'auto' searches the (tile_y, tile_x) surface for the grid that "
            "sends the least area to the GPU, which needs a tissue mask. Each "
            "tile is an independent model.eval, so cells crossing a tile edge "
            "are stitched together."
        ),
        "gpu_tile_overlap": "Outer-tile overlap in full-resolution pixels. Too small overlap might miss cells.",
        "gpu_tile_iou": "IoU for declaring two tiles' fragments the same cell. ",
        "gpu_tile_containment": "Intersection as a fraction of the smaller fragment's area to stitch cells. 0 disables this.",
        "gpu_merge": "How to resolve a merged cell's overlapping fragments on a shared z-plane."
                    " `union` dissolves them, `largest` keeps the biggest.",
        "n_gpu_threads": "Concurrent model.eval callers sharing 1 loaded model. More helpful for multi-GPU systems than 1 single GPU.",
        "n_prep_workers": "Threads reading and downscaling tile windows ahead of the GPU.",
        "n_post_workers": "Threads polygonising finished tiles behind the GPU. " 
                    "n_gpu_threads, n_prep_workers, and n_post_workers all compete for the job's CPUs.",
        "min_size": "Minimum polygon area, in full-resolution pixels squared, applied after upscaling. Distinct from --cp_min_size.",
        "min_z": "seg_mode=3d_* only: drop cells spanning fewer than this many z-planes.",
        "poly_tile_size": "Inner tiling for parallel polygon extraction. This is not the GPU tiling: it splits one label array "
                          "whose IDs are already globally unique, so fragments re-merge by ID.",
        "poly_tile_overlap": "Overlap for the inner polygon tiling.",
        "poly_tolerance": "Polygon simplification tolerance, in pixels.",
        "n_poly_workers": "Processes for polygon extraction. None = all available CPUs.",
    }

    # ------------------------------------------------------------------

    def validate(self) -> None:
        super().validate()

        if self.scale < 1:
            raise ConfigError(f"scale must be >= 1, got {self.scale}.")
        if self.n_z_planes < 1:
            raise ConfigError(f"n_z_planes must be >= 1, got {self.n_z_planes}.")
        if self.micron_per_z <= 0:
            raise ConfigError(f"micron_per_z must be > 0, got {self.micron_per_z}.")
        if self.poly_tile_overlap >= self.poly_tile_size:
            raise ConfigError(
                f"poly_tile_overlap ({self.poly_tile_overlap}) must be smaller "
                f"than poly_tile_size ({self.poly_tile_size})."
            )
        if self.min_z > self.n_z_planes:
            raise ConfigError(
                f"min_z ({self.min_z}) exceeds n_z_planes ({self.n_z_planes}); "
                f"every cell would be filtered out."
            )
        if isinstance(self.anisotropy, (int, float)) and self.anisotropy <= 0:
            raise ConfigError(f"anisotropy must be > 0, got {self.anisotropy}.")

        # mode-specific parameters that would otherwise be silently ignored
        if self.seg_mode != "3d_true" and self.anisotropy != "auto":
            raise ConfigError(
                f"anisotropy is only used when seg_mode='3d_true' "
                f"(got seg_mode={self.seg_mode!r}); leave it as 'auto'."
            )
        if self.micron_to_mosaic_path and (
            self.seg_mode != "3d_true" or self.anisotropy != "auto"
        ):
            # Warn rather than raise: unlike `anisotropy` — a numeric knob that
            # would silently change results if ignored — this is an inert path.
            logger.warning(
                "micron_to_mosaic_path=%s is set but unused: it only feeds the "
                "anisotropy='auto' derivation under seg_mode='3d_true'.",
                self.micron_to_mosaic_path,
            )
        # --- GPU tiling ---------------------------------------------------
        gpu_tile = parse_gpu_tile_size(self.gpu_tile_size)   # raises on garbage
        if gpu_tile is not None:
            if self.gpu_tile_overlap < 0:
                raise ConfigError(
                    f"gpu_tile_overlap must be >= 0, got {self.gpu_tile_overlap}."
                )
            if isinstance(gpu_tile, tuple):
                if min(gpu_tile) <= self.gpu_tile_overlap:
                    raise ConfigError(
                        f"gpu_tile_overlap ({self.gpu_tile_overlap}) must be "
                        f"smaller than every axis of gpu_tile_size {gpu_tile}; "
                        f"otherwise the grid never advances."
                    )
                if min(gpu_tile) < self.min_gpu_tile():
                    raise ConfigError(
                        f"gpu_tile_size {gpu_tile} is below the floor of "
                        f"{self.min_gpu_tile()} px "
                        f"(poly_tile_size={self.poly_tile_size} x "
                        f"scale={self.scale}): the tile would downscale to less "
                        f"than one inner polygon tile. Raise gpu_tile_size, or "
                        f"lower poly_tile_size if the small tile is deliberate."
                    )
            if not 0.0 < self.gpu_tile_iou < 1.0:
                raise ConfigError(
                    f"gpu_tile_iou must be in (0, 1), got {self.gpu_tile_iou}."
                )
            # 0 is allowed here, unlike gpu_tile_iou: it means "disabled".
            if not 0.0 <= self.gpu_tile_containment < 1.0:
                raise ConfigError(
                    f"gpu_tile_containment must be in [0, 1), got "
                    f"{self.gpu_tile_containment}."
                )
            if self.n_gpu_threads < 1:
                raise ConfigError(
                    f"n_gpu_threads must be >= 1, got {self.n_gpu_threads}."
                )
            for name in ("n_prep_workers", "n_post_workers"):
                value = getattr(self, name)
                if value is not None and value < 1:
                    raise ConfigError(f"{name} must be >= 1 or None, got {value}.")
            if self.apply_clahe:
                raise ConfigError(
                    "apply_clahe cannot be combined with gpu_tile_size: CLAHE is "
                    "a global contrast transform, and applying it per tile would "
                    "make adjacent tiles incomparable at their shared edge."
                )
        else:
            # Inert-parameter warning, same treatment as micron_to_mosaic_path.
            defaults = type(self)
            for name in ("gpu_tile_overlap", "gpu_tile_iou", "gpu_tile_containment",
                         "gpu_merge",
                         "n_gpu_threads", "n_prep_workers", "n_post_workers"):
                if getattr(self, name) != getattr(defaults, name):
                    logger.warning(
                        "%s=%s is set but unused: it only applies when "
                        "gpu_tile_size is given (the tiled path).",
                        name, getattr(self, name),
                    )

        # Mode-specific, NOT tiled-specific: keep it outside the tiled branch.
        if not 0 < self.max_size_fraction <= 1:
            raise ConfigError(f"max_size_fraction must be in (0, 1], got {self.max_size_fraction}.")
        if self.flow3d_smooth_z < 0:
            raise ConfigError(
                f"flow3d_smooth_z must be >= 0 (0 disables), got "
                f"{self.flow3d_smooth_z}."
            )
        # Warn only if EXPLICITLY set for a mode that cannot use it. Comparing
        # against the class default matters now that the default is non-zero:
        # otherwise every 3d_stitched and 2d run warns about a value it never
        # chose. Same treatment the tiled-only parameters get below.
        if (self.seg_mode != "3d_true"
                and self.flow3d_smooth_z != type(self).flow3d_smooth_z):
            logger.warning(
                "flow3d_smooth_z=%s is ignored for seg_mode=%r: cellpose applies "
                "flow3D_smooth only on the do_3D path. 3d_stitched has no ring "
                "artifact to fix.", self.flow3d_smooth_z, self.seg_mode,
            )

        if self.seg_mode == "2d" and self.min_z > 1:
            logger.info(
                "seg_mode=2d produces a single-plane mask, so the min_z=%d filter "
                "does not apply; 2D->3D replication across %d planes happens "
                "downstream in ingest_polygons.", self.min_z, self.n_z_planes,
            )

    # ------------------------------------------------------------------
    def min_gpu_tile(self) -> int:
        """
        Smallest sane outer tile edge, in full-resolution pixels.
        """
        return self.poly_tile_size * self.scale

    # ------------------------------------------------------------------
    def resolved_gpu_tile(self):
        """``None`` | ``"auto"`` | ``(tile_y, tile_x)``.

        ``None`` selects the single-call path; anything else selects the tiled
        one. Validated in :meth:`validate`, so this cannot raise here for a
        config that was built through ``from_kwargs``.
        """
        return parse_gpu_tile_size(self.gpu_tile_size)

    # ------------------------------------------------------------------
    def resolved_anisotropy(self, images_dir: str | Path) -> float | None:
        """The anisotropy to hand cellpose, or None when the mode ignores it.

        ``images_dir`` is the region's images directory, used to auto-discover the
        transform; ``micron_to_mosaic_path`` overrides it when set.
        """
        if self.seg_mode != "3d_true":
            return None
        if self.anisotropy != "auto":
            return float(self.anisotropy)
        source = self.micron_to_mosaic_path or images_dir
        native_xy = read_native_xy_um(source)
        value = derive_anisotropy(self.micron_per_z, native_xy, self.scale)
        logger.info(
            "anisotropy='auto' -> %.3f (z=%.3f um / (native_xy=%.4f um x scale=%d))",
            value, self.micron_per_z, native_xy, self.scale,
        )
        return value
