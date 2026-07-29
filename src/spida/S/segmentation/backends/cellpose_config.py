"""Cellpose run configuration + model resolution.

Two concerns, both cellpose-specific and both deliberately free of any
``cellpose`` import at module scope so this file stays importable from *any*
pixi env (``segment-region cellpose --list-params`` must work from
``preprocessing``, and the tests must run without a GPU):

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
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Literal

import numpy as np

from .config import BackendConfig, ConfigError

logger = logging.getLogger(__package__)

#: Environment variable cellpose itself honours for its model cache.
CELLPOSE_ENV_VAR = "CELLPOSE_LOCAL_MODELS_PATH"
#: SPIDA's config key for the same thing (config JSON / .env).
SPIDA_ENV_VAR = "CELLPOSE_MODEL_PATH"

#: Models SPIDA knows about, with the notes that should inform picking one.
KNOWN_MODELS: dict[str, str] = {
    "cpsam": "Cellpose-SAM (ViT-L). SPIDA production default.",
    "cpsam_v2": "Improved SAM, better on low contrast. Costless drop-in for "
                "cpsam (same VRAM and runtime, ~95% agreement).",
    "cpdino": "DINOv3 ViT-L backbone. Needs the `dinov3` package.",
    "cpdino-vitb": "DINOv3 ViT-B backbone. ~2.2x faster and lowest VRAM, but "
                   "finds ~21% more cells than cpsam and over-segments faint "
                   "regions — not biologically validated. Needs `dinov3`.",
}

#: Models whose weights cannot be loaded without the DINOv3 package.
REQUIRES_DINOV3: frozenset[str] = frozenset({"cpdino", "cpdino-vitb"})


def resolve_model_dir(explicit: str | Path | None = None) -> Path:
    """Resolve the directory holding cellpose weights.

    Precedence: explicit config field > ``CELLPOSE_MODEL_PATH`` (SPIDA's config
    key) > ``CELLPOSE_LOCAL_MODELS_PATH`` (cellpose's own) > cellpose's default
    ``~/.cellpose/models``.
    """
    for candidate in (explicit, os.getenv(SPIDA_ENV_VAR), os.getenv(CELLPOSE_ENV_VAR)):
        if candidate:
            return Path(candidate).expanduser()
    return Path.home() / ".cellpose" / "models"


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
                "now. Set %s to a directory holding the weights to avoid "
                "re-downloading on every job.",
                self.name, self.model_dir, SPIDA_ENV_VAR,
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


def read_native_xy_um(images_dir: str | Path) -> float:
    """Micron size of one full-resolution mosaic pixel, from the region's
    ``micron_to_mosaic_pixel_transform.csv`` (a 3x3 micron->pixel matrix, so
    element [0][0] is pixels per micron in x)."""
    path = Path(images_dir) / "micron_to_mosaic_pixel_transform.csv"
    if not path.is_file():
        raise ConfigError(
            f"cannot derive anisotropy: {path} not found. Pass an explicit "
            f"--anisotropy=<float> instead."
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
    micron_per_z: float = 1.5
    n_z_planes: int = 7

    # --- cellpose eval ----------------------------------------------------
    flow_threshold: float = 0.0
    cellprob_threshold: float = -2.0
    tile_norm_blocksize: int = 0
    diameter: float | None = None
    cp_min_size: int = 100
    batch_size: int = 8
    stitch_threshold: float = 0.2

    # --- inputs -----------------------------------------------------------
    nuc_stain_name: str = "DAPI"
    cyto_stain_name: str | None = None
    image_ext: str = ".tif"

    # --- optional preprocessing -------------------------------------------
    apply_clahe: bool = False
    clahe_clip_limit: float = 20.0
    clahe_tile_grid_size: int = 8

    # --- polygon post-processing ------------------------------------------
    min_size: int = 100
    min_z: int = 3
    poly_tile_size: int = 1000
    poly_tile_overlap: int = 200
    poly_tolerance: float = 0.5
    n_poly_workers: int | None = None

    ALIASES: ClassVar[dict[str, str]] = {
        "model_name": "model",
        "tile_size": "poly_tile_size",
        "overlap": "poly_tile_overlap",
        "clahe_clipLimit": "clahe_clip_limit",
        "clahe_tileGridSize": "clahe_tile_grid_size",
    }

    CHOICES: ClassVar[dict[str, tuple]] = {
        "seg_mode": ("2d", "3d_stitched", "3d_true"),
        "z_reduce": ("middle", "mip"),
        "model": tuple(KNOWN_MODELS),
    }

    HELP: ClassVar[dict[str, str]] = {
        "seg_mode": "2d = single z-slice broadcast across z; 3d_stitched = 2D "
                    "per slice linked by stitch_threshold (SPIDA production "
                    "default); 3d_true = joint 3D flows using anisotropy.",
        "z_reduce": "For seg_mode=2d only: take the middle z-slice, or a "
                    "maximum-intensity projection over z.",
        "model": "Cellpose model NAME (not a path). "
                 + " | ".join(f"{k}: {v}" for k, v in KNOWN_MODELS.items()),
        "model_dir": "Directory holding the weights. Defaults to "
                     f"${SPIDA_ENV_VAR}, then ${CELLPOSE_ENV_VAR}, then "
                     "~/.cellpose/models.",
        "use_bfloat16": "Half-precision inference. ~1.85x faster at mean "
                        "per-cell IoU 0.985 vs fp32.",
        "scale": "xy downsample factor applied before inference. Coupled to "
                 "anisotropy.",
        "anisotropy": "seg_mode=3d_true only. 'auto' derives z_step / "
                      "(native_xy * scale) from the region's "
                      "micron_to_mosaic transform, which is the physically "
                      "correct value. Lowering it is measurably all cost and "
                      "no benefit.",
        "micron_per_z": "Micron thickness of one z-step.",
        "n_z_planes": "Number of z-planes in the input stack.",
        "tile_norm_blocksize": "Cellpose normalisation block size (0 = global).",
        "diameter": "Expected cell diameter in pixels; None lets cellpose decide.",
        "cp_min_size": "Cellpose's INTERNAL minimum mask size, in downsampled "
                       "pixels. Distinct from --min_size.",
        "batch_size": "Measured to have no effect on speed or output; exposed "
                      "for reproducibility only.",
        "stitch_threshold": "seg_mode=3d_stitched only: IoU for linking masks "
                            "across adjacent z-slices.",
        "cyto_stain_name": "Second channel (e.g. PolyT). Cheap to add (~1% GPU) "
                           "but changes the target from nucleus to whole cell.",
        "apply_clahe": "Contrast-limited adaptive histogram equalisation before "
                       "inference.",
        "min_size": "Minimum POLYGON area, in full-resolution pixels squared, "
                    "applied after upscaling. Distinct from --cp_min_size.",
        "min_z": "seg_mode=3d_* only: drop cells spanning fewer than this many "
                 "z-planes.",
        "poly_tile_size": "Inner tiling for parallel polygon extraction. This "
                          "is NOT the GPU tiling: it splits one label array "
                          "whose IDs are already globally unique, so fragments "
                          "re-merge by ID.",
        "poly_tile_overlap": "Overlap for the inner polygon tiling.",
        "poly_tolerance": "Polygon simplification tolerance, in pixels.",
        "n_poly_workers": "Processes for polygon extraction. None = all "
                          "available CPUs.",
    }

    # ------------------------------------------------------------------
    @classmethod
    def from_kwargs(cls, **kwargs) -> "CellposeConfig":
        """As the base, but first folds the legacy mode flags into ``seg_mode``.

        ``do_3D`` and ``project_3d_to_2d`` are not simple renames — they select
        a mode — so they are translated here rather than via ``ALIASES``.
        """
        legacy = {}
        for old in ("do_3D", "project_3d_to_2d"):
            if old in kwargs:
                legacy[old] = kwargs.pop(old)

        if legacy and "seg_mode" in kwargs:
            raise ConfigError(
                f"seg_mode cannot be combined with the deprecated flag(s) "
                f"{', '.join(sorted(legacy))}; use seg_mode alone."
            )
        if legacy:
            if legacy.get("project_3d_to_2d"):
                kwargs["seg_mode"] = "2d"
                kwargs.setdefault("z_reduce", "mip")
            elif legacy.get("do_3D"):
                kwargs["seg_mode"] = "3d_true"
            else:
                kwargs["seg_mode"] = "3d_stitched"
            logger.warning(
                "CellposeConfig: %s is deprecated; interpreted as seg_mode=%r",
                ", ".join(sorted(legacy)), kwargs["seg_mode"],
            )

        return super().from_kwargs(**kwargs)  # type: ignore[return-value]

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
        if self.seg_mode == "2d" and self.min_z > 1:
            logger.info(
                "seg_mode=2d: the mask is broadcast across all %d z-planes, so "
                "min_z=%d filters nothing.", self.n_z_planes, self.min_z,
            )

    # ------------------------------------------------------------------
    def resolved_anisotropy(self, images_dir: str | Path) -> float | None:
        """The anisotropy to hand cellpose, or None when the mode ignores it."""
        if self.seg_mode != "3d_true":
            return None
        if self.anisotropy != "auto":
            return float(self.anisotropy)
        native_xy = read_native_xy_um(images_dir)
        value = derive_anisotropy(self.micron_per_z, native_xy, self.scale)
        logger.info(
            "anisotropy='auto' -> %.3f (z=%.3f um / (native_xy=%.4f um x scale=%d))",
            value, self.micron_per_z, native_xy, self.scale,
        )
        return value
