"""Per-tile transforms, selected by `--transform`.

Each is a module satisfying the Transform protocol in `base`: NAME, USES_GPU,
DEFAULTS, validate() and run(src, dst, ctx, **params). Adding one means writing the
module and listing it here.
"""
from . import dark_sectioning, deconwolf
from .base import TileContext, Transform, check_params

TRANSFORMS = {
    deconwolf.NAME: deconwolf,
    dark_sectioning.NAME: dark_sectioning,
}


def resolve_transform(name: str):
    """Look up a `--transform` name, or fail with what exists."""
    try:
        return TRANSFORMS[name]
    except KeyError:
        raise ValueError(
            f"unknown transform {name!r}; available: {sorted(TRANSFORMS)}"
        ) from None


__all__ = [
    "TRANSFORMS",
    "resolve_transform",
    "TileContext",
    "Transform",
    "check_params",
    "deconwolf",
    "dark_sectioning",
]
