"""
Deconvolved mosaic planes as written by ``S/image_prep`` (``transform-image``) next to the raw planes:
2D ``mosaic_<stain>_z<N>.decon.tif`` (from ``mosaic_<stain>_z<N>.tif``), 3D
``mosaic_<stain>_z_stack_z<i>.decon.tif`` (i = position in the stack, 0-based).
"""

from __future__ import annotations

SUFFIX = ".decon.tif"
Z_STACK_TAG = "z_stack"


def plane_name(
    stain: str,
    z: int | str,
    suffix: str = SUFFIX,
    z_stack: bool = False,
) -> str:
    """
    File name of one deconvolved plane.
    Parameters:
    stain (str): Stain, e.g. "DAPI".
    z (int | str): Plane (2D: the raw plane's z; 3D: the position in the stack).
    suffix (str): File suffix (default ".decon.tif").
    z_stack (bool): 3D (z-stack) naming.
    """
    return f"mosaic_{stain}_{Z_STACK_TAG}_z{z}{suffix}" if z_stack else f"mosaic_{stain}_z{z}{suffix}"
