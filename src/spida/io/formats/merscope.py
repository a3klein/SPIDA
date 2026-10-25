"""
MERSCOPE output: the micron -> mosaic-pixel transform and the mosaic TIFF naming
(one ``mosaic_<stain>_z<N>.tif`` per stain per plane, under the region's ``images/``).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from spida.io.images import probe_mosaic_shape

logger = logging.getLogger(__package__)

TRANSFORM_FILE = "micron_to_mosaic_pixel_transform.csv"


def read_transform(path: str | Path) -> np.ndarray:
    """
    The 3x3 micron -> mosaic-pixel affine (whitespace-delimited rows).
    Parameters:
    path (str | Path): The transform file, or the region's images dir.
    """
    path = Path(path)
    if path.is_dir():
        path = path / TRANSFORM_FILE
    rows = [list(map(float, ln.split())) for ln in path.read_text().splitlines() if ln.strip()]
    if len(rows) != 3 or any(len(r) != 3 for r in rows):
        raise ValueError(f"{path}: micron-to-mosaic transform must be a 3x3 matrix")
    return np.asarray(rows, dtype=float)


def find_stain_planes(
    images_dir: str | Path,
    stain: str,
    image_ext: str = ".tif"
) -> list[Path]:
    """
    z-ordered plane paths for one stain.

    Parameters: 
        images_dir - directory containing the mosaic TIFFs
        stain - substring to match in the stain name
        image_ext - file extension of the mosaic TIFFs (default: ".tif")
    Returns: 
        list of Path objects for the mosaic TIFFs, sorted in natural order
    """
    from natsort import natsorted

    directory = Path(images_dir)
    files = [
        f for f in directory.glob("*" + image_ext)
        if "_masks" not in f.name and "_flows" not in f.name and stain in f.name
    ]
    if not files:
        raise ValueError(
            f"no planes for stain {stain!r} with extension {image_ext!r} in {directory}"
        )
    return list(natsorted(files))


def raw_mosaic_shape(
    images_dir: str | Path,
    stain: str = "DAPI",
) -> tuple[int, int]:
    """
    (rows, cols) of a raw mosaic plane (decon outputs and z-stacks skipped), from the header only.
    Parameters:
    images_dir (str | Path): The region's images dir.
    stain (str): Stain to probe; any stain if it has no plane.
    """
    images_dir = Path(images_dir)
    raw = lambda p: ".decon." not in p.name and "_z_stack_" not in p.name
    cands = sorted(p for p in images_dir.glob(f"mosaic_{stain}_z*.tif") if raw(p))
    if not cands:
        cands = sorted(p for p in images_dir.glob("mosaic_*_z*.tif") if raw(p))
    return probe_mosaic_shape(cands[0])
