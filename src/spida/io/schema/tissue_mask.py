"""
The tissue mask: a boolean raster in image-pixel space (aligned to the mosaic), written by
``tissue-mask`` and read by segmentation and deconvolution to skip out-of-tissue tiles.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

FILENAME = "tissue_mask.tif"


def load(
    path: str | Path,
) -> np.ndarray:
    """
    The mask as a full-resolution boolean array. Written deflate-compressed, so it cannot be memmap'd.
    Parameters:
    path (str | Path): The mask file.
    """
    import tifffile

    mask = tifffile.imread(str(path))
    if mask.ndim != 2:
        raise ValueError(f"{path}: expected a 2D mask, got shape {mask.shape}")
    return mask > 0


def write(
    path: str | Path,
    mask: np.ndarray,
) -> Path:
    """
    Write a boolean mask as an 8-bit (0 / 255) deflate-compressed TIFF.
    Parameters:
    path (str | Path): The mask file.
    mask (np.ndarray): 2D boolean mask.
    """
    import tifffile

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, mask.astype(np.uint8) * 255, compression="deflate")
    return path
