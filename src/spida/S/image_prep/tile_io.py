"""Tile-level TIFF writing and z-projection."""
import logging
from pathlib import Path

import imageio.v3 as iio
import numpy as np
import tifffile

logger = logging.getLogger(__package__)


BIGTIFF_BYTE_THRESHOLD = (2**32) - 1


def write_tiff(
    output_path: str | Path,
    image: np.ndarray,
    *,
    bigtiff_threshold: int = BIGTIFF_BYTE_THRESHOLD,
):
    """Write a TIFF, using BigTIFF when the payload exceeds 4GB."""
    output_path = Path(output_path)
    bigtiff = image.nbytes > bigtiff_threshold
    tifffile.imwrite(output_path, image, bigtiff=bigtiff)


def project_down_2D(input_file, output_file: str | Path = None):
    """
    Project a 3D image to 2D by taking the maximum projection along the z-axis.
    """
    if output_file is None:
        output_file = input_file.with_suffix(".2d.tif")
    decon_img = iio.imread(input_file)

    if decon_img.ndim == 3 and decon_img.shape[0] > 1:
        # Remove first and last z slices if there are multiple z-slices
        if decon_img.shape[0] > 2:
            decon_img = decon_img[1:-1]
        # Z-project to 2D using max projection
        decon_img_2d = decon_img.max(axis=0)
    else:
        decon_img_2d = decon_img.squeeze()

    # Save the 2D projected tile
    iio.imwrite(output_file, decon_img_2d)

def project_to_3D(input_file, output_file: str | Path = None, z_stack_size: int = None):
    """
    Project a 3D padded image (so with repetitive edge slices) to 3D by taking the middle z_stack_size slices along the z-axis.
    """
    if z_stack_size is None:
        raise ValueError("project_to_3D requires z_stack_size")
    if output_file is None:
        output_file = input_file.with_suffix(".3d.tif")
    # Read only the needed pages one at a time to avoid loading the full padded volume.
    try:
        with tifffile.TiffFile(input_file) as tf:
            n_slices = len(tf.pages)
            start = (n_slices - z_stack_size) // 2
            end = start + z_stack_size
            if start < 0 or end > n_slices:
                raise ValueError(
                    f"TIFF has {n_slices} pages but z_stack_size={z_stack_size}"
                )
            decon_img_3d = np.stack([tf.pages[i].asarray() for i in range(start, end)])
    except Exception:
        logger.warning(
            "project_to_3D: page-level read failed for %s — falling back to full load",
            input_file,
        )
        decon_img = iio.imread(input_file)
        if decon_img.ndim == 3:
            start = (decon_img.shape[0] - z_stack_size) // 2
            end = start + z_stack_size
            decon_img_3d = decon_img[start:end].copy()
            del decon_img
        else:
            decon_img_3d = decon_img.squeeze()

    # Save the 3D projected tile
    iio.imwrite(output_file, decon_img_3d)
