import os
import sys
from dotenv import load_dotenv # type: ignore
load_dotenv()

import configparser
import subprocess
from pathlib import Path

import numpy as np
import imageio.v3 as iio


def _configure_deconwolf():
    """Configure Deconwolf."""
    config = configparser.ConfigParser()
    config_path = os.getenv("DECONWOLF_CONFIG", "/anvil/projects/x-mcb130189/aklein/BICAN/fish_supp/.config.ini")
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at {config_path}. Please create a .config.ini file.")
    config.read(config_path)
    dw_path = Path(config.get("Paths", "dw_path"))
    if not dw_path.exists():
        raise FileNotFoundError(
            f"Deconwolf not found at {dw_path}. To use Deconwolf, please install it and set the path in the .config.ini file."
        )
    psf_path = Path(config.get("Paths", "dw_psf_path"))
    if not psf_path.exists():
        raise FileNotFoundError(
            f"Deconwolf PSF not found at {psf_path}. To use Deconwolf, please compute PSFs using `dw_bw`"
        )
    return dw_path, psf_path


def deconwolf(
    img_path: str | Path,
    colors: list | np.ndarray | int, #TODO: This runs on only one color at a time. 
    decon_path: str | Path = None,
    z_step: float = 0.6,
    iter: int = 100,
    gpu: bool = False,
    tilesize: int = None,
    lock = None,
) -> np.ndarray:
    """Deconwolf deconvolution.

    Parameters
    ----------
    img
        A (Z,X,Y) or (C,Z,X,Y) image.
    colors
        A list of channel colors. Must be the same length as the number of channels in the image.
    z_step
        Z step size in microns.
    iter
        Number of iterations.
    gpu
        Whether to use the GPU.
    tilesize
        Size of tiles. If None, no tiling is used.

    Returns
    -------
    filtered
        a numpy array of the same shape and dtype as the input image.
    """
    # Setup
    dw_path, psf_path = _configure_deconwolf()
    if isinstance(img_path, str):
        img_path = Path(img_path)
    img = iio.imread(img_path)
    if decon_path is None:
        decon_path = img_path.with_suffix(".decon.tif")
    
    has_channels = len(img.shape) == 4
    if not has_channels:
        img = np.expand_dims(img, axis=0)        

    if isinstance(colors, int):
        colors = [colors]
    # Run deconwolf on each channel
    gpu = "--gpu" if gpu else ""
    tile = f"--tilesize {tilesize}" if tilesize else ""
    for i, color in enumerate(colors):
        color_psf = psf_path / f"{color}_z{int(z_step*1000)}_psf.tiff"
        if not color_psf.exists():
            raise FileNotFoundError(f"PSF for color {color} and z_step {z_step} not found at {color_psf}.")

        command = f"{dw_path} --out {decon_path} --iter {iter} {tile} {gpu} {img_path} {color_psf} --verbose 0"
        # if "lock" in globals():
        if lock is not None:
            with lock: 
                subprocess.run(command, check=True, shell=True)
        else: 
            subprocess.run(command, check=True, shell=True)

    return decon_path
