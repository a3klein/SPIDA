"""Richardson-Lucy deconvolution via the deconwolf binary."""
import logging
import subprocess
from pathlib import Path

from spida.S.image_prep.config import _configure_deconwolf
from spida.S.image_prep.transforms.base import TileContext, check_params

logger = logging.getLogger(__package__)

NAME = "deconwolf"
USES_GPU = True
DEFAULTS = {"iter": 100}


def validate(params: dict) -> dict:
    """Drop psigma when tiling is on: deconwolf ignores it there, and tilesize is
    what bounds VRAM, so keeping tilesize is the safer of the two."""
    params = check_params(run, params, NAME)
    if params.get("psigma") is not None and params.get("tilesize"):
        logger.warning(
            "deconwolf ignores --psigma when --tilesize is set: dropping psigma=%s, "
            "keeping tilesize=%s. To use psigma, omit tilesize.",
            params["psigma"], params["tilesize"],
        )
        params = {k: v for k, v in params.items() if k != "psigma"}
    return params


def run(
    src: Path,
    dst: Path,
    ctx: TileContext,
    iter: int = 100,
    tilesize: int = None,
    psigma: float = None,
) -> Path:
    """Deconvolve one tile to `dst`.

    tilesize
        deconwolf's internal tiling, within the tile it is given. None is untiled.
    psigma
        Regulariser for the Richardson-Lucy sparsity artifact. Requires no tilesize.
    """
    dw_path, psf_path = _configure_deconwolf()
    src, dst = Path(src), Path(dst)

    colors = [ctx.colors] if isinstance(ctx.colors, int) else ctx.colors
    gpu_flag = "--gpu" if ctx.gpu else ""
    tile = f"--tilesize {tilesize}" if tilesize else ""
    ps = f"--psigma {psigma}" if psigma is not None else ""

    for color in colors:
        color_psf = psf_path / f"{color}_z{int(ctx.z_step_um * 1000)}_psf.tiff"
        if not color_psf.exists():
            raise FileNotFoundError(
                f"PSF for color {color} and z_step {ctx.z_step_um} not found at {color_psf}."
            )

        command = (
            f"{dw_path} --out {dst} --iter {iter} {tile} {ps} {gpu_flag} "
            f"{src} {color_psf} --verbose 0"
        )
        if ctx.gpu_lock is not None:
            with ctx.gpu_lock:
                subprocess.run(command, check=True, shell=True)
        else:
            subprocess.run(command, check=True, shell=True)

    return dst
