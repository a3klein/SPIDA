"""Shared imaging config: the deconwolf ini and the optics it declares.

Read by every transform that needs physical parameters, so NA / pixel size /
wavelength have one source rather than one per filter.
"""
import configparser
from pathlib import Path

from spida._constants import OPTICS_DEFAULTS


def _configure_deconwolf(config_path: str | Path | None):
    """Deconwolf binary and PSF dir from the ``[Paths]`` section of ``config_path``."""
    if config_path is None:
        raise ValueError("no deconwolf config; pass --deconwolf_config_file or set DECONWOLF_CONFIG in the config")
    config = configparser.ConfigParser()
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found at {config_path}. Please create a .config.ini file."
        )
    config.read(str(config_path))
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


def get_optics(channel: str | None = None, config_path: str | Path | None = None) -> dict:
    """Imaging optics from the [Optics] section of the deconwolf ini.

    Any consumer needing physical parameters -- deconvolution, PSF generation, or a
    background-removal filter that takes wavelength/NA/pixel size -- should read them
    from here rather than hardcoding, so there is a single value to correct when one
    is wrong. Missing keys fall back to `OPTICS_DEFAULTS`; a missing section is not an
    error, so this does not break an existing deployment whose ini predates it.

    `channel` adds `emission_nm` for that channel when the ini declares it.

    Lives here rather than in `spida.config` because it reads the SAME ini file as
    `_configure_deconwolf` above (its path is `SpidaConfig.deconwolf_config_file`);
    `spida.config` handles SPIDA's own JSON config. No ``config_path`` gives the
    defaults, which are in `spida._constants` (no I/O).
    """
    config = configparser.ConfigParser()
    optics = dict(OPTICS_DEFAULTS)
    if config_path is not None and Path(config_path).is_file():
        config.read(str(config_path))
        if config.has_section("Optics"):
            for k in OPTICS_DEFAULTS:
                if config.has_option("Optics", k):
                    optics[k] = config.getfloat("Optics", k)
            if channel and config.has_option("Optics", f"emission_nm_{channel}"):
                optics["emission_nm"] = config.getfloat(
                    "Optics", f"emission_nm_{channel}")
    return optics
