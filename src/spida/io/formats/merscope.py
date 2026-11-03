"""
MERSCOPE output: the micron -> mosaic-pixel transform, the detected transcripts and the mosaic TIFF
naming (one ``mosaic_<stain>_z<N>.tif`` per stain per plane, under the region's ``images/``).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

from spida.io.images import probe_mosaic_shape

logger = logging.getLogger(__package__)

IMAGES_DIR = "images"
TRANSFORM_FILE = "micron_to_mosaic_pixel_transform.csv"
TRANSCRIPTS_PARQUET = "detected_transcripts.parquet"
TRANSCRIPTS_CSV = "detected_transcripts.csv"
# one raw plane image; decon outputs (``.decon.tif``) do not match
PLANE_RE = re.compile(r"mosaic_(?P<stain>[\w-]+)_z(?P<z>\d+)\.tif$")


def raw_transcripts(
    raw_dir: str | Path,
) -> Path:
    """
    The region's detected transcripts: the parquet if present, else the csv.
    Parameters:
    raw_dir (str | Path): The raw MERSCOPE region dir.
    """
    raw_dir = Path(raw_dir)
    parquet = raw_dir / TRANSCRIPTS_PARQUET
    return parquet if parquet.exists() else raw_dir / TRANSCRIPTS_CSV


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


# MERSCOPE ULTRA data organization file
def read_info(path: str | Path) -> pd.DataFrame:
    """
    Read the MERSCOPE ULTRA data organization file

    Parameters
    ----------
    path : str or Path
        Path to the data organization file (CSV format)

    Returns
    -------
    pd.DataFrame with columns channel, prefix, round, color, frame, zPos, file_pattern
    """
    import pandas as pd

    # Read in data and get used columns
    df_data_org = pd.read_csv(path)
    df_info = df_data_org[
        ["channelName", "imageType", "imagingRound", "color", "frame", "zPos"]
    ]
    df_info = df_info.rename(
        columns={
            "channelName": "channel",
            "imagingRound": "round",
            "imageType": "prefix",
        }
    )
    # assign types
    df_info = df_info.astype(
        {
            "channel": str,
            "prefix": str,
            "round": str,
            "color": int,
            "frame": "object",
            "zPos": "object",
        }
    )
    # format round
    df_info["round"] = "_" + df_info["round"]
    df_info.loc[df_info["round"] == "_-1", "round"] = ""
    # formate frame
    df_info["frame"] = df_info["frame"].apply(
        lambda x: np.asarray(x.strip("[]").split(","), dtype=int)
    )
    df_info["zPos"] = df_info["zPos"].apply(
        lambda x: np.asarray(x.strip("[]").split(","), dtype=float)
    )
    # get file pattern
    df_info["file_pattern"] = (
        df_info["prefix"] + df_info["round"] + "_{fov:04d}.jp2.zip"
    )

    return df_info
