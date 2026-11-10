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
from spida.io.schema.decon import SUFFIX as DECON_SUFFIX, Z_STACK_TAG

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

    Matches whole names only: raw ``mosaic_<stain>_z<N>.tif``, or decon (``io/schema/decon``)
    ``mosaic_<stain>[_z_stack]_z<N>.decon.tif``. The stain must match exactly (DAPI does not pick up
    DAPI2), and the name must end in ``image_ext``: ".tif" selects raw planes only, an ``image_ext``
    containing ".decon" decon planes only (e.g. "_z3.decon.tif" -> plane 3 alone).
    Raises if no plane matches, or if a z appears twice (e.g. 2D and z-stack decon names side by side).

    Parameters:
        images_dir - directory containing the mosaic TIFFs
        stain - stain name, matched exactly
        image_ext - required ending of the file name (default: ".tif")
    Returns:
        list of Path objects for the mosaic TIFFs, sorted by z
    """
    directory = Path(images_dir)
    if ".decon" in image_ext:
        pat = re.compile(rf"mosaic_{re.escape(stain)}(?:_{Z_STACK_TAG})?_z(\d+){re.escape(DECON_SUFFIX)}")
    else:
        pat = re.compile(rf"mosaic_{re.escape(stain)}_z(\d+)\.tif")
    by_z: dict[int, list[Path]] = {}
    for f in directory.iterdir():
        m = pat.fullmatch(f.name)
        if m and f.name.endswith(image_ext):
            by_z.setdefault(int(m[1]), []).append(f)
    if not by_z:
        hint = ""
        if ".decon" not in image_ext and any(directory.glob(f"mosaic_{stain}_*{DECON_SUFFIX}")):
            hint = f" (it holds deconvolved planes: pass image_ext={DECON_SUFFIX!r})"
        raise ValueError(
            f"no planes for stain {stain!r} with extension {image_ext!r} in {directory}{hint}"
        )
    dup = {z: sorted(p.name for p in fs) for z, fs in by_z.items() if len(fs) > 1}
    if dup:
        raise ValueError(f"{directory}: more than one {stain!r} file per plane: {dup}")
    return [by_z[z][0] for z in sorted(by_z)]


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
