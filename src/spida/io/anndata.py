"""
Reading and writing SPIDA's AnnData files: .h5ad, or a zarr-encoded AnnData directory (anndata only, no spatialdata).

read_adata reads a file, read_obs reads only its obs.
write_adata writes an .h5ad useing a temp-swap strategy.
"""

import logging
import os
from pathlib import Path

import anndata as ad
import pandas as pd

logger = logging.getLogger(__package__)

TMP_SUFFIX = "__tmp"

# safe read anndata.
def read_adata(
    path: str | Path,
    backed: str | None = None,
) -> ad.AnnData:
    """
    An AnnData file: ``.h5ad`` (optionally backed), or a zarr-encoded AnnData directory.
    Parameters:
    path (str | Path): The file or directory.
    backed (str, optional): Passed to anndata.read_h5ad (ignored for a zarr directory).
    """
    path = Path(path)
    if path.is_dir():
        return ad.read_zarr(path)
    return ad.read_h5ad(path, backed=backed)

# read just the obs dataframe
def read_obs(
    path: str | Path,
) -> pd.DataFrame:
    """
    Only the obs of an AnnData file (an .h5ad is opened backed and closed again).
    Parameters:
    path (str | Path): The file or directory.
    """
    path = Path(path)
    if path.is_dir():
        return ad.read_zarr(path).obs
    a = ad.read_h5ad(path, backed="r")
    try:
        return a.obs.copy()
    finally:
        a.file.close()

# safe write anndata using temp-swap. 
def write_adata(
    path: str | Path,
    adata: ad.AnnData,
) -> None:
    """
    Write ``adata`` as an .h5ad: to ``<path>__tmp.h5ad`` first, then swap, so a crash never leaves a
    truncated file at ``path``. Creates the parent directory. Replaces an existing file.
    Parameters:
    path (str | Path): The .h5ad to write.
    adata (AnnData): The table.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.stem}{TMP_SUFFIX}{path.suffix or '.h5ad'}")
    try:
        adata.write_h5ad(tmp)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
