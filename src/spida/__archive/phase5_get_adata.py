"""Archived in phase 5 (2026-10-25): replaced by ``spida.io.store.read_table`` (anndata-only read,
verified identical on every p3 PU / THM2 table). Was ``spida.utilities.sd_utils._get_adata``."""
import logging
import warnings

import anndata as ad

logger = logging.getLogger(__package__)


def _get_adata(zarr_path, table_key: str) -> ad.AnnData:
    """
    Read one table from the store at ``zarr_path``.
    Parameters:
    zarr_path (str | Path): The dataset's zarr store.
    table_key (str): The table to read.
    """
    ### TODO: Try, except block for when spatialdata is not installed in a given environment
    ### Then load the AnnData object directly from the zarr store

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd

    logger.info(f"Loading {table_key} from {zarr_path}")
    sdata = sd.read_zarr(zarr_path)
    adata = sdata[table_key].copy()
    return adata
