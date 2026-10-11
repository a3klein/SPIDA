"""I-module IO: tables as h5ad copies. Only ``backup_table`` needs spatialdata."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad

from spida.io.layout import table_key

logger = logging.getLogger(__name__)


def read_adata(path) -> ad.AnnData:
    """An AnnData file: ``.h5ad``, or a zarr-encoded AnnData directory."""
    path = Path(path)
    return ad.read_zarr(path) if path.is_dir() else ad.read_h5ad(path)


def load_table(
    dataset,
    prefix: str,
    suffix: str = "",
    adata_file=None,
    refresh: bool = False,
):
    """
    The table to annotate and the AnnData file it is kept in: ``(adata, h5ad_path)``.

    Parameters:
        dataset: The region's locations (Dataset).
        prefix: Segmentation name, e.g. cellpose_cell.
        suffix: Table suffix, e.g. _filt.
        adata_file: Used as-is; annotations are written back to it.
        refresh: Re-read from the zarr store even if the AnnData file exists (drops its annotations).
            Otherwise the region's file in the AnnData store is used, created from the zarr table if missing.
    """
    if adata_file is not None:
        return read_adata(adata_file), Path(adata_file)
    key = table_key(prefix, suffix)
    path = dataset.table_h5ad(key)
    if path.exists() and not refresh:
        logger.info(f"reading table copy {path}")
        return ad.read_h5ad(path), path
    logger.info(f"copying table {key} from {dataset.require('zarr')} to {path}")
    adata = ad.read_zarr(dataset.require("zarr") / "tables" / key)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return adata, path


def tag_columns(
    adata: ad.AnnData,
    columns,
    method: str,
    reference: str,
) -> dict:
    """
    Rename ``columns`` of ``adata.obs`` to ``<method>_<reference>_<rest>``, dropping a leading or
    embedded ``<method>_`` from the original name. Returns the old -> new map.

    Parameters:
        columns: obs columns to rename (the ones the annotation run added).
        method: e.g. allcools.
        reference: Reference name.
    """
    tag = f"{method}_{reference}_"
    token = f"{method}_"
    names = {c: tag + (c.replace(token, "", 1) if token in c else c) for c in columns}
    stale = [n for n in names.values() if n in adata.obs.columns]  # a re-run of the same reference
    adata.obs = adata.obs.drop(columns=stale).rename(columns=names)
    return names


def backup_table(
    dataset,
    adata: ad.AnnData,
    key: str,
) -> None:
    """
    Replace a table in the region's zarr store (needs spatialdata).

    Parameters:
        dataset: The region's locations (Dataset).
        key: Table name in the zarr store, e.g. cellpose_cell_table_filt.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.utilities.sd_utils import _backup_adata
    _backup_adata(dataset.require("zarr"), adata, key)
