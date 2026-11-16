"""I-module IO: tables as h5ad copies, and sync-annotations back into the zarr table's obs (no spatialdata)."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad

from spida.io.anndata import read_adata, read_obs, write_adata  # noqa: F401  (read_adata re-exported)
from spida.io.layout import table_key

logger = logging.getLogger(__name__)


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
        return read_adata(path), path
    if path.exists():  # --refresh_table: the rebuilt copy has none of this file's annotations
        lost = [c for c in read_obs(path).columns if c.startswith(("mmc_", "allcools_"))]
        if lost:
            logger.warning("--refresh_table: rebuilding %s drops %d annotation columns: %s",
                           path.name, len(lost), ", ".join(lost))
    logger.info(f"copying table {key} from {dataset.require('zarr')} to {path}")
    from spida.io.store import read_table
    adata = read_table(dataset.require("zarr"), key)
    write_adata(path, adata)
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


def sync_annotations(
    dataset,
    prefix: str,
    suffix: str = "_filt",
    adata_file=None,
    reference: str | None = None,
    columns=None,
    overwrite: bool = False,
) -> list[str]:
    """
    Copy annotation columns from the table's AnnData file into the zarr table's obs (obs only; X,
    layers, obsm and uns in the zarr store are not touched). Needs zarr + anndata, not spatialdata.
    Returns the synced column names.

    Parameters:
        dataset: The region's locations (Dataset).
        prefix: Segmentation name, e.g. cellpose_cell.
        suffix: Table suffix, e.g. _filt.
        adata_file: AnnData to sync from (default: the region's file in the AnnData store).
        reference: Only this reference's columns (mmc_<reference>_* / allcools_<reference>_*).
        columns: Exactly these obs columns instead (default: every mmc_* / allcools_* column).
        overwrite: Replace columns that already exist in the zarr table.
    """
    from spida.io.store import is_writeable, require, sync_obs

    key = table_key(prefix, suffix)
    zarr_path = dataset.require("zarr")
    src_path = Path(adata_file) if adata_file is not None else dataset.table_h5ad(key)
    require(zarr_path, [key], files=[src_path])
    src = read_obs(src_path)

    if columns:
        cols = list(columns)
        if missing := [c for c in cols if c not in src.columns]:
            raise KeyError(f"{src_path.name} has no columns {missing}")
    else:
        tags = (f"mmc_{reference}_", f"allcools_{reference}_") if reference else ("mmc_", "allcools_")
        cols = [c for c in src.columns if c.startswith(tags)]
    if not cols:
        raise ValueError(f"no annotation columns to sync in {src_path}")
    is_writeable(zarr_path, columns={key: cols}, overwrite=overwrite)
    return sync_obs(zarr_path, key, src, cols, source=src_path.name)
