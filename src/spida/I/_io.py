"""I-module IO: tables as h5ad copies. Only ``backup_table`` needs spatialdata; ``sync_annotations`` does not."""

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
    if path.exists():  # --refresh_table: the rebuilt copy has none of this file's annotations
        old = ad.read_h5ad(path, backed="r")
        lost = [c for c in old.obs.columns if c.startswith(("mmc_", "allcools_"))]
        old.file.close()
        if lost:
            logger.warning("--refresh_table: rebuilding %s drops %d annotation columns: %s",
                           path.name, len(lost), ", ".join(lost))
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
    import zarr
    from anndata.io import read_elem
    from spida.io.store import is_writeable, update_obs

    key = table_key(prefix, suffix)
    zarr_path = dataset.require("zarr")
    src_path = Path(adata_file) if adata_file is not None else dataset.table_h5ad(key)
    if src_path.is_dir():
        src = ad.read_zarr(src_path).obs
    else:
        a = ad.read_h5ad(src_path, backed="r")
        src = a.obs.copy()
        a.file.close()

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

    obs = read_elem(zarr.open_group(str(zarr_path / "tables" / key), mode="r")["obs"])
    if extra := len(src.index.difference(obs.index)):
        raise ValueError(f"{extra} cells in {src_path.name} are not in the zarr table {key}")
    if absent := len(obs.index.difference(src.index)):
        logger.warning("%d cells of %s are not in %s; their synced values are missing", absent, key, src_path.name)
    for c in cols:
        col = src[c].reindex(obs.index)
        obs[c] = col.astype("boolean") if absent and src[c].dtype == bool else col  # bool + missing
    update_obs(zarr_path, key, obs)
    logger.info("synced %d columns into %s: %s", len(cols), key, ", ".join(cols))
    return cols
