"""
Reading and writing a zarr store.

read_table reads one table with anndata only (no spatialdata needed); open_store opens the whole store.
is_writeable checks up front for overwrites (pathlib only).
write_element replaces an existing element by temp-swap.
update_obs rewrites only a table's (anndata) obs (no spatialdata needed).
"""

import logging
import warnings
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__package__)

ELEMENT_KINDS = ("images", "labels", "points", "shapes", "tables")
TMP_SUFFIX = "__tmp" # where to store temp files for overwrites. 


# Checks if an element is present on disk.
def element_on_disk(
    zarr_path: str | Path,
    key: str,
) -> str | None:
    """
    Where element ``key`` is in the store, e.g. "tables/cellpose_cell_table", or None.
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    key (str): Element name.
    """
    for kind in ELEMENT_KINDS:
        if (Path(zarr_path) / kind / key).is_dir():
            return f"{kind}/{key}"
    return None

# Return column names of a table in spatialdata.
def obs_columns(
    zarr_path: str | Path,
    key: str,
) -> list[str]:
    """
    Column names of table ``key``'s obs, read from the zarr attrs (no data loaded); [] if absent.
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    key (str): Table name.
    """
    import json

    meta = Path(zarr_path) / "tables" / key / "obs" / "zarr.json"
    if not meta.is_file():
        return []
    return list(json.loads(meta.read_text()).get("attributes", {}).get("column-order", []))


def _file_exists(path: Path) -> bool:
    # a directory counts only if it has contents
    return path.is_file() or (path.is_dir() and any(path.iterdir()))

# Main function to check if a command can write its outputs without overwriting existing files.
def is_writeable(
    zarr_path: str | Path | None = None,
    elements: Iterable[str] = (),
    files: Iterable[str | Path] = (),
    columns: dict[str, Iterable[str]] | None = None,
    overwrite: bool = False,
) -> bool:
    """
    Check before any compute that a command may write its outputs.
    
    Raises FileExistsError for all outputs that exist, unless ``overwrite``.
    Parameters:
    zarr_path (str | Path, optional): The zarr store (needed if ``elements`` given).
    elements (Iterable[str]): Element names the command writes to the store.
    files (Iterable[str | Path]): Files or directories the command writes outside the store.
    columns (dict[str, Iterable[str]], optional): {table: obs columns the command adds to that existing table}.
    overwrite (bool): The command's --overwrite flag.
    """
    elements = list(elements)
    files = [Path(f) for f in files]
    if (elements or columns) and zarr_path is None:
        raise ValueError("is_writeable: elements / columns given without a zarr store")

    found = [p for k in elements if (p := element_on_disk(zarr_path, k))] # existing elements
    found += [str(f) for f in files if _file_exists(f)] # existing files
    for key, cols in (columns or {}).items(): # existing columns in existing columns
        have = set(obs_columns(zarr_path, key))
        if hit := [c for c in cols if c in have]:
            found.append(f"tables/{key} obs: {', '.join(hit)}")
    # check if any of the outputs exist, and raise if so (unless overwrite)
    if found and not overwrite:
        listing = "\n  ".join(found)
        raise FileExistsError(f"already exists (pass --overwrite to replace):\n  {listing}")
    if found:
        logger.info("Overwriting: %s", ", ".join(found))
    return True


def require(
    zarr_path: str | Path | None = None,
    elements: Iterable[str] = (),
    files: Iterable[str | Path] = (),
    columns: dict[str, Iterable[str]] | None = None,
) -> bool:
    """
    Check before any compute that a command's inputs exist.

    Raises FileNotFoundError listing every missing input.
    Parameters:
    zarr_path (str | Path, optional): The zarr store (needed if ``elements`` / ``columns`` given).
    elements (Iterable[str]): Element names the command reads from the store.
    files (Iterable[str | Path]): Files or directories (non-empty) the command reads outside the store.
    columns (dict[str, Iterable[str]], optional): {table: obs columns the command reads from that table}.
    """
    elements = list(elements)
    files = [Path(f) for f in files]
    if (elements or columns) and zarr_path is None:
        raise ValueError("require: elements / columns given without a zarr store")

    missing = [k for k in elements if not element_on_disk(zarr_path, k)]  # missing elements
    missing += [str(f) for f in files if not _file_exists(f)]  # missing files
    for key, cols in (columns or {}).items():  # missing columns (or their whole table)
        if element_on_disk(zarr_path, key) != f"tables/{key}":
            missing.append(f"tables/{key}")
        elif miss := [c for c in cols if c not in set(obs_columns(zarr_path, key))]:
            missing.append(f"tables/{key} obs: {', '.join(miss)}")
    if missing:
        listing = "\n  ".join(dict.fromkeys(missing))
        raise FileNotFoundError(f"missing inputs:\n  {listing}")
    return True

# save spatialdata import. 
def _sd():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd
    return sd


def open_store(
    zarr_path: str | Path,
):
    """
    The region's whole store as a SpatialData object (elements read lazily where spatialdata does).
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    """
    sd = _sd()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return sd.read_zarr(zarr_path)


def read_table(
    zarr_path: str | Path,
    key: str,
):
    """
    One table of the store, read with anndata only (same AnnData as ``open_store(...)[key]``).
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    key (str): Table name.
    """
    import anndata as ad

    logger.info("Loading %s from %s", key, zarr_path)
    return ad.read_zarr(Path(zarr_path) / "tables" / key)


def _drop(sdata, key: str) -> None:
    # remove from memory first: a lazy element may be backed by the files being deleted
    if key in sdata:
        del sdata[key]
    sdata.delete_element_from_disk(key)

# Safe write element to spatialdata zarr store.
def write_element(
    zarr_path: str | Path,
    key: str,
    element,
    overwrite: bool = False,
) -> None:
    """
    Write ``element`` to the store as ``key``.
    An existing ``key`` is replaced only with ``overwrite``: written to ``key__tmp`` first, then
    swapped in, so a crash never leaves the store without a complete copy.
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    key (str): Element name.
    element: A parsed spatialdata element (image, labels, points, shapes or table).
    overwrite (bool): Replace ``key`` if it already exists.
    """
    sd = _sd()
    zarr_path = Path(zarr_path)
    tmp = f"{key}{TMP_SUFFIX}"
    sdata = sd.read_zarr(zarr_path)

    if element_on_disk(zarr_path, tmp):  # left by a crashed swap
        logger.warning("Removing %s left by an earlier failed write.", tmp)
        _drop(sdata, tmp)

    if not element_on_disk(zarr_path, key):
        sdata[key] = element
        sdata.write_element(key)
        return
    if not overwrite:
        raise FileExistsError(f"{element_on_disk(zarr_path, key)} already exists (pass --overwrite to replace)")

    # 1. write the new element under tmp
    sdata[tmp] = element
    sdata.write_element(tmp)
    # 2. replace key with the copy on disk
    sdata = sd.read_zarr(zarr_path)
    _drop(sdata, key)
    sdata[key] = sdata[tmp]
    sdata.write_element(key)
    # 3. remove tmp
    sdata = sd.read_zarr(zarr_path)
    _drop(sdata, tmp)


# For doing a faster update on an anndata object without needing to write the whole table. 
def update_obs(
    zarr_path: str | Path,
    key: str,
    obs,
) -> None:
    """
    Replace only the obs of table ``key``; X, layers, obsm and uns are left on disk untouched.
    Not crash-safe: a failure mid-write leaves obs (and spatialdata reads of the store) broken
    until a complete obs is written back (tested; see docs/io_redesign.md).
    Parameters:
    zarr_path (str | Path): The region's zarr store.
    key (str): Table name.
    obs (pd.DataFrame): The new obs; must have the table's index, in the same order.
    """
    import zarr
    from anndata.io import read_elem, write_elem

    zarr_path = Path(zarr_path)
    if element_on_disk(zarr_path, key) != f"tables/{key}":
        raise KeyError(f"no table {key} in {zarr_path}")
    group = zarr.open_group(str(zarr_path / "tables" / key), mode="a")
    old = read_elem(group["obs"])
    if not obs.index.equals(old.index):
        raise ValueError(f"update_obs: new obs index does not match table {key}")
    uns = group["uns"]
    attrs = read_elem(uns["spatialdata_attrs"]) if "spatialdata_attrs" in uns else {}
    missing = [attrs[k] for k in ("region_key", "instance_key") if k in attrs and attrs[k] not in obs.columns]
    if missing:
        raise ValueError(f"update_obs: obs lacks {missing}, which spatialdata needs to link the table")
    write_elem(group, "obs", obs)
    with warnings.catch_warnings():  # zarr notes the parquet elements it does not own
        warnings.filterwarnings("ignore", module="zarr")
        zarr.consolidate_metadata(str(zarr_path))
