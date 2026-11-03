"""Archived in phase 6 (2026-10-27) from ``spida.utilities.sd_utils``: ``_backup_adata`` (replaced by
``spida.io.store.write_table``), ``_backup_element_sd`` and ``_assign_new_table`` (no callers); and ``spida.I._io.backup_table`` (the
removed ``--backup_to_spatialdata`` path; sync-annotations replaces it)."""
import logging
import warnings

import anndata as ad

logger = logging.getLogger(__name__)


def _backup_adata(zarr_path, element: ad.AnnData, element_name: str):
    """
    Replace a table in the store at ``zarr_path`` (temp-swap: a crash never leaves it missing).
    Parameters:
    zarr_path (str | Path): The dataset's zarr store.
    element: The new element to write to the spatialdata object.
    element_name (str): The name of the element in the spatialdata object.
    """

    ### TODO: Verify this works correctly for all element types (Points / Shapes / Tables)
    # Main use case is going to be for tables
    # It Does not work for all elements :(
    from spida.io.store import write_element

    if element.raw:
        logger.info("Removing Raw data from AnnData object before writing to disk.")
        del element.raw

    write_element(zarr_path, element_name, element, overwrite=True)


def _backup_element_sd(sdata, element, element_name: str):
    """
    Backup the current state of a spatialdata element before modifying it.
    This function deletes the existing element from disk and writes the new element to the spatialdata object.
    Parameters:
    sdata (sd.SpatialData): The spatialdata object to modify.
    element: The new element to write to the spatialdata object.
    element_name (str): The name of the element in the spatialdata object.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd

    print(f"attempting to remove {element_name} from disk")
    sdata.delete_element_from_disk(element_name)  # Remove the old table from disk
    sdata[element_name] = element
    sdata.write_element(element_name)


def _assign_new_table(zarr_path, element: ad.AnnData, element_name: str):
    """
    Write a new table into the store at ``zarr_path``; raises if it already exists.
    Parameters:
    zarr_path (str | Path): The dataset's zarr store.
    element: The new element to write to the spatialdata object.
    element_name (str): The name of the element in the spatialdata object.
    """

    ### TODO: Verify this works correctly for all element types (Points / Shapes / Tables)
    # Main use case is going to be for tables
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd

    sdata = sd.read_zarr(zarr_path)

    full_name = element_name

    if f"tables/{full_name}" not in sdata.elements_paths_on_disk():
        sdata[full_name] = element
        sdata.write_element(full_name)
    else:
        raise ValueError(
            f"Element {full_name} already exists in the spatialdata object. Skipping write operation."
        )


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
