
import warnings
from pathlib import Path
import logging

import anndata as ad

from spida._constants import *

# Set up logging
logger = logging.getLogger(__name__)


def _gen_keys(prefix_name, exp_name, reg_name):
    """
    Generate keys for spatialdata objects based on the experiment and region names.

    Parameters:
    prefix_name (str): Prefix for the keys.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.

    Returns:
    dict: A dictionary containing the generated keys.
    """

    keys = {
        SHAPES_KEY: f"{prefix_name}_{exp_name}_{reg_name}_polygons",
        POINTS_KEY: f"{prefix_name}_{exp_name}_{reg_name}_transcripts",
        TABLE_KEY: f"{prefix_name}_table",
        IMAGE_KEY: f"default_{exp_name}_{reg_name}_z3",
    }

    return keys


def _region_to_donor(reg_name: str) -> str:
    """
    Convert a region name to a donor name.

    Parameters:
    reg_name (str): Name of the region.

    Returns:
    str: Donor name derived from the region name.
    """
    if "2424" in reg_name:
        return "UCI2424"
    elif "5224" in reg_name:
        return "UCI5224"
    elif "4723" in reg_name:
        return "UCI4723"
    elif "7648" in reg_name:
        return "UWA7648"
    else:
        logger.warning(
            f"Region name {reg_name} does not match any known donor names. Returning the region name as donor name."
        )
        return reg_name

def _validate_adata(adata):
    """
    Making sure that an anndata does not contain any case insensitive duplicates in obs_name
    """

    # From Copilot
    def remove_case_insensitive_duplicates(df):
        """
        Remove case-insensitive duplicate columns from a DataFrame, keeping uppercase versions.

        Parameters:
            df (pd.DataFrame): Input DataFrame

        Returns:
            pd.DataFrame: DataFrame with case-insensitive duplicates removed
        """
        # Create a mapping of lowercase column names to their original names
        col_mapping = {}
        for col in df.columns:
            col_lower = col.lower()
            if col_lower in col_mapping:
                # If current column name is uppercase and stored one is lowercase, replace it
                if col.isupper() and not col_mapping[col_lower].isupper():
                    col_mapping[col_lower] = col
            else:
                col_mapping[col_lower] = col

        # Get the list of columns to keep
        cols_to_keep = list(col_mapping.values())

        # Return DataFrame with only the selected columns
        return df[cols_to_keep]

    adata.obs = remove_case_insensitive_duplicates(adata.obs)
    return adata


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


def _write_adata(zarr_path, table_key: str, output_file: Path):
    """
    Write one table of the store at ``zarr_path`` to an H5AD file.
    Parameters:
    zarr_path (str | Path): The dataset's zarr store.
    table_key (str): The table to export.
    output_file (Path): The H5AD file to write.
    """
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_zarr(f"{zarr_path}/tables/{table_key}")
    adata.write_h5ad(Path(output_file))


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


def _get_obs_or_gene(adata, col, layer): 
    """
    Helper function to get the column from obs or var.
    """
    import pandas.api as papi
    # ensure cluster_by in obs
    if col in adata.obs.columns:  # check if col is in obs
        _drop_col = False
        if not papi.types.is_numeric_dtype(adata.obs[col].dtype): # make sure coord_base is numeric
            try: 
                adata.obs[col] = adata.obs[col].astype(float)
            except ValueError:
                raise ValueError(f"Column {col} is not numeric, please check the data type.")
    elif col in adata.var_names:  # check if col is in var_names
        _drop_col = True
        if col not in adata.var_names:
            raise ValueError(f"Color by {col} is not in adata.obs and is not a valid gene.")
        if layer is None: 
            adata.obs[col] = adata[:,col].X.toarray().flatten() 
        else: 
            adata.obs[col] = adata[:,col].layers[layer].toarray().flatten()
    else:
        raise ValueError(f"Column {col} not found in adata.obs or adata.var_names.")

    return adata, _drop_col