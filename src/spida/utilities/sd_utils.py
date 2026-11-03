
from pathlib import Path
import logging

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


def _write_adata(zarr_path, table_key: str, output_file: Path):
    """
    Write one table of the store at ``zarr_path`` to an H5AD file.
    Parameters:
    zarr_path (str | Path): The dataset's zarr store.
    table_key (str): The table to export.
    output_file (Path): The H5AD file to write.
    """
    from spida.io.store import read_table

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    adata = read_table(zarr_path, table_key)
    adata.write_h5ad(Path(output_file))


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