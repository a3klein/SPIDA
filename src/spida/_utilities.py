import sys
import warnings
from pathlib import Path

import anndata as ad

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

spida_path = "/ceph/cephatlas/aklein/spida/src"
sys.path.append(spida_path)
from spida._constants import *


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
        IMAGE_KEY : f"default_{exp_name}_{reg_name}_z3",
    }
    
    return keys


def _region_to_donor(reg_name:str) -> str: 
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
        raise ValueError(f"Unknown region name {reg_name}. Cannot determine donor name.")
    
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

def _backup_adata(sdata:sd.SpatialData, element:ad.AnnData, element_name:str): 
    """
    Backup the current state of a spatialdata element before modifying it.
    This function deletes the existing element from disk and writes the new element to the spatialdata object.
    Parameters:
    sdata (sd.SpatialData): The spatialdata object to modify.
    element: The new element to write to the spatialdata object.
    element_name (str): The name of the element in the spatialdata object.
    """

    ### TODO: Verify this works correctly for all element types (Points / Shapes / Tables)
    # Main use case is going to be for tables

    sdata.delete_element_from_disk(element_name)  # Remove the old table from disk
    sdata[element_name] = element
    sdata.write_element(element_name)


def _write_adata(exp_name, reg_name, prefix_name, output_path:Path): 
    """
    Write the AnnData object to an H5AD file.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    output_path (Path): Path to the output directory where the H5AD file will be saved.
    """
    # Get the donor name for the filename
    donor_name = _region_to_donor(reg_name)
    Path(output_path).mkdir(parents=True, exist_ok=True)  # Ensure the output directory exists

    sdata = sd.read_zarr(f"/data/aklein/bican_zarr/{exp_name}/{reg_name}")
    adata = sdata[f'{prefix_name}_table']
    adata.write_h5ad(Path(f"{output_path}/adata_{donor_name}.h5ad"))


