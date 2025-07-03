import os
from pathlib import Path
import glob
import warnings
import anndata as ad


def backup_adata_region(exp_name:str, reg_name:str, prefix_name:str, adata_path:str=None):
    """
    Backup function for AnnData objects. Due to env incompatibilities with 
    SpatialData and some of the annotation / integration libraries. 

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd
        from spida._utilities import _backup_adata, _gen_keys, _region_to_donor
        from spida._constants import TABLE_KEY 

    print("BACKING UP ADATA, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )


    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    _donor = _region_to_donor(reg_name)
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # Getting the sdata object (right now from a constant zarr store path)
        zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
        zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
        sdata = sd.read_zarr(zarr_path)

    if not adata_path:
        adata_path = os.getenv("ANNDATA_STORE_PATH")
        if not adata_path:
            raise ValueError("Please provide an adata_path or set the ANNDATA_STORE_PATH environment variable.")

    # getting the adata object
    adata = ad.read_h5ad(f"{adata_path}/{exp_name}/{prefix_name}/adata_{_donor}.h5ad")

    _backup_adata(sdata, adata, KEYS[TABLE_KEY])

def backup_adata_experiment(exp_name:str, prefix_name:str, adata_path:str=None):
    """
    Backup function for AnnData objects for an entire experiment.
    
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    """
    
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        backup_adata_region(exp_name, reg_name, prefix_name, adata_path=adata_path)