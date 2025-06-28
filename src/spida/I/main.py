### Integration / Annotations main file 

from pathlib import Path
import glob
import warnings

import fire # type: ignore
from dotenv import load_dotenv  # type: ignore
load_dotenv()

# import warnings
# with warnings.catch_warnings():
#     warnings.filterwarnings("ignore")
#     import spatialdata as sd
import anndata as ad

import os
import sys    
# sys.path.append(os.getenv("SPIDA_PATH"))
# from spida._utilities import _gen_keys, _region_to_donor, _write_adata, _backup_adata
# from spida._constants import TABLE_KEY 



def backup_adata(exp_name:str, reg_name:str, prefix_name:str, adata_path:str=None):
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

def backup_adata_exp(exp_name:str, prefix_name:str, adata_path:str=None):
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
        backup_adata(exp_name, reg_name, prefix_name, adata_path=adata_path)



def allcools(exp_name:str,
             reg_name:str, 
             prefix_name:str, 
             ref_path:str,
             anndata_store_path:str=None, 
             annotations_store_path:str=None,
             **kwargs):
    """
    Run ALLCools integration on a given experiment and region.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    ref_path (str): Path to the reference RNA AnnData object .
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
    to None.
    **kwargs: Additional keyword arguments for ALLCools integration.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.allcools import run_allcools_seurat

    print("RUNNING ALLCOOLS INTEGRATION, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )

    # # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")

    ref_adata = ad.read_h5ad(ref_path)

    adata = run_allcools_seurat(ref_adata, adata,
                                anndata_store_path,
                                annotations_store_path,
                                **kwargs)
    
    print("DONE WITH ALLCOOLS INTEGRATION")

def allcools_experiment(exp_name:str,
                        prefix_name:str, 
                        ref_path:str,
                        anndata_store_path:str=None, 
                        annotations_store_path:str=None,
                        **kwargs):
    """
    Run ALLCools integration for an entire experiment.
    
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    ref_path (str): Path to the reference RNA AnnData object .
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
    to None.
    **kwargs: Additional keyword arguments for ALLCools integration.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        allcools(exp_name, reg_name, prefix_name, ref_path,
                 anndata_store_path=anndata_store_path,
                 annotations_store_path=annotations_store_path,
                 **kwargs)


def mmc_setup(ref_path:str,
              heirarchy_list:list,
              BRAIN_REGION:str,
              CODEBOOK:str,
              codebook_path:str=None,
              mmc_store_path:str=None,
              ref_norm:str="log2CPM", 
              **kwargs
              ): 
    """
    Setup function for MapMyCells integration.
    
    Parameters:
    ref_path (str): Path to the reference data.
    heirarchy_list (list): List of hierarchy levels for the annotation.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    codebook_path (str, optional): Path to the codebook file. Defaults to None.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    ref_norm (str, optional): Normalization method for the reference AnnData.X. Defaults to "log2CPM".
    **kwargs: Additional keyword arguments.

    
    """
    from spida.I.mmc import setup_mmc

    setup_mmc(ref_path=ref_path,
              BRAIN_REGION=BRAIN_REGION,
              CODEBOOK=CODEBOOK,
              codebook_path=codebook_path,
              heirarchy_list=heirarchy_list,
              mmc_store_path=mmc_store_path,
              ref_norm=ref_norm,
              **kwargs)
    """
    Setup function for MapMyCells integration.
    """

    print("DONE")
    return 0


def mmc_run(
        exp_name:str, 
        reg_name:str, 
        prefix_name:str, 
        BRAIN_REGION:str, 
        CODEBOOK:str, 
        mmc_store_path:str=None, 
        anndata_store_path:str=None, 
        annotations_store_path:str=None,
        **kwargs
        ):
    """
    Run MapMyCells annotation on a given experiment and region.
    
    Parameters: 
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    
    from spida.I.mmc import run_mmc

    print("RUNNING MMC ANNOTATIONS, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )

    # # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")

    ### Need to get the adata_path from this function (it must have been written out beforehand)


    return run_mmc(
            adata,
            BRAIN_REGION,
            CODEBOOK,
            mmc_store_path,
            anndata_store_path, 
            annotations_store_path,
            kwargs=kwargs,
        )
    
def mmc_experiment(
        exp_name:str, 
        prefix_name:str, 
        BRAIN_REGION:str, 
        CODEBOOK:str, 
        mmc_store_path:str=None, 
        anndata_store_path:str=None, 
        annotations_store_path:str=None,
        **kwargs
        ):
    """
    Run MapMyCells annotation for an entire experiment.
    
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        mmc_run(exp_name, reg_name, prefix_name, 
                BRAIN_REGION=BRAIN_REGION, CODEBOOK=CODEBOOK,
                mmc_store_path=mmc_store_path,
                anndata_store_path=anndata_store_path,
                annotations_store_path=annotations_store_path,
                **kwargs)

def moscot(): 
    raise NotImplementedError("MOSCOT integration is not implemented yet.")



if __name__ == "__main__":
    fire.Fire({"backup_adata_region" : backup_adata,
               "backup_adata_experiment": backup_adata_exp,
                "allcools_integration" : allcools, 
                "allcools_integration_experiment": allcools_experiment,
               "setup_mmc": mmc_setup,
               "mmc_annotation_region": mmc_run,
               "mmc_annotation_experiment": mmc_experiment,
               "moscot_integration": moscot})