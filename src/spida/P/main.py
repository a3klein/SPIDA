import os
import sys
import glob
from pathlib import Path
from dotenv import load_dotenv # type: ignore
load_dotenv()
import fire # type: ignore
import warnings


warnings.filterwarnings('ignore', category=UserWarning, module='zarr')

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

sys.path.append(os.getenv("SPIDA_PATH"))
from spida.P.filtering import run_filtering
from spida.P.setup_adata import run_setup
from spida.plotting.P_plots import plot_filtering, plot_setup
from spida._utilities import _gen_keys, _region_to_donor, _write_adata, _backup_adata
from spida._constants import TABLE_KEY # type: ignore

from matplotlib.backends.backend_pdf import PdfPages

# Setting Logging for this module 
import logging
logging.basicConfig(filename="/ceph/cephatlas/aklein/spida/tests/spida.log", level=logging.INFO)


def filter_cells_region(exp_name:str, reg_name:str, prefix_name:str, cutoffs_path:Path=None, plot:bool=False, image_path:Path=None): 

    """
    Filter cells for a specific region in an experiment.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file. If None, uses a default path.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    Raises:
    NotImplementedError: If plotting is requested but not implemented.
    """

    print("FILTERING CELLS, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    
    # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    sdata = sd.read_zarr(zarr_path)
    
    # default cutoffs path
    if cutoffs_path is None: 
        cutoffs_path = Path("/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json")
    
    # determining donor from region name
    donor_name = _region_to_donor(reg_name)

    adata = sdata[KEYS[TABLE_KEY]].copy()
    adata = run_filtering(adata,exp_name,reg_name,prefix_name,donor_name,cutoffs_path)
    
    _backup_adata(sdata, adata, KEYS[TABLE_KEY])
    
    print("DONE")

    if plot: 
        if image_path is None: 
            image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_filt.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)

        pdf_file = PdfPages(image_path)
        fig, ax = plot_filtering(adata, exp_name, reg_name, prefix_name)
        pdf_file.savefig(fig)
        pdf_file.close()

def filter_cells_all(exp_name:str, prefix_name:str, cutoffs_path:Path=None, plot:bool=False, image_path:Path=None): 
    """
    Filter cells for all regions in an experiment (calls filter_cells_all).
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file. If None, uses a default path.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """


    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        filter_cells_region(exp_name, reg_name, prefix_name, cutoffs_path=cutoffs_path, plot=plot, image_path=image_path)


def write_adata(exp_name,
                reg_name:str=None,
                prefix_names:list[str]=None,
                output_path:Path="/ceph/cephatlas/aklein/bican/data/anndatas/",
                ): 
    """
    """
    print(exp_name, reg_name, prefix_names, output_path)
    # Iterating over all prefixes
    for p in prefix_names: 
        print("PREFIX:", p)
        # if a specific region is provided, write only that region
        if reg_name: 
            fout = f"{output_path}/{exp_name}/{p}"
            _write_adata(exp_name, reg_name, p, fout)
        # if no region is provided, write all regions
        else: 
            region_list = glob.glob(f"/data/aklein/bican_zarr/{exp_name}/region_*")
            for reg in region_list: 
                rname = reg.split("/")[-1]
                fout = f"{output_path}/{exp_name}/{p}"
                _write_adata(exp_name, rname, p, fout)

def setup_adata_region(exp_name:str, reg_name:str, prefix_name:str, plot=False, image_path:Path=None):
    """
    Setup the AnnData object for downstream analysis.
    This involves normalizing data, calculating PCA, umap, tsne, and leiden clusters 

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """

    print("SETTING UP ADATA, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )

    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    sdata = sd.read_zarr(zarr_path)
    adata = sdata[KEYS[TABLE_KEY]].copy()

    # determining donor from region name
    donor_name = _region_to_donor(reg_name)

    adata = run_setup(adata, exp_name, reg_name, prefix_name, donor_name)

    # backup the adata object
    _backup_adata(sdata, adata, KEYS[TABLE_KEY])

    if plot: 
        if image_path is None: 
            image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_setup.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)

        pdf_file = PdfPages(image_path)
        plot_setup(adata, exp_name, reg_name, prefix_name, pdf_file=pdf_file)
        pdf_file.close()
    

def setup_adata_all(exp_name:str, prefix_name:str, plot:bool=False, image_path:Path=None): 
    """
    Setup the AnnData objects for all regions in an experiment.
    This function iterates over all regions and calls setup_adata for each region.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        setup_adata_region(exp_name, reg_name, prefix_name, plot=plot, image_path=image_path)


if __name__ == "__main__": 
    fire.Fire({"filter_cells_region": filter_cells_region,
               "filter_cells_all": filter_cells_all,
               "write_adata": write_adata,
               "setup_adata_region" : setup_adata_region, 
               "setup_adata_all" : setup_adata_all, 
               # "backup_adata"
    })
