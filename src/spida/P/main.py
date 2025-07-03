import os
import sys
import glob
from pathlib import Path
from dotenv import load_dotenv # type: ignore
load_dotenv()
import fire # type: ignore
import warnings


warnings.filterwarnings('ignore', category=UserWarning, module='zarr')

from spida.P.filtering import run_filtering
from spida.P.setup_adata import run_setup
from spida.pl import plot_filtering, plot_setup
from spida._utilities import _gen_keys, _region_to_donor, _write_adata, _backup_adata, _get_adata
from spida._constants import TABLE_KEY # type: ignore

from matplotlib.backends.backend_pdf import PdfPages

# Setting Logging for this module 
import logging
logging.basicConfig(filename="/ceph/cephatlas/aklein/spida/tests/spida.log", level=logging.INFO)

class P_cli(): 
    """
    Command line interface for filtering and setting up AnnData objects in the BICAN project.
    This class provides methods to filter cells in specific regions, setup AnnData objects, and write them to disk.
    
    Methods:
    - filter_cells_region: Filter cells for a specific region in an experiment.
    - filter_cells_all: Filter cells for all regions in an experiment.
    - write_adata: Write AnnData objects to disk for a specific experiment and region.
    - setup_adata_region: Setup the AnnData object for downstream analysis for a specific region.
    - setup_adata_all: Setup the AnnData objects for all regions in an experiment.
    """

    
    def filter_cells_region(self, 
                            exp_name:str,
                            reg_name:str, 
                            prefix_name:str, 
                            cutoffs_path:Path=None, 
                            plot:bool=False, 
                            image_path:Path=None): 

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

        # default cutoffs path
        if cutoffs_path is None: 
            cutoffs_path = Path("/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json")
        
        # determining donor from region name
        donor_name = _region_to_donor(reg_name)
        
        # Get KEYS
        KEYS = _gen_keys(prefix_name, exp_name, reg_name)
        # Get the AnnData object
        adata = _get_adata(exp_name, reg_name, prefix_name)
        # Run the filtering 
        adata = run_filtering(adata,exp_name,reg_name,prefix_name,donor_name,cutoffs_path)
        # backup the AnnData object
        _backup_adata(exp_name, reg_name, adata, KEYS[TABLE_KEY])
        
        print("DONE")
        if plot:
            self.plot_filtering_region(exp_name, reg_name, prefix_name, image_path=image_path)

    def plot_filtering_region(self, exp_name:str, reg_name:str, prefix_name:str, image_path:Path=None):
        """
        Plot the filtering results for a specific region in an experiment.
        Parameters:
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        image_path (Path, optional): Path to save the plot. If None, uses a default path.
        """

        adata = _get_adata(exp_name, reg_name, prefix_name)

        if image_path is None: 
            image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_filt.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)

        pdf_file = PdfPages(image_path)
        fig, ax = plot_filtering(adata, exp_name, reg_name, prefix_name)
        pdf_file.savefig(fig)
        pdf_file.close()

    def filter_cells_all(self, 
                         exp_name:str, 
                         prefix_name:str, 
                         cutoffs_path:Path=None, 
                         plot:bool=False, 
                         image_path:Path=None): 
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
            self.filter_cells_region(exp_name, reg_name, prefix_name, cutoffs_path=cutoffs_path, plot=plot, image_path=image_path)


    def write_adata(self, 
                    exp_name,
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

    def setup_adata_region(self, 
                           exp_name:str, 
                           reg_name:str, 
                           prefix_name:str,
                           plot=False, 
                           image_path:Path=None):
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

        # determining donor from region name
        donor_name = _region_to_donor(reg_name)
        # Get KEYS
        KEYS = _gen_keys(prefix_name, exp_name, reg_name)
        # Get the AnnData object
        adata = _get_adata(exp_name, reg_name, prefix_name)
        # Run the setup
        adata = run_setup(adata, exp_name, reg_name, prefix_name, donor_name)
        # backup the adata object
        _backup_adata(exp_name, reg_name, adata, KEYS[TABLE_KEY])
        
        print("DONE SETUP")
        if plot: 
            self.plot_setup_region(exp_name, reg_name, prefix_name, image_path=image_path)

    def plot_setup_region(self, exp_name:str, reg_name:str, prefix_name:str, image_path:Path=None): 
        """
        Plot the setup results for a specific region in an experiment.
        This function generates a PDF with the setup plots.
        
        Parameters: 
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        image_path (Path, optional): Path to save the plot. If None, uses a default path.
        """

        adata = _get_adata(exp_name, reg_name, prefix_name)

        if image_path is None: 
            image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_setup.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)

        pdf_file = PdfPages(image_path)
        plot_setup(adata, exp_name, reg_name, prefix_name, pdf_file=pdf_file)
        pdf_file.close()


    def setup_adata_all(self, 
                        exp_name:str, 
                        prefix_name:str,
                        plot:bool=False,
                        image_path:Path=None): 
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
            self.setup_adata_region(exp_name,
                                    reg_name, 
                                    prefix_name,
                                    plot=plot,
                                    image_path=image_path)


if __name__ == "__main__": 
    fire.Fire(P_cli)