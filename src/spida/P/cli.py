#!/usr/bin/env python3
"""
CLI interface for SPIDA P module (Pipeline) using argparse.
"""
import os
import sys
import argparse
import logging
from pathlib import Path
import inspect
import glob
import warnings

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    # dotenv not available, continue without it
    pass

warnings.filterwarnings('ignore', category=UserWarning, module='zarr')

from spida.utilities.script_utils import ParseKwargs
from spida.P.filtering import run_filtering
from spida.P.setup_adata import run_setup
from spida.pl import plot_filtering, plot_setup
from spida._utilities import _gen_keys, _region_to_donor, _write_adata, _backup_adata, _get_adata
from spida._constants import TABLE_KEY
from matplotlib.backends.backend_pdf import PdfPages

# logger = logging.getLogger(__name__)
import spida.P as P

DESCRIPTION = """ SPIDA P Module (Pipeline) CLI - Filter cells and setup AnnData objects for BICAN project """
EPILOGUE = """ """

def setup_logging(**kwargs): 
    logging.basicConfig(level=logging.INFO,
                       format='[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s - %(message)s',
                       datefmt="%Y-%m-%dT%H:%M:%S%z")


# P Module Functions
def filter_cells_region(exp_name: str, reg_name: str, prefix_name: str, 
                       cutoffs_path: str = None, plot: bool = False, 
                       image_path: str = None):
    """Filter cells for a specific region in an experiment."""
    print("FILTERING CELLS, EXPERIMENT %s, REGION %s, PREFIX %s" % (exp_name, reg_name, prefix_name))
    
    # default cutoffs path
    if cutoffs_path is None: 
        cutoffs_path = os.getenv("DEF_CUTOFFS_PATH", "/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json")
    
    # determining donor from region name
    donor_name = _region_to_donor(reg_name)
    
    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    adata = _get_adata(exp_name, reg_name, prefix_name)
    # Run the filtering 
    adata = run_filtering(adata, exp_name, reg_name, prefix_name, donor_name, cutoffs_path)
    # backup the AnnData object
    _backup_adata(exp_name, reg_name, adata, KEYS[TABLE_KEY])
    
    print("DONE")
    if plot:
        plot_filtering_region(exp_name, reg_name, prefix_name, image_path=image_path)


def plot_filtering_region(exp_name: str, reg_name: str, prefix_name: str, image_path: str = None):
    """Plot the filtering results for a specific region in an experiment."""
    adata = _get_adata(exp_name, reg_name, prefix_name)

    if image_path is None: 
        image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
        image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_filt.pdf")
        image_path.parent.mkdir(parents=True, exist_ok=True)

    pdf_file = PdfPages(image_path)
    fig, ax = plot_filtering(adata, exp_name, reg_name, prefix_name)
    pdf_file.savefig(fig)
    pdf_file.close()


def filter_cells_all(exp_name: str, prefix_name: str, cutoffs_path: str = None, 
                    plot: bool = False, image_path: str = None):
    """Filter cells for all regions in an experiment."""
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        filter_cells_region(exp_name, reg_name, prefix_name, 
                           cutoffs_path=cutoffs_path, plot=plot, image_path=image_path)


def write_adata(exp_name: str, reg_name: str = None, prefix_names: list = None, 
               output_path: str = "/ceph/cephatlas/aklein/bican/data/anndatas/"):
    """Write AnnData objects to disk for a specific experiment and region."""
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
            zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
            region_list = glob.glob(f"{zarr_store}/{exp_name}/region_*")
            for reg in region_list: 
                rname = reg.split("/")[-1]
                fout = f"{output_path}/{exp_name}/{p}"
                _write_adata(exp_name, rname, p, fout)


def setup_adata_region(exp_name: str, reg_name: str, prefix_name: str,
                      plot: bool = False, image_path: str = None):
    """Setup the AnnData object for downstream analysis for a specific region."""
    print("SETTING UP ADATA, EXPERIMENT %s, REGION %s, PREFIX %s" % (exp_name, reg_name, prefix_name))

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
        plot_setup_region(exp_name, reg_name, prefix_name, image_path=image_path)


def plot_setup_region(exp_name: str, reg_name: str, prefix_name: str, image_path: str = None):
    """Plot the setup results for a specific region in an experiment."""
    adata = _get_adata(exp_name, reg_name, prefix_name)

    if image_path is None: 
        image_store = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
        image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_setup.pdf")
        image_path.parent.mkdir(parents=True, exist_ok=True)

    pdf_file = PdfPages(image_path)
    plot_setup(adata, exp_name, reg_name, prefix_name, pdf_file=pdf_file)
    pdf_file.close()


def setup_adata_all(exp_name: str, prefix_name: str, plot: bool = False, image_path: str = None):
    """Setup the AnnData objects for all regions in an experiment."""
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        setup_adata_region(exp_name, reg_name, prefix_name, plot=plot, image_path=image_path)


# Subparser registration functions
def filter_cells_region_register_subparser(subparser):
    """Register subparser for filter-cells-region command."""
    parser = subparser.add_parser(
        'filter-cells-region',
        help='Filter cells for a specific region in an experiment',
        description='Filter cells for a specific region in an experiment'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('reg_name', help='Name of the region')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--cutoffs-path', type=str, help='Path to the cutoffs JSON file')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return


def filter_cells_all_register_subparser(subparser):
    """Register subparser for filter-cells-all command."""
    parser = subparser.add_parser(
        'filter-cells-all',
        help='Filter cells for all regions in an experiment',
        description='Filter cells for all regions in an experiment'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--cutoffs-path', type=str, help='Path to the cutoffs JSON file')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return


def write_adata_register_subparser(subparser):
    """Register subparser for write-adata command."""
    parser = subparser.add_parser(
        'write-adata',
        help='Write AnnData objects to disk',
        description='Write AnnData objects to disk for a specific experiment and region'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('--reg-name', type=str, help='Name of the region (optional, if not provided writes all regions)')
    parser.add_argument('--prefix-names', type=str, nargs='+', required=True, help='List of prefix names')
    parser.add_argument('--output-path', type=str, default='/ceph/cephatlas/aklein/bican/data/anndatas/', 
                       help='Output directory path')
    return


def setup_adata_region_register_subparser(subparser):
    """Register subparser for setup-adata-region command."""
    parser = subparser.add_parser(
        'setup-adata-region',
        help='Setup the AnnData object for downstream analysis for a specific region',
        description='Setup the AnnData object for downstream analysis for a specific region'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('reg_name', help='Name of the region')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return


def setup_adata_all_register_subparser(subparser):
    """Register subparser for setup-adata-all command."""
    parser = subparser.add_parser(
        'setup-adata-all',
        help='Setup the AnnData objects for all regions in an experiment',
        description='Setup the AnnData objects for all regions in an experiment'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return


def plot_filtering_region_register_subparser(subparser):
    """Register subparser for plot-filtering-region command."""
    parser = subparser.add_parser(
        'plot-filtering-region',
        help='Plot the filtering results for a specific region',
        description='Plot the filtering results for a specific region in an experiment'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('reg_name', help='Name of the region')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return


def plot_setup_region_register_subparser(subparser):
    """Register subparser for plot-setup-region command."""
    parser = subparser.add_parser(
        'plot-setup-region',
        help='Plot the setup results for a specific region',
        description='Plot the setup results for a specific region in an experiment'
    )
    parser.add_argument('exp_name', help='Name of the experiment')
    parser.add_argument('reg_name', help='Name of the region')
    parser.add_argument('prefix_name', help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--image-path', type=str, help='Path to save the plot')
    return 



def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOGUE,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='{filter-cells-region,filter-cells-all,write-adata,setup-adata-region,setup-adata-all,plot-filtering-region,plot-setup-region}'
    )

    # Adding all the subparsers from above
    current_module = sys.modules[__name__]
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if "register_subparser" in name:
            register_subparser_func(subparsers)

    # initiate args: 
    args = None
    if len(sys.argv) > 1: 
        if sys.argv[1] in ["-v", "--version"]: 
            print(P.__version__)
            exit()
        else: 
            args = parser.parse_args()
    else:
        args = parser.parse_args(['--help'])
        exit()
    
    # setup logging here: 
    if not logging.root.handlers: 
        setup_logging()

    # collect args
    args_var = vars(args)
    for key, value in args_var.items():
        logging.info(f"Argument: {key} = {value}, Type: {type(value)}")

    command = args_var.pop("command").lower().replace("_", "-")
    
    # Function mapping
    function_map = {
        "filter-cells-region": filter_cells_region,
        "filter-cells-all": filter_cells_all,
        "write-adata": write_adata,
        "setup-adata-region": setup_adata_region,
        "setup-adata-all": setup_adata_all,
        "plot-filtering-region": plot_filtering_region,
        "plot-setup-region": plot_setup_region,
    }
    
    if command in function_map:
        func = function_map[command]
    else:
        logging.error(f"Unknown command: {command}")
        parser.print_help()
        sys.exit(1)

    logging.info(f"Running command: {command}")
    func(**args_var)
    logging.info(f"Command {command} completed successfully.")
    return

if __name__ == "__main__":
    main()
