#!/usr/bin/env python3
"""
CLI interface for SPIDA P module (PreProcessing) using argparse.
"""
import sys
import argparse
import logging
import inspect
import warnings

from spida.utilities.script_utils import ParseKwargs, parse_path, parse_list
import spida.P as P

logger = logging.getLogger(__package__)
warnings.filterwarnings('ignore', category=UserWarning, module='zarr')

DESCRIPTION = """
    Command line interface for preprocessing pipeline steps of spatial data
    This interface provides methods to filter spatial data, setup an AnnData object for downstream analysis, 
    run spatially aware clustering and alignment algorithms, and visualize the results.
    
    Methods:

    [Utilities] 
    write-adata                        - Write AnnData objects to a .h5ad file for a specific experiment and region.
    
    [Plotting] 
    plot-filtering-region              - Plot function for filtering results for a specific region.
    plot-filtering-experiment          - Plot function for filtering results for an entire experiment. # NOT IMPLEMENTED
    plot-setup-region                  - Plot function for setup results for a specific region.
    plot-setup-experiment              - Plot function for setup results for an entire experiment.# NOT IMPLEMENTED

    [Filtering]
    filter-cells-region               - Filter cells for a specific region in an experiment.
    filter-cells-all                  - Filter cells for all regions in an experiment.
    
    [Setup]
    setup-adata-region                - Setup the AnnData object for downstream analysis for a specific region.
    setup-adata-all                   - Setup the AnnData objects for all regions in an experiment.
    remove-doublets-region            - Remove doublets from the AnnData object. (gpu required)
    remove-doublets-all               - Remove doublets from all AnnData objects in an experiment. (gpu required)
    resolvi-region                    - Run the Resolvi algorithm for a specific region. (gpu required)
    resolvi-all                       - Run the Resolvi algorithm for all regions in an experiment. (gpu required)
    """
EPILOGUE = """
Author: Amit Klein
Documentation:
"""

def config_warnings(): 
    """
    Turning off warnings for some internal libraries that are not relevant for the user.
    """
    warnings.filterwarnings('ignore', category=UserWarning, module='zarr')
    warnings.filterwarnings('ignore', category=UserWarning, module='anndata')
    warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')
    warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
    warnings.filterwarnings('ignore', category=UserWarning, module='xarray_schema')
    warnings.filterwarnings('ignore', category=FutureWarning, module='dask')
    warnings.filterwarnings('ignore', category=UserWarning, module='ome_zarr')
    warnings.filterwarnings('ignore', category=SyntaxWarning, module='leidenalg')

def setup_logging(stdout=False, quiet=False, **kwargs): 

    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setLevel(logging.INFO if not quiet else logging.WARNING)
    stream_handler.setFormatter(logging.Formatter(
        '[%(levelname)s | %(name)s | %(module)s | L%(lineno)d] %(asctime)s - %(message)s',
        datefmt="%Y-%m-%dT%H:%M:%S%z"
    ))
    
    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO if not quiet else logging.WARNING)


# Subparser registration functions
def filter_cells_region_register_subparser(subparser):
    """Register subparser for filter-cells-region command."""
    parser = subparser.add_parser(
        'filter-cells-region',
        aliases=['filter_cells_region'],
        help='Filter cells for a specific region in an experiment',
        description='Filter cells for a specific region in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--cutoffs_path', type=parse_path, help='Path to the cutoffs JSON file')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    return


def filter_cells_all_register_subparser(subparser):
    """Register subparser for filter-cells-all command."""
    parser = subparser.add_parser(
        'filter-cells-all',
        aliases=['filter_cells_all'],
        help='Filter cells for all regions in an experiment',
        description='Filter cells for all regions in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--cutoffs_path', type=parse_path, help='Path to the cutoffs JSON file')
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    return


def write_adata_register_subparser(subparser):
    """Register subparser for write-adata command."""
    parser = subparser.add_parser(
        'write-adata',
        aliases=['write_adata'],
        help='Write AnnData objects to disk',
        description='Write AnnData objects to disk for a specific experiment and region'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('--reg_name', type=str, help='Name of the region (optional, if not provided writes all regions)')
    parser.add_argument('--prefix_names', type=parse_list, required=True, help='List of prefix names')
    parser.add_argument('--output_path', type=parse_path, default=None, 
                       help='Output directory path (default is None)')
    return

def setup_adata_region_register_subparser(subparser):
    """Register subparser for setup-adata-region command."""
    parser = subparser.add_parser(
        'setup-adata-region',
        aliases=['setup_adata_region'],
        help='Setup the AnnData object for downstream analysis for a specific region',
        description='Setup the AnnData object for downstream analysis for a specific region'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--suffix', type=str, default="_filt", help="Suffix for the output table key in the spatialdata object (default: '_filt')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    return


def setup_adata_all_register_subparser(subparser):
    """Register subparser for setup-adata-all command."""
    parser = subparser.add_parser(
        'setup-adata-all',
        aliases=['setup_adata_all'],
        help='Setup the AnnData objects for all regions in an experiment',
        description='Setup the AnnData objects for all regions in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--suffix', type=str, default="_filt", help="Suffix for the output table key in the spatialdata object (default: '_filt')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    return

def resolvi_region_register_subparser(subparser):
    """Register subparser for resolvi-region command."""
    parser = subparser.add_parser(
        'resolvi-region',
        aliases=['resolvi_region'],
        help='Run the Resolvi algorithm for a specific region',
        description='Run the Resolvi algorithm for a specific region'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--suffix', type=str, default="_filt", help="Suffix for the table key in the spatialdata object (default: '_filt')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    parser.add_argument('--model_kwargs', nargs='*', action=ParseKwargs, default={},
                        help='Additional keyword arguments for the Resolvi model (default: {})')
    return


def resolvi_all_register_subparser(subparser):
    """Register subparser for resolvi-all command."""
    parser = subparser.add_parser(
        'resolvi-all',
        aliases=['resolvi_all'],
        help='Run the Resolvi algorithm for all regions in an experiment',
        description='Run the Resolvi algorithm for all regions in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--suffix', type=str, default="_filt", help="Suffix for the table key in the spatialdata object (default: '_filt')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    parser.add_argument('-k', '--model_kwargs', nargs='*', action=ParseKwargs, default={},
                        help='Additional keyword arguments for the Resolvi model (default: {})')
    return

def remove_doublets_region_register_subparser(subparser):
    """Register subparser for remove-doublets-region command."""
    parser = subparser.add_parser(
        'remove-doublets-region',
        aliases=['remove_doublets_region'],
        help='Remove doublets from the AnnData object for a specific region',
        description='Remove doublets from the AnnData object for a specific region'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--threshold', type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
    parser.add_argument('--suffix', type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    parser.add_argument('--model_kwargs', nargs='*', action=ParseKwargs, default={},
                        help='Additional keyword arguments for the SOLO model (default: {})')
    return

def remove_doublets_all_register_subparser(subparser):
    """Register subparser for remove-doublets-all command."""
    parser = subparser.add_parser(
        'remove-doublets-all',
        aliases=['remove_doublets_all'],
        help='Remove doublets from all AnnData objects in an experiment',
        description='Remove doublets from all AnnData objects in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--threshold', type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
    parser.add_argument('--suffix', type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
    parser.add_argument('--plot', action='store_true', help='Whether to plot the results')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    parser.add_argument('--model_kwargs', nargs='*', action=ParseKwargs, default={},
                        help='Additional keyword arguments for the SOLO model (default: {})')
    return

def plot_filtering_region_register_subparser(subparser):
    """Register subparser for plot-filtering-region command."""
    parser = subparser.add_parser(
        'plot-filtering-region',
        aliases=['plot_filtering_region'],
        help='Plot the filtering results for a specific region',
        description='Plot the filtering results for a specific region in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
    return


def plot_setup_region_register_subparser(subparser):
    """Register subparser for plot-setup-region command."""
    parser = subparser.add_parser(
        'plot-setup-region',
        aliases=['plot_setup_region'],
        help='Plot the setup results for a specific region',
        description='Plot the setup results for a specific region in an experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--image_path', type=parse_path, help='Path to save the plot')
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

    config_warnings()

    # collect args
    args_var = vars(args)
    for key, value in args_var.items():
        logger.info(f"Argument: {key} = {value}, Type: {type(value)}")

    command = args_var.pop("command").lower().replace("_", "-")
    if command in ["filter-cells-region"]: 
        from .main import filter_cells_region as func
    elif command in ["filter-cells-all"]:
        from .main import filter_cells_all as func
    elif command in ["write-adata"]:
        from .main import write_adata as func
    elif command in ["setup-adata-region"]:
        from .main import setup_adata_region as func
    elif command in ["setup-adata-all"]:
        from .main import setup_adata_all as func
    elif command in ["resolvi-region"]:
        from .main import resolvi_cluster_region as func
    elif command in ["resolvi-all"]:
        from .main import resolvi_cluster_all as func
    elif command in ["remove-doublets-region"]:
        from .main import remove_doublets_region as func
    elif command in ["remove-doublets-all"]:
        from .main import remove_doublets_all as func
    elif command in ["plot-filtering-region"]:
        from .main import plot_filtering_region as func
    elif command in ["plot-setup-region"]:
        from .main import plot_setup_region as func
    else:
        logger.error(f"Unknown command: {command}")
        parser.print_help()
        sys.exit(1)

    # validate environment (if needed)

    logger.info(f"Running command: {command}")
    func(**args_var)
    logger.info(f"Command {command} completed successfully.")
    return

if __name__ == "__main__":
    main()
