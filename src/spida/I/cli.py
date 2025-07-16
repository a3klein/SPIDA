#!/usr/bin/env python3
"""
CLI interface for SPIDA integration and annotation pipeline using argparse.
"""
import os
import sys
import argparse
import logging
from pathlib import Path
import inspect

from spida.utilities.script_utils import ParseKwargs, parse_list, parse_path

# logger = logging.getLogger(__name__)
import spida.I as I

DESCRIPTION = """
    Command line interface for integration and annotation tasks in the BICAN project.
    This interface provides methods to backup AnnData objects, run ALLCools integration, 
    setup MapMyCells integration, and run MapMyCells annotation for specific regions or entire experiments.
    
    Methods:

    [Utilities] 
    backup-adata-region                - Backup function for AnnData objects for a specific region.
    backup-adata-experiment            - Backup function for AnnData objects for an entire experiment.
    
    [ALLCools Integration]
    allcools-integration-region        - Run ALLCools integration on a given experiment and region.
    allcools-integration-experiment    - Run ALLCools integration for an entire experiment.
    
    [MapMyCells Annotation]
    mmc-setup                          - Setup function for MapMyCells integration.
    mmc-annotation-region              - Run MapMyCells annotation on a given experiment and region.
    mmc-annotation-experiment          - Run MapMyCells annotation for an entire experiment.
    
    [MOSCOT Integration]
    moscot-integration                 - Placeholder for MOSCOT integration (not implemented yet).
    """
EPILOGUE = """
Author: Amit Klein
Documentation:
"""

def setup_logging(**kwargs): 
    
    logging.basicConfig(level=logging.INFO,
                       format='[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s - %(message)s',
                       datefmt="%Y-%m-%dT%H:%M:%S%z")


def backup_adata_region_register_subparser(subparser):
    """Register subparser for backup_adata_region command."""
    parser = subparser.add_parser(
        'backup-adata-region',
        help='Backup AnnData objects for a specific region',
        description='Backup function for AnnData objects for a specific region due to env incompatibilities'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--adata-path', type=parse_path, default=None, help='Path to the store of AnnData objects (default: None)')
    return

def backup_adata_experiment_register_subparser(subparser):
    """Register subparser for backup_adata_experiment command."""
    parser = subparser.add_parser(
        'backup-adata-experiment',
        help='Backup AnnData objects for an entire experiment',
        description='Backup function for AnnData objects for an entire experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('--adata-path', type=parse_path, default=None, help='Path to the store of AnnData objects (default: None)')
    return

def allcools_integration_region_register_subparser(subparser):
    """Register subparser for allcools_integration_region command."""
    parser = subparser.add_parser(
        'allcools-integration-region',
        help='Run ALLCools integration on a given experiment and region',
        description='Run ALLCools integration on a given experiment and region. ENVIRONMENT = "preprocessing"'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('ref_path', type=parse_path, help='Path to the reference RNA AnnData object')
    parser.add_argument('--anndata-store-path', type=parse_path, default=None, help='Path to the store of AnnData objects (default: None)')
    parser.add_argument('--annotations-store-path', type=parse_path, default=None, help='Path to the store of annotation specific files (default: None)')
    parser.add_argument('-k', '--kwargs', nargs='*', action=ParseKwargs,
                        help='Additional keyword arguments for ALLCools integration in key=value format')
    return

def allcools_integration_experiment_register_subparser(subparser):
    """Register subparser for allcools_integration_experiment command."""
    parser = subparser.add_parser(
        'allcools-integration-experiment',
        help='Run ALLCools integration for an entire experiment',
        description='Run ALLCools integration for an entire experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('ref_path', type=str, help='Path to the reference RNA AnnData object')
    parser.add_argument('--anndata-store-path', type=parse_path, default=None, help='Path to the store of AnnData objects (default: None)')
    parser.add_argument('--annotations-store-path', type=parse_path, default=None, help='Path to the store of annotation specific files (default: None)')
    parser.add_argument('-k', '--kwargs', nargs='*', action=ParseKwargs,
                        help='Additional keyword arguments for ALLCools integration in key=value format')
    return

def mmc_setup_register_subparser(subparser):
    """Register subparser for mmc_setup command."""
    parser = subparser.add_parser(
        'mmc-setup',
        help='Setup function for MapMyCells integration',
        description='Setup function for MapMyCells integration'
    )
    parser.add_argument('ref_path', type=parse_path, help='Path to the reference data')
    parser.add_argument('hierarchy_list', type=parse_list, help='List of hierarchy levels for the annotation')
    parser.add_argument('brain_region', type=str, help='Brain region for the annotation')
    parser.add_argument('codebook', type=str, help='Codebook for the annotation')
    parser.add_argument('--codebook-path', type=parse_path, default=None, help='Path to the codebook file (default: None)')
    parser.add_argument('--mmc-store-path', type=parse_path, default=None, help='Path to the store of MapMyCells markers (default: None)')
    parser.add_argument('--ref-norm', type=str, default='log2CPM', help='Normalization method for the reference AnnData.X (default: log2CPM)')
    parser.add_argument('-k', '--kwargs', nargs='*', action=ParseKwargs,
                        help='Additional keyword arguments in key=value format')
    return

def mmc_annotation_region_register_subparser(subparser):
    """Register subparser for mmc_annotation_region command."""
    parser = subparser.add_parser(
        'mmc-annotation-region',
        help='Run MapMyCells annotation on a given experiment and region',
        description='Run MapMyCells annotation on a given experiment and region'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('reg_name', type=str, help='Name of the region')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('brain_region', type=str, help='Brain region for the annotation')
    parser.add_argument('codebook', type=str, help='Codebook for the annotation')
    parser.add_argument('--mmc-store-path', type=parse_path, default=None, help='Path to the store of MapMyCells markers (default: None)')
    parser.add_argument('--anndata-store-path', type=parse_path, default=None, help='Path to the store of AnnData objects (default: None)')
    parser.add_argument('--annotations-store-path', type=parse_path, default=None, help='Path to the store of annotation specific files (default: None)')
    parser.add_argument('-k', '--kwargs', nargs='*', action=ParseKwargs,
                        help='Additional keyword arguments in key=value format')
    return

def mmc_annotation_experiment_register_subparser(subparser):
    """Register subparser for mmc_annotation_experiment command."""
    parser = subparser.add_parser(
        'mmc-annotation-experiment',
        help='Run MapMyCells annotation for an entire experiment',
        description='Run MapMyCells annotation for an entire experiment'
    )
    parser.add_argument('exp_name', type=str, help='Name of the experiment')
    parser.add_argument('prefix_name', type=str, help='Prefix for the keys in the spatialdata object')
    parser.add_argument('brain_region', type=str, help='Brain region for the annotation')
    parser.add_argument('codebook', type=str, help='Codebook for the annotation')
    parser.add_argument('--mmc-store-path', type=str, default=None, help='Path to the store of MapMyCells markers (default: None)')
    parser.add_argument('--anndata-store-path', type=str, default=None, help='Path to the store of AnnData objects (default: None)')
    parser.add_argument('--annotations-store-path', type=str, default=None, help='Path to the store of annotation specific files (default: None)')
    parser.add_argument('-k', '--kwargs', nargs='*', action=ParseKwargs,
                        help='Additional keyword arguments in key=value format')
    return

def moscot_integration_register_subparser(subparser):
    """Register subparser for moscot_integration command."""
    parser = subparser.add_parser(
        'moscot-integration',
        help='Placeholder for MOSCOT integration (not implemented yet)',
        description='Placeholder for MOSCOT integration (not implemented yet)'
    )
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
        metavar='{backup-adata-region,backup-adata-experiment,allcools-integration-region,allcools-integration-experiment,mmc-setup,mmc-annotation-region,mmc-annotation-experiment,moscot-integration}'
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
            print(I.__version__)
            exit()
        else: 
            args = parser.parse_args()
    else:
        args = parser.parse_args(['--help'])
        exit()
    
    # setup logging here: 
    if not logging.root.handlers: 
        setup_logging(stdout=True, quiet=False)

    # collect args
    args_var = vars(args)
    for key, value in args_var.items():
        logging.info(f"Argument: {key} = {value}, Type: {type(value)}")

    command = args_var.pop("command").lower().replace("_", "-")
    if command in ["backup-adata-region"]:
        from main import backup_adata_region as func
    elif command in ["backup-adata-experiment"]:
        from main import backup_adata_experiment as func
    elif command in ["allcools-integration-region"]:
        from main import allcools_integration_region as func
    elif command in ["allcools-integration-experiment"]:
        from main import allcools_integration_experiment as func
    elif command in ["mmc-setup"]:
        from main import mmc_setup as func
    elif command in ["mmc-annotation-region"]:
        from main import mmc_annotation_region as func
    elif command in ["mmc-annotation-experiment"]:
        from main import mmc_annotation_experiment as func
    elif command in ["moscot-integration"]:
        from main import moscot_integration as func
    else:
        logging.error(f"Unknown command: {command}")
        parser.print_help()
        sys.exit(1)

    # validate environment (if needed)

    logging.info(f"Running command: {command}")
    func(**args_var)
    logging.info(f"Command {command} completed successfully.")
    return 

if __name__ == "__main__":
    main()
