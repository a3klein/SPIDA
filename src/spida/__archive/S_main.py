#!/usr/bin/env python3
"""
CLI interface for SPIDA segmentation pipeline using argparse.
"""

import sys
import argparse
import logging
import warnings
import inspect

from spida.utilities.script_utils import ParseKwargs, parse_path, parse_list, parse_dict

import spida.S as S

logger = logging.getLogger(__package__)

DESCRIPTION = """ 
    Command line interface for SPIDA S (Spatial Processing) Module.
    This interface provides methods to ingest raw data, segment spatial data using a variety of algorithms,
    deconvolve large images for better segmentation, and align cell based segmentations to the nuclei segmentation
    they are rooted in to remove dubious artifacts. 
    
    Methods:

    [Utilities] 
    decon-image                        - Deconvolve large image files in tiles using DeconWolf algorithm.
    load-decon-images                  - Load deconvolved images into spatialdata objects.

    [Segmentations] 
    run / run-segmentation-region      - Run segmentation on a single region using specified algorithm.
    experiment / segment-experiment    - Run segmentation for all regions in an experiment using specified algorithm.

    [Alignment]    
    align / align-proseg               - Align cell-based segmentations to the nuclei segmentation.

    [io]
    ingest-region                      - Ingest a specific region of an experiment into a spatialdata object.
    ingest-all                         - Ingest all regions of an experiment into spatialdata objects.
    load-segmentation-region           - Load segmentation data for a specific region into a spatialdata object.
    load-segmentation-all              - Load segmentation data for all regions of an experiment into spatialdata objects.
"""

EPILOGUE = """
Author: Amit Klein
Documentation:
"""


def config_warnings():
    """
    Turning off warnings for some internal libraries that are not relevant for the user.
    """
    warnings.filterwarnings("ignore", category=UserWarning, module="zarr")
    warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
    warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")
    warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")
    warnings.filterwarnings("ignore", category=FutureWarning, module="dask")
    warnings.filterwarnings("ignore", category=UserWarning, module="spatialdata")


def setup_logging(stdout=False, quiet=False, **kwargs):
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setLevel(logging.INFO if not quiet else logging.WARNING)
    stream_handler.setFormatter(
        logging.Formatter(
            "[%(levelname)s | %(name)s | %(module)s | L%(lineno)d] %(asctime)s - %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S%z",
        )
    )

    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO if not quiet else logging.WARNING)

    # logger.basicConfig(level=logging.INFO,
    #                    format='[%(levelname)s | %(module)s | L%(lineno)d] %(asctime)s - %(message)s',
    #                    datefmt="%Y-%m-%dT%H:%M:%S%z")


def ingest_region_register_subparser(subparser):
    """Register subparser for ingest-region command."""
    parser = subparser.add_parser(
        "ingest-region",
        aliases=["ingest_region"],
        help="Ingest a specific region of an experiment into a spatialdata object",
        description="Ingest a specific region of an experiment into a spatialdata object",
    )
    parser.add_argument("exp_name", type=str, help="Name of the experiment")
    parser.add_argument("reg_name", type=str, help="Name of the region")
    parser.add_argument(
        "--type",
        type=str,
        default="merscope",
        help="Type of the data to ingest (default: merscope)",
    )
    parser.add_argument(
        "--prefix_name",
        type=str,
        default="default",
        help="Prefix for the keys in the spatialdata object (default: default)",
    )
    parser.add_argument(
        "--source",
        type=str,
        default="machine",
        help="Source of the data (default: machine)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Plot results after ingestion"
    )
    parser.add_argument(
        "--root_path", type=parse_path, default=None,
        help="Root path for the processed data (default: None, uses PROCESSED_ROOT_PATH in .env)"
    )
    parser.add_argument(
        "--zarr_store", type=parse_path, default=None,
        help="Path for the Zarr store (default: None, uses ZARR_STORAGE_PATH in .env)"
    )
    parser.add_argument(
        "--image_store", type=parse_path, default=None,
        help="Path for the image store (default: None, uses IMAGE_STORE_PATH in .env)"
    )
    return parser

def ingest_all_register_subparser(subparser):
    """Register subparser for ingest-all command."""
    parser = subparser.add_parser(
        "ingest-all",
        aliases=["ingest_all"],
        help="Ingest all regions of an experiment into spatialdata objects",
        description="Ingest all regions of an experiment into spatialdata objects",
    )
    parser.add_argument("exp_name", type=str, help="Name of the experiment")
    parser.add_argument(
        "--type",
        type=str,
        default="merscope",
        help="Type of the data to ingest (default: merscope)",
    )
    parser.add_argument(
        "--prefix_name",
        type=str,
        default="default",
        help="Prefix for the keys in the spatialdata object (default: default)",
    )
    parser.add_argument(
        "--source",
        type=str,
        default="machine",
        help="Source of the data (default: machine)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Plot results after ingestion"
    )
    parser.add_argument(
        "--root_path", type=parse_path, default=None,
        help="Root path for the processed data (default: None, uses PROCESSED_ROOT_PATH in .env)"
    )
    parser.add_argument(
        "--zarr_store", type=parse_path, default=None,
        help="Path for the Zarr store (default: None, uses ZARR_STORAGE_PATH in .env)"
    )
    parser.add_argument(
        "--image_store", type=parse_path, default=None,
        help="Path for the image store (default: None, uses IMAGE_STORE_PATH in .env)"
    )
    return parser


def load_segmentation_region_register_subparser(subparser):
    """Register subparser for load-segmentation-region command."""
    parser = subparser.add_parser(
        "load-segmentation-region",
        aliases=["load_segmentation_region"],
        help="Load segmentation data into spatialdata objects for a specific region",
        description="Load segmentation data into spatialdata objects for a specific region",
    )
    parser.add_argument("exp_name", type=str, help="Name of the experiment")
    parser.add_argument("reg_name", type=str, help="Name of the region")
    parser.add_argument(
        "seg_dir", type=parse_path, help="Directory containing segmentation data"
    )
    parser.add_argument(
        "--type",
        type=str,
        default="vpt",
        help="Type of the segmentation data to load (default: vpt)",
    )
    parser.add_argument(
        "--prefix_name",
        type=str,
        default="default",
        help="Prefix for the keys in the spatialdata object (default: default)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Plot results after loading segmentation"
    )
    # parser.add_argument('--cell_metadata_fname', type=str, default='cell_metadata.csv')
    # parser.add_argument('--cell_by_gene_fname', type=str, default='cell_by_gene.csv')
    # parser.add_argument('--detected_transcripts_fname', type=str, default='detected_transcripts.csv')
    # parser.add_argument('--cellpose_micron_space_fname', type=str, default='cellpose_micron_space.parquet')
    parser.add_argument(
        "-k",
        "--load_kwargs",
        nargs="*",
        action=ParseKwargs,
        default={},
        help="Additional keyword arguments for loading functions in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return


def load_segmentation_all_register_subparser(subparser):
    """Register subparser for load-segmentation-all command."""
    parser = subparser.add_parser(
        "load-segmentation-all",
        aliases=["load_segmentation_all"],
        help="Load segmentation data for all regions of an experiment into spatialdata objects",
        description="Load segmentation data for all regions of an experiment into spatialdata objects",
    )
    parser.add_argument("exp_name", type=str, help="Name of the experiment")
    parser.add_argument(
        "seg_dir", type=parse_path, help="Directory containing segmentation data"
    )
    parser.add_argument(
        "--type",
        type=str,
        default="vpt",
        help="Type of the segmentation data to load (default: vpt)",
    )
    parser.add_argument(
        "--prefix_name",
        type=str,
        default="default",
        help="Prefix for the keys in the spatialdata object (default: default)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Plot results after loading segmentation"
    )
    parser.add_argument(
        "-k",
        "--load_kwargs",
        nargs="*",
        action=ParseKwargs,
        default={},
        help="Additional keyword arguments for loading functions in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return


def deconwolf_register_subparser(subparser):
    """Get parser for decon_image"""
    parser = subparser.add_parser(
        "decon-image",
        aliases=["decon_image"],
        help="Deconvolve large image files in tiles",
        description="Deconvolve large image files in tiles using DeconWolf algorithm",
    )

    parser.add_argument(
        "-i",
        "--image_path",
        type=parse_path,
        required=True,
        help="Path to the image file or directory",
    )
    parser.add_argument(
        "--data_org_path",
        type=str,
        required=True,
        help="Path to data organization file",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=parse_path,
        default="tiles_output",
        help="Output directory for tiles",
    )
    parser.add_argument(
        "--channels",
        type=parse_list,
        required=True,
        help="Channel for segmentation (e.g., DAPI or PolyT,DAPI)",
    )
    parser.add_argument(
        "-ts",
        "--tile_size",
        type=int,
        default=2960,
        help="Tile size in pixels (default: 2960)",
    )
    parser.add_argument(
        "--overlap",
        type=int,
        default=100,
        help="Overlap between tiles in pixels (default: 100)",
    )
    parser.add_argument(
        "--visualize_grid", action="store_true", help="Visualize the tiling grid"
    )

    parser.add_argument(
        "--z_step", type=float, default=1.5, help="axial(z) step size in micrometers"
    )
    parser.add_argument(
        "--filter",
        type=str,
        default=None,
        help="Filter to apply to the image before segmentation",
    )
    parser.add_argument(
        "--filter_args",
        type=parse_dict,
        default={},
        help="Additional filter arguments (e.g., key1=val1,key2=val2)",
    )
    parser.add_argument("--gpu", type=bool, default=False, help="Use GPU")
    parser.add_argument(
        "--continue_stalled",
        type=bool,
        default=False,
        help="Continue processing if some tiles already processed",
    )
    parser.add_argument(
        "--plot_thr", type=bool, default=False, help="Plot thresholding histogram"
    )
    parser.add_argument(
        "--match_pre",
        type=bool,
        default=False,
        help="Match deconvolved tile intensity histogram to pre-deconvolved tiles",
    )
    return


def load_decon_images_register_subparser(subparser):
    """Register subparser for loading deconvolved images."""
    parser = subparser.add_parser(
        "load-decon-images",
        aliases=["load_decon_images"],
        help="Load deconvolved images into spatialdata objects",
        description="Load deconvolved images into spatialdata objects",
    )
    parser.add_argument("exp_name", type=str, help="Name of the experiment")
    parser.add_argument("reg_name", type=str, help="Name of the region")
    parser.add_argument(
        "--image_dir",
        type=parse_path,
        default=None,
        help="Directory containing deconvolved images (default: None)",
    )
    parser.add_argument(
        "--image_name",
        type=str,
        default="decon_images",
        help="Name of the deconvolved images in the spatialdata object(default: decon_images)",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default=".decon.tif",
        help="Suffix of the deconvolved images (default: .decon.tif)",
    )
    parser.add_argument("--plot", action="store_true", help="Plot results")
    parser.add_argument(
        "--load_kwargs",
        nargs="*",
        action=ParseKwargs,
        default={},
        help="Additional keyword arguments for loading functions in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return parser


def run_segmentation_region_register_subparser(subparser):
    # Subcommand: run_segmentation
    run_parser = subparser.add_parser(
        "run",
        aliases=["run-segmentation-region", "run_segmentation_region"],
        help="Run segmentation on a single region",
        description="Run an implemented segmentation algorithm on a given region",
    )
    run_parser.add_argument(
        "type",
        choices=["proseg", "vpt", "cellpose", "mesmer"],
        help="Type of segmentation to run",
    )
    run_parser.add_argument("exp_name", help="Name of the experiment")
    run_parser.add_argument("reg_name", help="Name of the region")
    run_parser.add_argument(
        "--input_dir",
        type=str,
        help="Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)",
    )
    run_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)",
    )
    run_parser.add_argument(
        "--config_path",
        type=str,
        default="/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json",
        help="Configuration file path (for VPT segmentation)",
    )
    run_parser.add_argument(
        "-k",
        "--kwargs",
        nargs="*",
        action=ParseKwargs,
        help="Additional keyword arguments to segmentation algorithms in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return


def run_segmentation_experiment_register_subparser(subparser):
    # Subcommand: segment_experiment
    exp_parser = subparser.add_parser(
        "experiment",
        aliases=["segment-experiment", "segment_experiment"],
        help="Run segmentation for all regions in an experiment",
        description="Run segmentation for all regions in an experiment",
    )
    exp_parser.add_argument(
        "type",
        choices=["proseg", "vpt", "cellpose", "mesmer"],
        help="Type of segmentation to run",
    )
    exp_parser.add_argument("exp_name", help="Name of the experiment")
    exp_parser.add_argument(
        "--input_dir",
        type=str,
        help="Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)",
    )
    exp_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)",
    )
    exp_parser.add_argument(
        "--config_path",
        type=str,
        default="/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json",
        help="Configuration file path (for VPT segmentation)",
    )
    exp_parser.add_argument(
        "-k",
        "--kwargs",
        type=str,
        nargs="*",
        help="Additional keyword arguments in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return


def run_proseg_alignment_register_subparser(subparser):
    # Subcommand: align_proseg
    align_parser = subparser.add_parser(
        "align",
        aliases=["align_proseg", "align-proseg"],
        help="Align Proseg transcripts to seed transcripts",
        description="Align Proseg transcripts to seed transcripts",
    )
    align_parser.add_argument("exp_name", help="Name of the experiment")
    align_parser.add_argument("reg_name", help="Name of the region")
    align_parser.add_argument(
        "--seed_prefix_name",
        type=str,
        default="default",
        help="Seed prefix name (default: default)",
    )
    align_parser.add_argument(
        "--prefix_name",
        type=str,
        default="proseg",
        help="Prefix name (default: proseg)",
    )
    align_parser.add_argument(
        "--out_prefix_name",
        type=str,
        default="proseg_aligned",
        help="Output prefix name (default: proseg_aligned)",
    )
    align_parser.add_argument(
        "--input_dir",
        type=str,
        help="Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)",
    )
    align_parser.add_argument(
        "--seg_dir",
        type=str,
        help="Segmentation directory (default: uses SEGMENTATION_OUT_PATH env var)",
    )
    align_parser.add_argument(
        "--x", type=str, default="x", help="X coordinate column name (default: x)"
    )
    align_parser.add_argument(
        "--y", type=str, default="y", help="Y coordinate column name (default: y)"
    )
    align_parser.add_argument(
        "--z",
        type=str,
        default="global_z",
        help="Z coordinate column name (default: global_z)",
    )
    align_parser.add_argument(
        "--cell_column",
        type=str,
        default="cell_id",
        help="Cell column name (default: cell_id)",
    )
    align_parser.add_argument(
        "--barcode_column",
        type=str,
        default="barcode_id",
        help="Barcode column name (default: barcode_id)",
    )
    align_parser.add_argument(
        "--gene_column",
        type=str,
        default="gene",
        help="Gene column name (default: gene)",
    )
    align_parser.add_argument(
        "--fov_column", type=str, default="fov", help="FOV column name (default: fov)"
    )
    align_parser.add_argument(
        "--cell_missing",
        type=str,
        default="-1",
        help="Cell missing value (default: -1)",
    )
    align_parser.add_argument(
        "--min_jaccard",
        type=float,
        default=0.4,
        help="Minimum Jaccard index (default: 0.4)",
    )
    align_parser.add_argument(
        "--min_prob", type=float, default=0.5, help="Minimum probability (default: 0.5)"
    )
    align_parser.add_argument(
        "--filter_blank", action="store_true", help="Filter blank genes"
    )
    align_parser.add_argument(
        "--cell_metadata_fname",
        type=str,
        default="merged_cell_metadata.csv",
        help="Cell metadata filename (default: merged_cell_metadata.csv)",
    )
    align_parser.add_argument(
        "--cell_by_gene_fname",
        type=str,
        default="merged_cell_by_gene.csv",
        help="Cell by gene filename (default: merged_cell_by_gene.csv)",
    )
    align_parser.add_argument(
        "--detected_transcripts_fname",
        type=str,
        default="merged_transcript_metadata.csv",
        help="Detected transcripts filename (default: merged_transcript_metadata.csv)",
    )
    align_parser.add_argument(
        "--cell_polygons_fname",
        type=str,
        default="merged_cell_polygons.geojson",
        help="Cell polygons filename (default: merged_cell_polygons.geojson)",
    )
    align_parser.add_argument(
        "-k",
        "--kwargs",
        type=str,
        nargs="*",
        help="Additional keyword arguments in key=value format (e.g., --kwargs param1=value1 param2=value2)",
    )
    return


def filter_to_ids_register_subparser(subparser):
    # Subcommand: filt-to-ids
    parser = subparser.add_parser(
        "filt-to-ids",
        aliases=["filt_to_ids", "filter-to-ids", "filter_to_ids"],
        help="Filter cellpose segmentation outputs to aligned proseg segmentation outputs",
        description="Filter cellpose segmentation outputs to aligned proseg segmentation outputs",
    )
    parser.add_argument(
        "--meta_path", type=parse_path, required=True, help="Path to the metadata file"
    )
    parser.add_argument(
        "--geom_path", type=parse_path, required=True, help="Path to the geometry file"
    )
    parser.add_argument(
        "--tz_path", type=parse_path, required=True, help="Path to the time zone file"
    )
    parser.add_argument(
        "--cbg_path",
        type=parse_path,
        required=True,
        help="Path to the cell by gene file",
    )
    return


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOGUE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Available commands", metavar="{run,experiment,align}"
    )

    # Adding all the subparsers from above
    current_module = sys.modules[__name__]
    for name, register_subparser_func in inspect.getmembers(
        current_module, inspect.isfunction
    ):
        if "register_subparser" in name:
            register_subparser_func(subparsers)

    # initiate args:
    args = None
    if len(sys.argv) > 1:
        if sys.argv[1] in ["-v", "--version"]:
            print(S.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        args = parser.parse_args(["--help"])
        exit()

    config_warnings()

    # setup logging here:
    if not logging.root.handlers:
        setup_logging(stdout=True, quiet=False)

    # collect args
    args_var = vars(args)
    for key, value in args_var.items():
        logger.info(f"Argument: {key} = {value}, Type: {type(value)}")

    command = args_var.pop("command").lower().replace("_", "-")
    if command in ["ingest-region"]:
        from .io.main import ingest_region as func
    elif command in ["ingest-all"]:
        from .io.main import ingest_all as func
    elif command in ["load-segmentation-region"]:
        from .io.main import load_segmentation_region as func
    elif command in ["load-segmentation-all"]:
        from .io.main import load_segmentation_all as func
    elif command in ["decon-image"]:
        from .decon_script import decon_image as func
    elif command in ["load-decon-images"]:
        from .io.main import load_deconvolution_region as func
    elif command in ["run", "run-segmentation-region"]:
        from .segmentation.main import run_segmentation as func
    elif command in ["experiment", "segment-experiment"]:
        from .segmentation.main import segment_experiment as func
    elif command in ["align", "align-proseg"]:
        from .segmentation.main import align_proseg as func
    elif command in ["filt-to-ids"]:
        from .segmentation.proseg import filt_to_ids as func
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


# import argparse
# import spida.S.decon_script as decon_script

# def main():
#     """ Main Entry Point for the S module of SPIDA"""

#     parser = argparse.ArgumentParser(prog="spida.S", description="SPIDA S module CLI")
#     subparsers = parser.add_subparsers(dest="command", required=True)

#     # Subcommand: decon_image
#     decon_parser = subparsers.add_parser(
#         'decon_image', parents=[decon_script.get_parser()], help='Deconvolve large image files in tiles'
#     )

#     args = parser.parse_args()
#     func = args.func
#     kwargs = vars(args)
#     func(**kwargs)
