#!/usr/bin/env python3
"""
CLI interface for SPIDA segmentation pipeline using argparse.
"""
import os
import sys
import argparse
import logging
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Import the main functions
from main import run_segmentation, segment_experiment, align_proseg
from spida.utilities.script_utils import ParseKwargs, parse_kwargs

# logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s - %(message)s',
                    datefmt="%Y-%m-%dT%H:%M:%S%z")



def create_parser():
    """Create the argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        description="SPIDA Segmentation Pipeline CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='{run,experiment,align}'
    )
    
    # Subcommand: run_segmentation
    run_parser = subparsers.add_parser(
        'run',
        help='Run segmentation on a single region',
        description='Run an implemented segmentation algorithm on a given region'
    )
    run_parser.add_argument('type', choices=['proseg', 'vpt', 'cellpose', 'mesmer'], help='Type of segmentation to run')
    run_parser.add_argument('exp_name', help='Name of the experiment')
    run_parser.add_argument('reg_name', help='Name of the region')
    run_parser.add_argument('--input_dir',type=str, help='Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)')
    run_parser.add_argument('--output_dir',type=str,help='Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)')
    run_parser.add_argument(
        '--config_path',type=str,default="/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json",
        help='Configuration file path (for VPT segmentation)'
    )
    run_parser.add_argument('-k', '--kwargs',nargs='*',action=ParseKwargs,
                            help='Additional keyword arguments to segmentation algorithms in key=value format (e.g., --kwargs param1=value1 param2=value2)')
    
    
    # Subcommand: segment_experiment
    exp_parser = subparsers.add_parser(
        'experiment',
        help='Run segmentation for all regions in an experiment',
        description='Run segmentation for all regions in an experiment'
    )
    exp_parser.add_argument('type',choices=['proseg', 'vpt', 'cellpose', 'mesmer'],help='Type of segmentation to run')
    exp_parser.add_argument('exp_name',help='Name of the experiment')
    exp_parser.add_argument('--input_dir',type=str,help='Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)')
    exp_parser.add_argument('--output_dir',type=str,help='Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)')
    exp_parser.add_argument('--config_path',type=str,default="/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json",
                            help='Configuration file path (for VPT segmentation)')
    exp_parser.add_argument('-k','--kwargs',type=str,nargs='*',
        help='Additional keyword arguments in key=value format (e.g., --kwargs param1=value1 param2=value2)')
    
    # Subcommand: align_proseg
    align_parser = subparsers.add_parser(
        'align',
        help='Align Proseg transcripts to seed transcripts',
        description='Align Proseg transcripts to seed transcripts'
    )
    align_parser.add_argument('exp_name',help='Name of the experiment')
    align_parser.add_argument('reg_name',help='Name of the region')
    align_parser.add_argument('--seed-prefix-name',type=str,default='default',help='Seed prefix name (default: default)')
    align_parser.add_argument('--prefix-name',type=str,default='proseg',help='Prefix name (default: proseg)')
    align_parser.add_argument('--input-dir',type=str,help='Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)')
    align_parser.add_argument('--seg-dir',type=str,help='Segmentation directory (default: uses SEGMENTATION_OUT_PATH env var)')
    align_parser.add_argument('--x',type=str,default='x',help='X coordinate column name (default: x)')
    align_parser.add_argument('--y',type=str,default='y',help='Y coordinate column name (default: y)')
    align_parser.add_argument('--z',type=str,default='global_z',help='Z coordinate column name (default: global_z)')
    align_parser.add_argument('--cell-column',type=str,default='cell_id',help='Cell column name (default: cell_id)')
    align_parser.add_argument('--barcode-column',type=str,default='barcode_id',help='Barcode column name (default: barcode_id)')
    align_parser.add_argument('--gene-column',type=str,default='gene',help='Gene column name (default: gene)')
    align_parser.add_argument('--fov-column',type=str,default='fov',help='FOV column name (default: fov)')
    align_parser.add_argument('--cell-missing',type=int,default=-1,help='Cell missing value (default: -1)')
    align_parser.add_argument('--min-jaccard',type=float,default=0.4,help='Minimum Jaccard index (default: 0.4)')
    align_parser.add_argument('--min-prob',type=float,default=0.5,help='Minimum probability (default: 0.5)')
    align_parser.add_argument('--filter-blank',action='store_true',help='Filter blank genes')
    align_parser.add_argument('--cell-metadata-fname',type=str,default='merged_cell_metadata.csv',
                              help='Cell metadata filename (default: merged_cell_metadata.csv)')
    align_parser.add_argument('--cell-by-gene-fname',type=str,default='merged_cell_by_gene.csv',
                              help='Cell by gene filename (default: merged_cell_by_gene.csv)')
    align_parser.add_argument('--detected-transcripts-fname',type=str,default='merged_transcript_metadata.csv',
                              help='Detected transcripts filename (default: merged_transcript_metadata.csv)')
    align_parser.add_argument('--cell-polygons-fname',type=str,default='merged_cell_polygons.geojson',
                              help='Cell polygons filename (default: merged_cell_polygons.geojson)')
    align_parser.add_argument('-k', '--kwargs',type=str,nargs='*',
                              help='Additional keyword arguments in key=value format (e.g., --kwargs param1=value1 param2=value2)')
    
    return parser


def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.command == 'run':
            # Convert input/output dirs to Path if provided
            input_dir = Path(args.input_dir) if args.input_dir else None
            output_dir = Path(args.output_dir) if args.output_dir else None
            
            kwargs = args.kwargs
            # Log the KWARGS: 
            for key, value in kwargs.items():
                logging.info(f"KWARG: {key} = {value}, Type: {type(value)}")
            
            if hasattr(args, 'config_path') and args.config_path:
                kwargs['config_path'] = args.config_path
            
            # # Parse additional kwargs from command line
            # if hasattr(args, 'kwargs') and args.kwargs:
            #     additional_kwargs = parse_kwargs(args.kwargs)
            #     kwargs.update(additional_kwargs)
            
            # # Log the KWARGS: 
            # for key, value in kwargs.items():
            #     logging.info(f"KWARG: {key} = {value}")
            
            run_segmentation(
                type=args.type,
                exp_name=args.exp_name,
                reg_name=args.reg_name,
                input_dir=input_dir,
                output_dir=output_dir,
                **kwargs
            )
            
        elif args.command == 'experiment':
            # Convert input/output dirs to Path if provided
            input_dir = Path(args.input_dir) if args.input_dir else None
            output_dir = Path(args.output_dir) if args.output_dir else None
            
            # Prepare kwargs for additional parameters
            kwargs = {}
            if hasattr(args, 'config_path') and args.config_path:
                kwargs['config_path'] = args.config_path
            
            # Parse additional kwargs from command line
            if hasattr(args, 'kwargs') and args.kwargs:
                additional_kwargs = parse_kwargs(args.kwargs)
                kwargs.update(additional_kwargs)
            
            segment_experiment(
                type=args.type,
                exp_name=args.exp_name,
                input_dir=input_dir,
                output_dir=output_dir,
                **kwargs
            )
            
        elif args.command == 'align':
            # Convert dirs to Path if provided
            input_dir = Path(args.input_dir) if args.input_dir else None
            seg_dir = Path(args.seg_dir) if args.seg_dir else None
            
            # Parse additional kwargs from command line
            additional_kwargs = {}
            if hasattr(args, 'kwargs') and args.kwargs:
                additional_kwargs = parse_kwargs(args.kwargs)
            
            align_proseg(
                exp_name=args.exp_name,
                reg_name=args.reg_name,
                seed_prefix_name=args.seed_prefix_name,
                prefix_name=args.prefix_name,
                input_dir=input_dir,
                seg_dir=seg_dir,
                x=args.x,
                y=args.y,
                z=args.z,
                cell_column=args.cell_column,
                barcode_column=args.barcode_column,
                gene_column=args.gene_column,
                fov_column=args.fov_column,
                cell_missing=args.cell_missing,
                min_jaccard=args.min_jaccard,
                min_prob=args.min_prob,
                filter_blank=args.filter_blank,
                cell_metadata_fname=args.cell_metadata_fname,
                cell_by_gene_fname=args.cell_by_gene_fname,
                detected_transcripts_fname=args.detected_transcripts_fname,
                cell_polygons_fname=args.cell_polygons_fname,
                **additional_kwargs
            )
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
