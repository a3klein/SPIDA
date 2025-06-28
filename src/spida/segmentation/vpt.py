import os
from pathlib import Path
import subprocess

from dotenv import load_dotenv # type: ignore
load_dotenv()

import pandas as pd


def _add_vpt_binary():
    """
    Add the path to the vpt binary to the PATH environment variable.
    This is necessary to ensure that the vpt command can be found when running the script.
    """
    # Check if the vpt binary exists
    vpt_path = os.getenv("VPT_BIN_PATH", "/gale/netapp/home2/aklein/miniconda3/envs/vpt/bin/vpt")
    if not os.path.exists(vpt_path):
        raise FileNotFoundError(f"VPT binary not found at {vpt_path}. Please check the installation.")

    # Need to add the vpt path to the PATH variable
    os.environ['PATH'] += os.pathsep + '/gale/netapp/home2/aklein/miniconda3/envs/vpt/bin'

    # Testing VPT works
    segmentation_command = """
        vpt --help
        """
    ret = subprocess.run(segmentation_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to run. Check the installation and the PATH variable.")
    
    print("VPT binary added to PATH successfully.")


def _cli_segmentation(CONFIG_FILE:str,
                      root_dir:str,
                      output_dir:str, 
                      region:str,
                      tile_size:int=2400,
                      tile_overlap:int=200):
    """
    Run the VPT command to segment images in a specified region.
    """
    
    segmentation_command = f"""
        vpt --verbose --processes 8 \
            run-segmentation \
            --segmentation-algorithm {CONFIG_FILE} \
            --input-images {root_dir}/{region}/images/ \
            --input-micron-to-mosaic {root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv \
            --output-path {output_dir}/{region} \
            --tile-size {tile_size} \
            --tile-overlap {tile_overlap} \
        """
    ret = subprocess.run(segmentation_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to segment images.")
    return ret
    

def _cli_partition_transcripts(root_dir:str,
                               output_dir:str,
                               region:str
                               ):                  
    """
    Run the VPT command to partition transcripts in a specified region.
    """

    ### VPT command to partition transcripts
    partition_transcripts_command = f"""
        vpt partition-transcripts \
            --input-boundaries {output_dir}/{region}/cellpose_micron_space.parquet \
            --input-transcripts {root_dir}/{region}/detected_transcripts.parquet \
            --output-entity-by-gene {output_dir}/{region}/cell_by_gene.csv \
            --output-transcripts {output_dir}/{region}/detected_transcripts.csv \
        """
    ret = subprocess.run(partition_transcripts_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to partition transcripts.")
    return ret
        
    
def _cli_get_metadata(root_dir:str,
                      output_dir:str,
                      region:str):
    """
    Run the VPT command to derive metadata from the segmented images and partitioned transcripts.
    """    

    metadata_command = f"""
        vpt derive-entity-metadata \
            --input-boundaries {output_dir}/{region}/cellpose_micron_space.parquet \
            --input-entity-by-gene {output_dir}/{region}/cell_by_gene.csv \
            --output-metadata {output_dir}/{region}/cell_metadata.csv \
            --overwrite
        """
    ret = subprocess.run(metadata_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to derive entity metadata.")

    ### VPT sum the signals
    sum_signals_command = f"""
        vpt sum-signals \
            --input-images {root_dir}/{region}/images/ \
            --input-boundaries {output_dir}/{region}/cellpose_micron_space.parquet \
            --input-micron-to-mosaic {root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv \
            --output-csv {output_dir}/{region}/sum_signals.csv
        """
    ret = subprocess.run(sum_signals_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to sum signals.")
    
    #### COMBINE METADATA AND SIGNALS
    metadata = pd.read_csv(f"{output_dir}/{region}/cell_metadata.csv", index_col=0)
    signals = pd.read_csv(f"{output_dir}/{region}/sum_signals.csv", index_col=0)
    signals.index.name = "EntityID"
    pd.merge(metadata, signals, left_on="EntityID", right_on="EntityID").to_csv(f"{output_dir}/{region}/cell_metadata.csv", index=True)
    

def run_vpt(root_dir:str,
            output_dir:str,
            region:str,
            config_path:Path):
    """
    Run the VPT command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    CONFIG_FILE_PATH (Path): Path to the VPT configuration file for segmentation.
    """
    
    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    # Add VPT binary to PATH
    _add_vpt_binary()

    # Run VPT segmentation
    _cli_segmentation(CONFIG_FILE=config_path,
                      root_dir=root_dir,
                      output_dir=output_dir,
                      region=region)
    _cli_partition_transcripts(root_dir=root_dir,
                              output_dir=output_dir,
                              region=region)
    _cli_get_metadata(root_dir=root_dir,
                      output_dir=output_dir,
                      region=region)
    


def _convert_geomerty(root_dir:str, output_dir:str, region:str):
    """
    Convert the geometry of the segmented images to VPT format.
    """

    ### VPT command to convert the geometry
    convert_geometry_command = f"""
        vpt convert-geometry \
            --input-boundaries {output_dir}/{region}/polygons.parquet \
            --output-boundaries {output_dir}/{region}/cellpose_micron_space.parquet \
            --convert-to-3D \
            --input-micron-to-mosaic {root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv \
            --number-z-planes 7 \
            --spacing-z-planes 1.5 \
        """

    ret = subprocess.run(convert_geometry_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to convert geometry.")

def seg_to_vpt(root_dir:str, 
               seg_out_dir:str,
               region:str): 
    """
    Convert other segmentation outputs to VPT format.
    Parameters:
    seg_out_dir (str): The directory containing the segmentation output files (gdf parquet file).
    region (str): The name of the region to process.
    """

    assert Path(f"{seg_out_dir}/{region}/polygons.parquet").exists(), "polygons.parquer file does not exist in the specified directory."

    _add_vpt_binary()

    _convert_geomerty(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region)
    print(f"Geometry conversion completed for region {region}.")
    print(f"Converted geometry saved to {seg_out_dir}/{region}/cellpose_micron_space.parquet.")
    
    _cli_partition_transcripts(root_dir=root_dir,
                              output_dir=seg_out_dir,
                              region=region)
    _cli_get_metadata(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region)






    
