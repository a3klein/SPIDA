import os
from pathlib import Path
import subprocess
import logging

from dotenv import load_dotenv # type: ignore
import pandas as pd

load_dotenv()
logger = logging.getLogger(__package__)


def _add_vpt_binary():
    """
    Add the path to the vpt binary to the PATH environment variable.
    This is necessary to ensure that the vpt command can be found when running the script.
    """
    # Check if the vpt binary exists
    vpt_path = os.getenv("VPT_BIN_PATH", "/gale/netapp/home2/aklein/miniconda3/envs/vpt/bin")
    if not os.path.exists(vpt_path):
        raise FileNotFoundError(f"VPT binary not found at {vpt_path}. Please check the installation.")

    # Need to add the vpt path to the PATH variable
    os.environ['PATH'] += os.pathsep + vpt_path

    # Testing VPT works
    segmentation_command = """
        vpt --help
        """
    ret = subprocess.run(segmentation_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to run. Check the installation and the PATH variable.")

    logger.info("VPT binary added to PATH successfully.")


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
            --overwrite
        """
    logger.info(f"Running VPT segmentation command: {segmentation_command}")

    ret = subprocess.run(segmentation_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to segment images.")
    return ret
    

def _cli_partition_transcripts(
        root_dir:str,
        output_dir:str,
        region:str,
        input_boundaries:str = "cellpose_micron_space.parquet",
        input_transcripts:str = "detected_transcripts.parquet",
        output_entity_by_gene:str = "cell_by_gene.csv",
        output_transcripts:str = "detected_transcripts.csv"
        ):                  
    """
    Run the VPT command to partition transcripts in a specified region.
    """

    ### VPT command to partition transcripts
    partition_transcripts_command = f"""
        vpt --verbose --processes 8 \
            partition-transcripts \
            --input-boundaries {output_dir}/{region}/{input_boundaries} \
            --input-transcripts {root_dir}/{region}/{input_transcripts} \
            --output-entity-by-gene {output_dir}/{region}/{output_entity_by_gene} \
            --output-transcripts {output_dir}/{region}/{output_transcripts} \
            --overwrite
        """
    logger.info(f"Running VPT partition transcripts command: {partition_transcripts_command}")

    ret = subprocess.run(partition_transcripts_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to partition transcripts.")
    return ret
        
    
def _cli_get_metadata(
        root_dir:str,
        output_dir:str,
        region:str, 
        input_boundaries:str = "cellpose_micron_space.parquet",
        input_entity_by_gene:str = "cell_by_gene.csv",
        output_metadata:str = "cell_metadata.csv",
        output_signals:str = "sum_signals.csv",

        ):
    """
    Run the VPT command to derive metadata from the segmented images and partitioned transcripts.
    """    

    metadata_command = f"""
        vpt --verbose --processes 8 \
            derive-entity-metadata \
            --input-boundaries {output_dir}/{region}/{input_boundaries} \
            --input-entity-by-gene {output_dir}/{region}/{input_entity_by_gene} \
            --output-metadata {output_dir}/{region}/{output_metadata} \
            --overwrite
        """
    
    logger.info(f"Running VPT derive metadata command: {metadata_command}")

    ret = subprocess.run(metadata_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to derive entity metadata.")

    ### VPT sum the signals
    sum_signals_command = f"""
        vpt --verbose --processes 8 \
            sum-signals \
            --input-images {root_dir}/{region}/images/ \
            --input-boundaries {output_dir}/{region}/{input_boundaries} \
            --input-micron-to-mosaic {root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv \
            --output-csv {output_dir}/{region}/{output_signals} \
            --overwrite

        """
    
    logger.info(f"Running VPT sum signals command: {sum_signals_command}")

    ret = subprocess.run(sum_signals_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to sum signals.")
    
    #### COMBINE METADATA AND SIGNALS
    metadata = pd.read_csv(f"{output_dir}/{region}/{output_metadata}", index_col=0)
    signals = pd.read_csv(f"{output_dir}/{region}/{output_signals}", index_col=0)
    signals.index.name = "EntityID"
    pd.merge(metadata, signals, left_on="EntityID", right_on="EntityID").to_csv(f"{output_dir}/{region}/{output_metadata}", index=True)

def run_vpt(root_dir:str,
            output_dir:str,
            region:str,
            config_path:Path,
            **vpt_filepaths):
    """
    Run the VPT command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    CONFIG_FILE_PATH (Path): Path to the VPT configuration file for segmentation.
    vpt_filepaths (dict): Additional file paths required for VPT processing, such as:
                    - input_boundaries
                    - input_transcripts
                    - output_entity_by_gene
                    - output_transcripts
                    - output_metadata
                    - output_signals
    """
    
    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    # Add VPT binary to PATH
    _add_vpt_binary()

    # Run VPT segmentation
    _cli_segmentation(CONFIG_FILE=config_path,
                      root_dir=root_dir,
                      output_dir=output_dir,
                      region=region,
                      **vpt_filepaths)
    
    _cli_partition_transcripts(root_dir=root_dir,
                              output_dir=output_dir,
                              region=region,
                              **vpt_filepaths)
    
    _cli_get_metadata(root_dir=root_dir,
                      output_dir=output_dir,
                      region=region,
                      **vpt_filepaths)
    


def _convert_geometry(
        root_dir:str,
        output_dir:str,
        region:str,
        input_boundaries:str = "polygons.parquet",
        output_boundaries:str = "cellpose_micron_space.parquet",
        convert_micron:bool = True,
        ):
    """
    Convert the geometry of the segmented images to VPT format.
    """

    ### VPT command to convert the geometry
    convert_geometry_command = f"""
        vpt --verbose --processes 8 \
            convert-geometry \
            --input-boundaries {output_dir}/{region}/{input_boundaries} \
            --output-boundaries {output_dir}/{region}/{output_boundaries} \
            --convert-to-3D \
            --number-z-planes 7 \
            --spacing-z-planes 1.5 \
            --overwrite \
        """
    
    if convert_micron:
        convert_geometry_command += f" --input-micron-to-mosaic {root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv"


    ret = subprocess.run(convert_geometry_command.split(), capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("VPT failed to convert geometry.")

def seg_to_vpt(root_dir:str, 
               seg_out_dir:str,
               region:str,
               **vpt_filepaths): 
    """
    Convert other segmentation outputs to VPT format.
    Parameters:
    seg_out_dir (str): The directory containing the segmentation output files (gdf parquet file).
    region (str): The name of the region to process.
    vpt_filepaths (dict): Additional file paths required for VPT processing, such as:
        - input_boundaries: Path to the input boundaries file (e.g., polygons.parquet).
        - output_boundaries: Path to save the converted boundaries file (e.g., cellpose_micron_space.parquet).
        - input_transcripts: Path to the input transcripts file (e.g., detected_transcripts.parquet).
        - output_transcripts: Path to save the converted transcripts file (e.g., detected_transcripts.csv).
        - output_entity_by_gene: Path to save the entity by gene file (e.g., cell_by_gene.csv).
        - output_metadata: Path to save the metadata file (e.g., cell_metadata.csv).
        - output_signals: Path to save the signals file (e.g., sum_signals.csv).

    """

    assert Path(f"{seg_out_dir}/{region}/{vpt_filepaths.get('input_boundaries', 'polygons.parquet')}").exists(), "polygons.parquet file does not exist in the specified directory."

    _add_vpt_binary()

    logger.info(f"Converting Geometry for region {region}...")

    _convert_geometry(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region,
                      input_boundaries=vpt_filepaths.get("input_boundaries", "polygons.parquet"),
                      output_boundaries=vpt_filepaths.get("output_boundaries", "cellpose_micron_space.parquet"),
                      )

    logger.info(f"Geometry conversion completed for region {region}.")
    logger.info(f"Converted geometry saved to {seg_out_dir}/{region}/cellpose_micron_space.parquet.")

    # Run VPT partition transcripts
    _cli_partition_transcripts(root_dir=root_dir,
                              output_dir=seg_out_dir,
                              region=region,
                              input_boundaries=vpt_filepaths.get("output_boundaries", "cellpose_micron_space.parquet"),
                              input_transcripts=vpt_filepaths.get("input_transcripts", "detected_transcripts.parquet"),
                              output_entity_by_gene=vpt_filepaths.get("output_entity_by_gene", "cell_by_gene.csv"),
                              output_transcripts=vpt_filepaths.get("output_transcripts", "detected_transcripts.csv"),
                                )

    logger.info(f"Partitioned transcripts for region {region}.")

    # Run VPT get metadata
    _cli_get_metadata(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region,
                      input_boundaries=vpt_filepaths.get("output_boundaries", "cellpose_micron_space.parquet"),
                      input_entity_by_gene=vpt_filepaths.get("output_entity_by_gene", "cell_by_gene.csv"),
                      output_metadata=vpt_filepaths.get("output_metadata", "cell_metadata.csv"),
                      output_signals=vpt_filepaths.get("output_signals", "sum_signals.csv"),
                      )

    logger.info(f"Metadata and signals derived for region {region}.")

def generate_metadata(root_dir:str,
                      seg_out_dir:str,
                      region:str,
                      input_boundaries:str = "cellpose_micron_space.parquet",
                      output_boundaries:str = "cellpose_micron_space.parquet",
                      input_entity_by_gene:str = "cell_by_gene.csv",
                      output_entity_by_gene:str = "cell_by_gene.csv",
                      input_transcripts:str = "detected_transcripts.parquet",
                      output_transcripts:str = "detected_transcripts.csv",
                      output_metadata:str = "cell_metadata.csv",
                      output_signals:str = "sum_signals.csv",
                      ):
    """
    Generate metadata from the segmented images and partitioned transcripts.
    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    seg_out_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    input_boundaries (str): Path to the input boundaries file (default is "cellpose_micron_space.parquet").
    output_boundaries (str): Path to save the converted boundaries file (default is "cellpose_micron_space.parquet").
    input_entity_by_gene (str): Path to the input entity by gene file (default is "cell_by_gene.csv").
    output_metadata (str): Path to save the metadata file (default is "cell_metadata.csv").
    output_signals (str): Path to save the signals file (default is "sum_signals.csv").
    """

    _add_vpt_binary()

    _convert_geometry(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region,
                      input_boundaries=input_boundaries,
                      output_boundaries=output_boundaries,
                      convert_micron=False
                      )
    logger.info(f"Geometry conversion completed for region {region}.")
    logger.info(f"Converted geometry saved to {seg_out_dir}/{region}/{output_boundaries}.")
    
    _cli_partition_transcripts(root_dir=seg_out_dir,
                              output_dir=seg_out_dir,
                              region=region,
                              input_boundaries=output_boundaries,
                              input_transcripts=input_transcripts,
                              output_entity_by_gene=output_entity_by_gene,
                              output_transcripts=output_transcripts,
                                )
    logger.info(f"Partitioned transcripts for region {region}.")
    logger.info(f"Partitioned transcripts saved to {seg_out_dir}/{region}/{output_transcripts}.")

    _cli_get_metadata(root_dir=root_dir,
                      output_dir=seg_out_dir,
                      region=region,
                      input_boundaries=output_boundaries,
                      input_entity_by_gene=input_entity_by_gene,
                      output_metadata=output_metadata,
                      output_signals=output_signals)

    logger.info(f"Metadata and signals derived for region {region}.")



    






    
