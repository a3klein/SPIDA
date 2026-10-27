import os
from pathlib import Path
import subprocess
import logging

import pandas as pd

from spida.io.formats.merscope import IMAGES_DIR, TRANSCRIPTS_CSV, TRANSCRIPTS_PARQUET, TRANSFORM_FILE

logger = logging.getLogger(__package__)


def _add_vpt_binary(vpt_bin_path: str = None):
    """
    Add the path to the vpt binary to the PATH environment variable.
    This is necessary to ensure that the vpt command can be found when running the script.
    """
    # Check if the vpt binary exists
    _native_hint = (
        "The `vpt` (vizgen-postprocessing) binary is required only for the legacy "
        "`--backend vpt` fallback. Install VPT and set VPT_BIN_PATH in the config (or pass "
        "--vpt_bin_path), or use the default `--backend native` pure-Python path, "
        "which needs no external binary."
    )
    if vpt_bin_path is None:
        raise FileNotFoundError(f"No VPT binary path given. {_native_hint}")
    vpt_path = str(vpt_bin_path)
    logger.info(f"Using VPT binary path: {vpt_path}")
    if not os.path.exists(vpt_path):
        raise FileNotFoundError(
            f"VPT binary not found at {vpt_path}. {_native_hint}"
        )

    # Need to add the vpt path to the PATH variable
    os.environ["PATH"] += os.pathsep + vpt_path

    # Testing VPT works
    try:
        subprocess.run(["vpt", "--help"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError) as e:
        raise RuntimeError(
            f"Found the VPT path at {vpt_path} but `vpt --help` failed to run "
            f"({e}). {_native_hint}"
        ) from e

    logger.info("VPT binary added to PATH successfully.")


def _cli_segmentation(
    CONFIG_FILE: str,
    raw_dir: str,
    seg_dir: str,
    tile_size: int = 2400,
    tile_overlap: int = 200,
):
    """
    Run the VPT command to segment images in a specified region.
    """

    segmentation_command = f"""
        vpt --verbose --processes 32 \n
            run-segmentation \n
            --segmentation-algorithm {CONFIG_FILE} \n
            --input-images {raw_dir}/{IMAGES_DIR}/ \n
            --input-micron-to-mosaic {raw_dir}/{IMAGES_DIR}/{TRANSFORM_FILE} \n
            --output-path {seg_dir} \n
            --tile-size {tile_size} \n
            --tile-overlap {tile_overlap} \n
            --overwrite
        """
    logger.info(f"Running VPT segmentation command: {segmentation_command}")

    ret = subprocess.run(segmentation_command.split(), capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError("VPT failed to segment images.")
    return ret


def _cli_partition_transcripts(
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "cellpose_micron_space.parquet",
    input_transcripts: str = TRANSCRIPTS_PARQUET,
    output_entity_by_gene: str = "cell_by_gene.csv",
    output_transcripts: str = TRANSCRIPTS_CSV,
):
    """
    Run the VPT command to partition transcripts in a specified region.
    """

    ### VPT command to partition transcripts
    partition_transcripts_command = f"""
        vpt --verbose --processes 32 \n
            partition-transcripts \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --input-transcripts {raw_dir}/{input_transcripts} \n
            --output-entity-by-gene {seg_dir}/{output_entity_by_gene} \n
            --output-transcripts {seg_dir}/{output_transcripts} \n
            --overwrite
        """
    logger.info(
        f"Running VPT partition transcripts command: {partition_transcripts_command}"
    )

    ret = subprocess.run(
        partition_transcripts_command.split(), capture_output=True, check=True
    )
    if ret.returncode != 0:
        raise ValueError("VPT failed to partition transcripts.")
    return ret


def _cli_get_metadata(
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "cellpose_micron_space.parquet",
    input_entity_by_gene: str = "cell_by_gene.csv",
    output_metadata: str = "cell_metadata.csv",
    output_signals: str = "sum_signals.csv",
):
    """
    Run the VPT command to derive metadata from the segmented images and partitioned transcripts.
    """

    metadata_command = f"""
        vpt --verbose --processes 32 \n
            derive-entity-metadata \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --input-entity-by-gene {seg_dir}/{input_entity_by_gene} \n
            --output-metadata {seg_dir}/{output_metadata} \n
            --overwrite
        """

    logger.info(f"Running VPT derive metadata command: {metadata_command}")

    ret = subprocess.run(metadata_command.split(), capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError("VPT failed to derive entity metadata.")

    ### VPT sum the signals
    sum_signals_command = f"""
        vpt --verbose --processes 32 \n
            sum-signals \n
            --input-images {raw_dir}/{IMAGES_DIR}/ \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --input-micron-to-mosaic {raw_dir}/{IMAGES_DIR}/{TRANSFORM_FILE} \n
            --output-csv {seg_dir}/{output_signals} \n
            --overwrite

        """

    logger.info(f"Running VPT sum signals command: {sum_signals_command}")

    ret = subprocess.run(sum_signals_command.split(), capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError("VPT failed to sum signals.")

    #### COMBINE METADATA AND SIGNALS
    metadata = pd.read_csv(f"{seg_dir}/{output_metadata}", index_col=0)
    signals = pd.read_csv(f"{seg_dir}/{output_signals}", index_col=0)
    signals.index.name = "EntityID"
    pd.merge(metadata, signals, left_on="EntityID", right_on="EntityID").to_csv(
        f"{seg_dir}/{output_metadata}", index=True
    )


def _cli_sum_signals(
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "cellpose_micron_space.parquet",
    output_signals: str = "sum_signals.csv",
): 
    """
    Run the VPT command to sum signals in a specified region.
    """

    sum_signals_command = f"""
        vpt --verbose --processes 32 \n
            sum-signals \n
            --input-images {raw_dir}/{IMAGES_DIR}/ \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --input-micron-to-mosaic {raw_dir}/{IMAGES_DIR}/{TRANSFORM_FILE} \n
            --output-csv {seg_dir}/{output_signals} \n
            --overwrite

        """

    logger.info(f"Running VPT sum signals command: {sum_signals_command}")

    ret = subprocess.run(sum_signals_command.split(), capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError("VPT failed to sum signals.")
    signals = pd.read_csv(f"{seg_dir}/{output_signals}", index_col=0)
    
    return signals


def run_vpt(
    raw_dir: str, seg_dir: str, config_path: Path, vpt_bin_path: str = None, **vpt_filepaths
):
    """
    Run the VPT command to process detected transcripts in a specified region.

    Parameters:
    raw_dir (str): The raw MERSCOPE region dir (images/, detected_transcripts).
    seg_dir (str): This segmentation's output dir for the region.
    CONFIG_FILE_PATH (Path): Path to the VPT configuration file for segmentation.
    vpt_bin_path (str): Path to the VPT binary.
    vpt_filepaths (dict): Additional file paths required for VPT processing, such as:
                    - input_boundaries
                    - input_transcripts
                    - output_entity_by_gene
                    - output_transcripts
                    - output_metadata
                    - output_signals
    """

    # Ensure the output directory exists
    Path(f"{seg_dir}").mkdir(parents=True, exist_ok=True)

    # Add VPT binary to PATH
    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    # Run VPT segmentation
    _cli_segmentation(
        CONFIG_FILE=config_path,
        raw_dir=raw_dir,
        seg_dir=seg_dir,
        **vpt_filepaths,
    )

    _cli_partition_transcripts(
        raw_dir=raw_dir, seg_dir=seg_dir, **vpt_filepaths
    )

    _cli_get_metadata(
        raw_dir=raw_dir, seg_dir=seg_dir, **vpt_filepaths
    )


def _convert_geometry(
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "polygons.parquet",
    output_boundaries: str = "cellpose_micron_space.parquet",
    convert_micron: bool = True,
    is_3d: bool = False,
):
    """
    Convert the geometry of the segmented images to VPT format.
    """

    ### VPT command to convert the geometry
    convert_geometry_command = f"""
        vpt --verbose --processes 16 \n
            convert-geometry \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --output-boundaries {seg_dir}/{output_boundaries} \n
            --convert-to-3D \n
            --number-z-planes 7 \n
            --spacing-z-planes 1.5 \n
            --overwrite \n
        """
    if is_3d: 
        convert_geometry_command = f"""
        vpt --verbose --processes 16 \n
            convert-geometry \n
            --input-boundaries {seg_dir}/{input_boundaries} \n
            --output-boundaries {seg_dir}/{output_boundaries} \n
            --overwrite \n
        """

    if convert_micron:
        convert_geometry_command += f" --input-micron-to-mosaic {raw_dir}/{IMAGES_DIR}/{TRANSFORM_FILE}"

    logger.info(f"Running VPT convert geometry command: {convert_geometry_command}")

    ret = subprocess.run(
        convert_geometry_command.split(), capture_output=True, check=True
    )
    if ret.returncode != 0:
        raise ValueError("VPT failed to convert geometry.")


def seg_to_vpt(
    raw_dir: str,
    seg_dir: str,
    vpt_bin_path: str = None,
    is_3d : bool = False,
    spacing_z: float = 1.5,
    **vpt_filepaths):
    """
    Convert other segmentation outputs to VPT format.
    Parameters:
    raw_dir (str): The raw MERSCOPE region dir (images/, detected_transcripts).
    seg_dir (str): This segmentation's output dir for the region (holds the gdf parquet).
    is_3d (bool): Whether the data is 3D or not. Default is False (2D).
    spacing_z (float): Micron spacing between adjacent z-planes; used for the 3D
        ZLevel = (z + 1) * spacing_z (1-based) derivation. Only used when is_3d. Default 1.5.
    vpt_filepaths (dict): Additional file paths required for VPT processing, such as:
        - input_boundaries: Path to the input boundaries file (e.g., polygons.parquet).
        - output_boundaries: Path to save the converted boundaries file (e.g., cellpose_micron_space.parquet).
        - input_transcripts: Path to the input transcripts file (e.g., detected_transcripts.parquet).
        - output_transcripts: Path to save the converted transcripts file (e.g., detected_transcripts.csv).
        - output_entity_by_gene: Path to save the entity by gene file (e.g., cell_by_gene.csv).
        - output_metadata: Path to save the metadata file (e.g., cell_metadata.csv).
        - output_signals: Path to save the signals file (e.g., sum_signals.csv).

    """

    input_boundaries = vpt_filepaths.get('input_boundaries', 'polygons.parquet')
    assert Path(
        f"{seg_dir}/{input_boundaries}"
    ).exists(), f"{input_boundaries} file does not exist in the specified directory: \n{seg_dir}/"

    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    logger.info(f"Converting Geometry in {seg_dir}...")

    output_boundaries=vpt_filepaths.get(
        "output_boundaries", "cellpose_micron_space.parquet"
    )
    if not is_3d: 
        _convert_geometry(
            raw_dir=raw_dir,
            seg_dir=seg_dir,
            input_boundaries=input_boundaries,
            output_boundaries=output_boundaries,
        )

    else: 
        from .vpt_utils import _3d_convert_geometry
        _3d_convert_geometry(
            raw_dir=raw_dir,
            seg_dir=seg_dir,
            input_boundaries=input_boundaries,
            output_boundaries=output_boundaries,
            spacing_z=spacing_z,
        )
        

    logger.info(f"Geometry conversion completed in {seg_dir}.")
    logger.info(
        f"Converted geometry saved to {seg_dir}/{output_boundaries}."
    )

    dt_fname = vpt_filepaths.get("input_transcripts", TRANSCRIPTS_PARQUET)
    dt_parquet = Path(f"{raw_dir}/{dt_fname}")
    logger.info(f"Checking for detected transcripts file at {dt_parquet}")
    if not dt_parquet.exists():
        dt_fname = TRANSCRIPTS_CSV
        assert Path(
            f"{raw_dir}/{dt_fname}"
        ).exists(), f"{TRANSCRIPTS_CSV} file does not exist in the specified root directory."

    # Run VPT partition transcripts
    _cli_partition_transcripts(
        raw_dir=raw_dir,
        seg_dir=seg_dir,
        input_boundaries=output_boundaries,
        input_transcripts=dt_fname,
        output_entity_by_gene=vpt_filepaths.get(
            "output_entity_by_gene", "cell_by_gene.csv"
        ),
        output_transcripts=vpt_filepaths.get(
            "output_transcripts", TRANSCRIPTS_CSV
        ),
    )

    logger.info(f"Partitioned transcripts in {seg_dir}.")

    # Run VPT get metadata
    _cli_get_metadata(
        raw_dir=raw_dir,
        seg_dir=seg_dir,
        input_boundaries=output_boundaries,
        input_entity_by_gene=vpt_filepaths.get(
            "output_entity_by_gene", "cell_by_gene.csv"
        ),
        output_metadata=vpt_filepaths.get("output_metadata", "cell_metadata.csv"),
        output_signals=vpt_filepaths.get("output_signals", "sum_signals.csv"),
    )

    logger.info(f"Metadata and signals derived in {seg_dir}.")


def generate_metadata(
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "cellpose_micron_space.parquet",
    output_boundaries: str = "cellpose_micron_space.parquet",
    input_entity_by_gene: str = "cell_by_gene.csv",
    output_entity_by_gene: str = "cell_by_gene.csv",
    input_transcripts: str = TRANSCRIPTS_PARQUET,
    output_transcripts: str = TRANSCRIPTS_CSV,
    output_metadata: str = "cell_metadata.csv",
    output_signals: str = "sum_signals.csv",
    vpt_bin_path: str = None,
    is_3d : bool = False,
):
    """
    Generate metadata from the segmented images and partitioned transcripts.
    Parameters:
    raw_dir (str): The raw MERSCOPE region dir (images/, detected_transcripts).
    seg_dir (str): This segmentation's output dir for the region.
    input_boundaries (str): Path to the input boundaries file (default is "cellpose_micron_space.parquet").
    output_boundaries (str): Path to save the converted boundaries file (default is "cellpose_micron_space.parquet").
    input_entity_by_gene (str): Path to the input entity by gene file (default is "cell_by_gene.csv").
    output_metadata (str): Path to save the metadata file (default is "cell_metadata.csv").
    output_signals (str): Path to save the signals file (default is "sum_signals.csv").
    """

    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    _convert_geometry(
        raw_dir=raw_dir,
        seg_dir=seg_dir,
        input_boundaries=input_boundaries,
        output_boundaries=output_boundaries,
        convert_micron=False,
        is_3d=is_3d
    )
    logger.info(f"Geometry conversion completed in {seg_dir}.")
    logger.info(
        f"Converted geometry saved to {seg_dir}/{output_boundaries}."
    )

    _cli_partition_transcripts(
        raw_dir=seg_dir,
        seg_dir=seg_dir,
        input_boundaries=output_boundaries,
        input_transcripts=input_transcripts,
        output_entity_by_gene=output_entity_by_gene,
        output_transcripts=output_transcripts,
    )
    logger.info(f"Partitioned transcripts in {seg_dir}.")
    logger.info(
        f"Partitioned transcripts saved to {seg_dir}/{output_transcripts}."
    )

    _cli_get_metadata(
        raw_dir=raw_dir,
        seg_dir=seg_dir,
        input_boundaries=output_boundaries,
        input_entity_by_gene=input_entity_by_gene,
        output_metadata=output_metadata,
        output_signals=output_signals,
    )

    logger.info(f"Metadata and signals derived in {seg_dir}.")
