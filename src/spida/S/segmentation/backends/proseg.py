import os

import logging
import warnings
import pathlib
import subprocess
import re
import gzip

import pandas as pd
import geopandas as gpd

logger = logging.getLogger(__package__)


def _determine_proseg_version():
    """
    Determine the version of the proseg binary installed.
    Returns:
        str: The version of the proseg binary.
    """
    try:
        ret = subprocess.run(["proseg", "-V"], capture_output=True, check=True)
        if ret.returncode != 0:
            raise ValueError("Proseg failed to run. Check the installation and the PATH variable.")
        proseg_version = ret.stdout.decode('utf-8').split(" ")[-1].strip()
        return proseg_version
    except FileNotFoundError:
        raise FileNotFoundError("Proseg binary not found. Please ensure it is installed and in your PATH.")

# click? 
def get_proseg_version():
    """
    Get the version of the proseg binary installed.
    Returns:
        str: The version of the proseg binary.
    """
    try:
        proseg_version = _determine_proseg_version()
        logger.info(f"Proseg version: {proseg_version}")
        return proseg_version
    except Exception as e:
        logger.error(f"Error determining Proseg version: {e}")
        raise

def _add_proseg_binary(rust_bin_path : str | None = None):
    """
    Add the path to the proseg binary to the PATH environment variable.
    This is necessary to ensure that the proseg command can be found when running the script.
    Returns:
        bool: True if the proseg version is 3.x, False otherwise.
    """
    # Need to add the rust path to the PATH variable
    if rust_bin_path is None:
        raise ValueError("no proseg binary dir; pass --rust_bin_path or set RUST_BIN_PATH in the config")
    rust_path = pathlib.Path(rust_bin_path)
    if not rust_path.exists():
        raise ValueError(f"RUST_BIN_PATH {rust_path} does not exist. Please set the correct path.")
    os.environ["PATH"] += os.pathsep + str(rust_path)

    # Test that proseg can be run:
    ret = subprocess.run(["proseg", "--help"], capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError(
            "Proseg failed to run. Check the installation and the PATH variable."
        )
    proseg_v3 = _determine_proseg_version()[0] == '3'
    logger.info(f"Proseg version 3.x detected: {proseg_v3}")

    proseg_out = ret.stdout.decode()
    proseg_options = proseg_out.split("Options:")[-1]

    # Extract the options from the proseg output
    global proseg_args
    escaped_start = re.escape("--")
    escaped_end = re.escape("\n")
    pattern = f"{escaped_start}(.*?){escaped_end}"
    proseg_args = re.findall(pattern, proseg_options.strip())
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    proseg_args = [ansi_escape.sub("", arg) for arg in proseg_args]
    proseg_args = [arg.split()[0] for arg in proseg_args]

    logger.info("proseg binary added to PATH successfully.")

    return proseg_v3




def _execute_cli_proseg_v3(root_dir: str, output_dir: str, region: str, **proseg_params):
    """
    Run the proseg command for the 3.0.1 version of proseg
    """
    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \n
        --output-path {output_dir} \n
        --output-spatialdata proseg_outputs.zarr \n
        --output-cell-polygons cell-polygons.geojson.gz \n
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \n
        --output-cell-metadata cell-metadata.csv.gz \n
        """
    for _key, _val in proseg_params.items():
        if _key not in proseg_args:
            logger.info("skipping unknown proseg parameter: %s" % _key)
            continue
        if isinstance(_val, bool):
            if _val:
                logger.info(f"Adding boolean proseg parameter: --{_key}")
                proseg_run_cmd += f"--{_key} \n"
        else:
            logger.info(f"Adding proseg parameter: --{_key} {_val}")
            proseg_run_cmd += f"--{_key} {_val} \n"

    # proseg_run_cmd += f"--nthreads {os.cpu_count()} \\\n"
    logger.info(f"Running proseg with command:\n{proseg_run_cmd}")

    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret



def run_proseg(root_dir: str, output_dir: str, region: str, rust_bin_path=None, **proseg_params):
    """
    Run the proseg command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    """

    for key, value in proseg_params.items():
        logger.info(f"Proseg parameter: {key} = {value}")

    # Ensure the output directory exists
    pathlib.Path(f"{output_dir}").mkdir(parents=True, exist_ok=True)

    # Add proseg path to environment
    v3_flag = _add_proseg_binary(rust_bin_path=rust_bin_path)

    # Prepare the command to run proseg
    if v3_flag: 
        ret = _execute_cli_proseg_v3(root_dir, output_dir, region, **proseg_params)
    else:
        raise RuntimeError("proseg 2.x can no longer be run (runner archived 2026-09-28); install proseg 3. "
                           "Loading existing v2 output is still supported.")

    print(ret.returncode)
    print(ret.stdout.decode("utf-8"))


def _proseg_3d_to_vpt_parquet(shapes_3d: gpd.GeoDataFrame, micron_per_z: float, output_path: str):
    shapes_3d = shapes_3d.rename_geometry("Geometry")
    shapes_3d['z'] = shapes_3d['layer']
    shapes_3d['ID'] = shapes_3d['cell']
    shapes_3d['ZIndex'] = shapes_3d['layer']
    # 1-based ZLevel (every plane owns a full micron_per_z slab), consistent with
    # cellpose/convert-geometry; column name standardized to "ZLevel".
    shapes_3d['ZLevel'] = (shapes_3d['ZIndex'] + 1) * micron_per_z
    shapes_3d['EntityID'] = pd.factorize(shapes_3d['ID'])[0] + 1
    shapes_3d.to_parquet(output_path, index=False)

def add_signals_meta_to_proseg(
    raw_dir: str | pathlib.Path,
    proseg_dir: str | pathlib.Path,
    shapes_filename: str = "cell-polygons-layers.geojson.gz",
    cell_meta_filename: str = "cell-metadata.csv.gz",
    micron_per_z: float = 1.5,
    vpt_bin_path: str | pathlib.Path | None = None,
    **kwargs
):
    """
    Add the sum signals metadata to the proseg cell metadata output.
    ``raw_dir`` is the raw MERSCOPE region dir; ``proseg_dir`` this region's proseg output dir.
    This is necessary for downstream analysis and visualization in VPT.
    """
    # from spida.S.segmentation.vpt_utils import read_micron_to_mosaic_transform, transform_geoms, update_geometry
    # import geopandas as gpd
    # import pandas as pd
    from spida.S.segmentation.vpt import _add_vpt_binary, _cli_sum_signals    
    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    if isinstance(proseg_dir, str):
        proseg_dir = pathlib.Path(proseg_dir)

    # Get paths and make sure they exist
    shapes_3d_path = proseg_dir / shapes_filename
    cell_meta_path = proseg_dir / cell_meta_filename
    assert shapes_3d_path.exists(), f"Shapes file not found at {shapes_3d_path}"
    assert cell_meta_path.exists(), f"Cell metadata file not found at {cell_meta_path}"
    # Read in the shapes and cell metadata
    with gzip.open(shapes_3d_path, "rb") as f:
        shapes_3d = gpd.read_file(f)
    logger.info(f"Read shapes from {shapes_3d_path} with shape {shapes_3d.shape}")
    with gzip.open(cell_meta_path, "rt") as f:
        cell_meta = pd.read_csv(f)
    logger.info(f"Read cell metadata from {cell_meta_path} with shape {cell_meta.shape}")

    # Convert proseg 3D shapes to vpt parquet format
    _proseg_3d_to_vpt_parquet(
        shapes_3d=shapes_3d,
        micron_per_z=micron_per_z,
        output_path= proseg_dir / "proseg_polygons.parquet"
    )
    logger.info(f"Converted proseg 3D shapes to vpt parquet format at {proseg_dir / 'proseg_polygons.parquet'}")
    
    # Sum signals 
    signals = _cli_sum_signals(
       raw_dir = raw_dir,
       seg_dir = proseg_dir,
       input_boundaries = "proseg_polygons.parquet",
       output_signals = "sum_signals.csv",
    )

    # Merge the signals with the cell metadata and save
    signals.index.name = 'EntityID'
    if (cell_meta.index.astype(int) == 0).any(): 
        signals.index = signals.index.astype(int) - 1
    cell_meta = cell_meta.merge(signals, left_on="cell", right_index=True).copy()
    cell_meta.to_csv(proseg_dir / "cell_metadata_with_signals.csv", index=False)
    logger.info(f"Merged signals with cell metadata and saved to {proseg_dir / 'cell_metadata_with_signals.csv'}")


