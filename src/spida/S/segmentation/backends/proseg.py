import os

import logging
import pathlib
import subprocess
import re

from spida.io.formats.segmenters import (
    PROSEG_CELL_METADATA, PROSEG_POLYGON_LAYERS, PROSEG_POLYGONS, PROSEG_SPATIALDATA,
)


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


def _execute_cli_proseg_v3(transcripts_path: str, output_dir: str, **proseg_params):
    """
    Run the proseg command for the 3.0.1 version of proseg
    """
    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {transcripts_path} \n
        --output-path {output_dir} \n
        --output-spatialdata {PROSEG_SPATIALDATA} \n
        --output-cell-polygons {PROSEG_POLYGONS} \n
        --output-cell-polygon-layers {PROSEG_POLYGON_LAYERS} \n
        --output-cell-metadata {PROSEG_CELL_METADATA} \n
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



def run_proseg(transcripts_path: str, output_dir: str, rust_bin_path=None, **proseg_params):
    """
    Run proseg on a MERSCOPE-format transcripts table whose ``cell_id`` column is the starting segmentation.

    Parameters:
    transcripts_path (str): The prior segmentation's ``detected_transcripts.csv`` (e.g. cellpose_cell's).
    output_dir (str): The directory where the output files will be saved.
    """

    for key, value in proseg_params.items():
        logger.info(f"Proseg parameter: {key} = {value}")

    # Ensure the output directory exists
    pathlib.Path(f"{output_dir}").mkdir(parents=True, exist_ok=True)

    # Add proseg path to environment
    v3_flag = _add_proseg_binary(rust_bin_path=rust_bin_path)

    # Prepare the command to run proseg
    if v3_flag: 
        ret = _execute_cli_proseg_v3(transcripts_path, output_dir, **proseg_params)
    else:
        raise RuntimeError("proseg 2.x can no longer be run (runner archived 2026-09-28); install proseg 3. "
                           "Loading existing v2 output is still supported.")

    print(ret.returncode)
    print(ret.stdout.decode("utf-8"))
