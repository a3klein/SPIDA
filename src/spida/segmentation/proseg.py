import os
from dotenv import load_dotenv # type: ignore
load_dotenv()

import pathlib
import subprocess


def _add_proseg_binary(): 

    """
    Add the path to the proseg binary to the PATH environment variable.
    This is necessary to ensure that the proseg command can be found when running the script.
    """
    # Need to add the rust path to the PATH variable
    os.environ['PATH'] += os.pathsep + os.getenv("RUST_BIN_PATH")

    #Test that proseg can be run: 
    ret = subprocess.run(["proseg", "-h"], capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("Proseg failed to run. Check the installation and the PATH variable.")
    
    print("proseg binary added to PATH successfully.")

def _execute_cli_proseg(root_dir:str, output_dir:str, region:str): 
    """
    Run the proseg command with the prespecified parameters using subprocess. 
    """

    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \
        --output-path {output_dir}/{region} \
        --output-expected-counts expected-counts.csv.gz \
        --output-cell-metadata cell-metadata.csv.gz \
        --output-transcript-metadata transcript-metadata.csv.gz \
        --output-cell-polygons cell-polygons.geojson.gz \
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \
        --output-union-cell-polygons union-cell-polygons.geojson.gz \
        --voxel-layers 4 \
        --nthreads 64 \
        """
    
    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret
    

def run_proseg(root_dir:str, output_dir:str, region:str):
    """
    Run the proseg command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    """
    
    # Ensure the output directory exists
    pathlib.Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)
    
    # Add proseg path to environment
    _add_proseg_binary()

    # Prepare the command to run proseg
    ret = _execute_cli_proseg(root_dir, output_dir, region)

    print(ret.returncode)
    print(ret.stdout.decode('utf-8'))



