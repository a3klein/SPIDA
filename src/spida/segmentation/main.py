import os
import sys
from dotenv import load_dotenv # type: ignore
load_dotenv()
from pathlib import Path
import fire # type: ignore
import glob

# sys.path.append(os.getenv("SPIDA_PATH"))


# The controller that calls the specified segmentation algorithms 
def run_segmentation(type:str, exp_name:str, reg_name:str, input_dir:Path=None, output_dir:Path=None, **kwargs):
    """
    Run an implemented segmentation algorithm on a given region

    Parameters:
    type (str): Type of segmentation to run (e.g., "proseg", "vpt", "cellpose", "mesmer").
    input_dir (Path): Directory containing the input data.
    output_dir (Path): Directory to save the output data.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    """
    
    if input_dir is None: 
        processed_root_path =  os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"
    if output_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        output_dir = f"{seg_out_path}/{exp_name}/{type}"

    # Printing out info # Implement logging here. 
    print(f"Running segmentation of type {type} on region {reg_name} in experiment {exp_name}.")
    print("Input directory:", input_dir)
    print("Output directory:", output_dir)

    if type == "proseg":
        from proseg import run_proseg
        run_proseg(input_dir, output_dir, reg_name)
    elif type == "vpt":
        from vpt import run_vpt
        config_path = kwargs.get("config_path", "/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json")
        run_vpt(input_dir, output_dir, reg_name, config_path=config_path)
    elif type == "cellpose":
        from cellposeSAM import run_cellposeSAM
        from vpt import seg_to_vpt
        run_cellposeSAM(input_dir, output_dir, reg_name, **kwargs)
        seg_to_vpt(input_dir, output_dir, reg_name)
    elif type == "mesmer":
        from mesmer import run_mesmer
        from vpt import seg_to_vpt
        run_mesmer(input_dir, output_dir, reg_name, **kwargs)
        seg_to_vpt(input_dir, output_dir, reg_name)

    else:
        raise ValueError(f"Unknown segmentation type: {type}")
    


def segment_experiment(type:str, exp_name:str, input_dir:Path=None, output_dir:Path=None, **kwargs): 
    """
    Run segmentation for all regions in an experiment.
    Parameters:
    type (str): Type of segmentation to run (e.g., "proseg", "vpt").
    exp_name (str): Name of the experiment.
    input_dir (Path, optional): Directory containing the input data. Defaults to None.
    output_dir (Path, optional): Directory to save the output data. Defaults to None.
    **kwargs: Additional keyword arguments to pass to the segmentation function.
    """

    root_path = os.getenv("ZARR_STORAGE_PATH")
    exp_path = Path(f"{root_path}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        run_segmentation(type=type, exp_name=exp_name, reg_name=reg_name,
                         input_dir=input_dir, output_dir=output_dir, **kwargs)



if __name__ == "__main__":
    fire.Fire({"run_segmentation": run_segmentation, 
               "segment_experiment": segment_experiment})