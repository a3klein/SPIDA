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
        run_proseg(input_dir, output_dir, reg_name, **kwargs)
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

def align_proseg(
        exp_name:str, 
        reg_name:str, 
        seed_prefix_name:str="default",
        prefix_name:str="proseg",
        input_dir:str|Path = None,
        seg_dir:str|Path = None, 
        x : str = "x",
        y : str = "y",
        z : str = "global_z",
        cell_column : str = "cell_id",
        barcode_column : str = "barcode_id",
        gene_column : str = "gene",
        fov_column : str = "fov",
        cell_missing : int = '-1',
        min_jaccard : float = 0.4,
        min_prob : float = 0.5,
        filter_blank : bool = False,
        cell_metadata_fname:str="merged_cell_metadata.csv",
        cell_by_gene_fname:str="merged_cell_by_gene.csv",
        detected_transcripts_fname:str="merged_transcript_metadata.csv",
        cell_polygons_fname:str="merged_cell_polygons.geojson",
        ):
    """
    Align Proseg transcripts to seed transcripts.
    """

    from proseg import align_proseg_transcripts
    from vpt import generate_metadata
    from spida.io import load_segmentation_region

    ### TODO: 
    # 1. Get the segmentation directory >>
    # 1. Call the align_proseg_transcripts function from the proseg module.
    # 2. Call the generate_metadata function from the vpt module.
    # 3. Call the load_vpt_segmentation function from the io module.
    # I am expecting issues with loading the transcript and table spatialdata_io build in functions

    
    if input_dir is None: 
        processed_root_path =  os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"

    if seg_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        seg_dir = Path(f"{seg_out_path}/{exp_name}/{prefix_name}")
    if isinstance(seg_dir, str):
        seg_dir = Path(seg_dir)

    print (f"Aligning Proseg transcripts for region {reg_name} in experiment {exp_name}.")
    # Aligning the proseg transcripts to the seed transcripts
    align_proseg_transcripts(
        exp_name=exp_name,
        reg_name=reg_name,
        seed_prefix_name=seed_prefix_name,
        prefix_name=prefix_name,
        x=x,
        y=y,
        z=z,
        cell_column=cell_column,
        barcode_column=barcode_column,
        gene_column=gene_column,
        fov_column=fov_column,
        cell_missing=cell_missing,
        min_jaccard=min_jaccard,
        min_prob=min_prob,
        filter_blank=filter_blank,
        merged_transcripts_fname= detected_transcripts_fname,
        merged_cell_by_gene_fname= cell_by_gene_fname,
        merged_cell_polygons_fname= cell_polygons_fname,
        save_dir=seg_dir / reg_name,
    )

    print (f"Generating metadata for aligned transcripts in region {reg_name} of experiment {exp_name}.")
    # Generating the metadata for the aligned transcripts
    generate_metadata(
        root_dir=input_dir,
        seg_out_dir=seg_dir,
        region=reg_name,
        input_boundaries=cell_polygons_fname,
        output_boundaries="merged_converted_boundaries.parquet",
        input_entity_by_gene=cell_by_gene_fname,
        output_metadata=cell_metadata_fname,
        output_signals="merged_sum_signals.csv"
    )

    print (f"Loading segmentation data for region {reg_name} in experiment {exp_name}.")
    # Loading the new segmentation data into the spatialdata object 
    load_segmentation_region(
        exp_name=exp_name,
        reg_name=reg_name,
        seg_dir=seg_dir,
        type="vpt",
        prefix_name="prefix_name",
        plot=False,
        cell_metadata_fname=cell_metadata_fname,
        cell_by_gene_fname=cell_by_gene_fname,
        detected_transcripts_fname=detected_transcripts_fname,
        cellpose_micron_space_fname=cell_polygons_fname,
        )


if __name__ == "__main__":
    fire.Fire({"run_segmentation": run_segmentation, 
               "segment_experiment": segment_experiment,
               "align_proseg": align_proseg,})