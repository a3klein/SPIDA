import os
import sys
from dotenv import load_dotenv # type: ignore
from pathlib import Path
import fire # type: ignore
import glob
import logging

load_dotenv()
logger = logging.getLogger(__package__)

# The controller that calls the specified segmentation algorithms 
def run_segmentation(type:str, exp_name:str, reg_name:str, input_dir:str|Path=None, output_dir:str|Path=None, **kwargs):
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
        processed_root_path = os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"
    if output_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        output_dir = f"{seg_out_path}/{exp_name}/{type}"

    # Printing out info # Implement logging here. 
    logger.info(f"Running segmentation of type {type} on region {reg_name} in experiment {exp_name}.")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    if type == "proseg":
        from .proseg import run_proseg
        run_proseg(input_dir, output_dir, reg_name, **kwargs)
    elif type == "vpt":
        from .vpt import run_vpt
        config_path = kwargs.get("config_path", "/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json")
        run_vpt(input_dir, output_dir, reg_name, config_path=config_path)
    elif type == "cellpose":
        from .cellposeSAM import run_cellposeSAM
        from .vpt import seg_to_vpt
        # TODO: change other segmentation methods to take in the direct image dir
        run_cellposeSAM(input_dir, output_dir, reg_name, **kwargs)
        input_dir_vpt = Path(input_dir).parents[1]
        seg_to_vpt(input_dir_vpt, output_dir, reg_name)
    elif type == "mesmer":
        from .mesmer import run_mesmer
        from .vpt import seg_to_vpt
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
        
def vpt_on_segmentation(exp_name:str, reg_name:str, input_dir:Path=None, output_dir:Path=None, **kwargs):
    """
    Convert segmentation results to VPT format.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    input_dir (Path, optional): Directory containing the input data. Defaults to None.
    output_dir (Path, optional): Directory to save the output data. Defaults to None.
    **kwargs: Additional keyword arguments to pass to the conversion function.
    """

    from vpt import seg_to_vpt
    if input_dir is None: 
        processed_root_path =  os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"
    if output_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        output_dir = f"{seg_out_path}/{exp_name}/vpt"
        
    seg_to_vpt(input_dir, output_dir, reg_name, **kwargs)

def align_proseg(
        exp_name:str, 
        reg_name:str, 
        seed_prefix_name:str="default",
        prefix_name:str="proseg",
        out_prefix_name:str="proseg_aligned",
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
        **kwargs
        ):
    """
    Align Proseg transcripts to seed transcripts.
    """

    from .proseg import align_proseg_transcripts
    from .vpt import generate_metadata, seg_to_vpt
    # from spida.S.io import load_segmentation_region

    if input_dir is None:
        processed_root_path = os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"

    if seg_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        seg_dir = Path(f"{seg_out_path}/{exp_name}/{prefix_name}")
    if isinstance(seg_dir, str):
        seg_dir = Path(seg_dir)

    logger.info(f"Aligning Proseg transcripts for region {reg_name} in experiment {exp_name}.")
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

    logger.info(f"Generating metadata for aligned transcripts in region {reg_name} of experiment {exp_name}.")
    # Generating the metadata for the aligned transcripts
    generate_metadata(
        root_dir=input_dir,
        seg_out_dir=seg_dir,
        region=reg_name,
        input_boundaries=cell_polygons_fname,
        output_boundaries="merged_converted_boundaries.parquet",
        input_transcripts=detected_transcripts_fname,
        input_entity_by_gene=cell_by_gene_fname,
        output_metadata=cell_metadata_fname,
        output_signals="merged_sum_signals.csv",
        output_entity_by_gene="cell_by_gene.csv", 
        output_transcripts="detected_transcripts.csv"
    )
    
    # generate_metadata(
    #     root_dir=input_dir,
    #     seg_out_dir=seg_dir,
    #     region=reg_name,
    #     input_boundaries=cell_polygons_fname,
    #     output_boundaries="merged_converted_boundaries.parquet",
    #     input_entity_by_gene=cell_by_gene_fname,
    #     output_metadata=cell_metadata_fname,
    #     output_signals="merged_sum_signals.csv"
    # )

    # logging.info(f"Loading segmentation data for region {reg_name} in experiment {exp_name}.")
    # # Loading the new segmentation data into the spatialdata object 
    # load_segmentation_region(
    #     exp_name=exp_name,
    #     reg_name=reg_name,
    #     seg_dir=seg_dir,
    #     type="vpt",
    #     prefix_name=out_prefix_name,
    #     plot=False,
    #     cell_metadata_fname=cell_metadata_fname,
    #     cell_by_gene_fname="cell_by_gene.csv", 
    #     detected_transcripts_fname="detected_transcripts.csv",
    #     cellpose_micron_space_fname="merged_converted_boundaries.parquet",
    #     )


if __name__ == "__main__":
    fire.Fire({"run_segmentation": run_segmentation, 
               "segment_experiment": segment_experiment,
               "align_proseg": align_proseg,})
