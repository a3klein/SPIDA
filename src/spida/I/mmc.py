
import os 
from dotenv import load_dotenv  # type: ignore
load_dotenv()

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from datetime import datetime
date = datetime.today().strftime("%Y%m%d")


def _write_empty_adata(ref_path:str, gene_list:list, qout_path:str): 
    """
    Write an empty anndata object with the intersection of the reference and query gene names.
    
    Parameters:
    ref_path (str): Path to the reference anndata file.
    gene_list (list): List of gene names included in the query anndata object.
    qout_path (str): Path to save the empty query anndata object.
    """
    # Get the reference index 
    adata = ad.read_h5ad(ref_path, backed='r')
    var_index = adata.var.index

    # Get the query index list (from codebook)
    query_index = var_index.intersection(pd.Index(gene_list))
    # make the anndata object (query)
    query = ad.AnnData(
        X=sparse.csr_matrix(np.zeros((1, len(query_index)))),
        var=pd.DataFrame(index=query_index),
    )

    query.write_h5ad(qout_path)

def _setup_mmc(ref_path:str,
              mmc_store_path:str, 
              BRAIN_REGION:str,
              CODEBOOK:str,
              codebook_path:str,
              heirarchy_list:list,
              ref_norm:str="log2CPM", 
              **kwargs): 
    
    out_dir = f"{mmc_store_path}/{BRAIN_REGION}-{CODEBOOK}" ### An entry in the reference file table
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Intermediate filepaths 
    marker_dir = f"{out_dir}/marker_dir" # mmc_data_dir / store path
    query_path = f"{out_dir}/query.h5ad"

    gene_list = pd.read_csv(codebook_path)['name']
    gene_list = gene_list[~gene_list.str.startswith("Blank-")]
    _write_empty_adata(ref_path, gene_list, query_path) 
        
    ### PARAMETERS:  --> KWARGS? 
    n_cpu = kwargs.get('n_cpu', 8)
    n_per_utility = kwargs.get('n_per_utility', 10)
    n_valid = kwargs.get('n_valid', 10)
    
    from cell_type_mapper.cli.precompute_stats_scrattch import ( # type: ignore
            PrecomputationScrattchRunner)
    from cell_type_mapper.cli.reference_markers import ( # type: ignore
        ReferenceMarkerRunner
    )
    from cell_type_mapper.cli.query_markers import ( # type: ignore
        QueryMarkerRunner
    )


    ### Precompute the stats for the reference adata files        
    print('\n====starting precomputed_stats====\n')
    _temp_path = ref_path.split("/")[-1].split(".")[0]
    precomp_output_path = f"{out_dir}/{_temp_path}_precomputed_stats.h5"
    print(precomp_output_path)

    # doing the actual calculations
    precomputation_config = {
        'hierarchy': heirarchy_list,
        'h5ad_path': ref_path,
        'output_path': precomp_output_path,
        'n_processors': n_cpu,
        'normalization': ref_norm,
        'clobber': True,
    }
    precomputation_runner = PrecomputationScrattchRunner(
        args=[], input_data=precomputation_config)

    precomputation_runner.run()
    print('\n====done with precomputed_stats====\n')


    # calculating reference markers    
    print('\n====starting reference_markers====\n')
    reference_config = {
        'precomputed_path_list': [precomp_output_path],
        'n_valid': n_valid,
        'output_dir': marker_dir,
        'clobber': True,
        'query_path' : query_path, # TODO 
    }

    reference_runner = ReferenceMarkerRunner(
        args=[], input_data=reference_config)

    reference_runner.run()
    print('\n====done with reference_markers====\n')


    # Downsampling to query markers 
    print("======starting query markers======\n")    
    query_config = {
        'output_path': f'{marker_dir}/calc_markers.json',
        'reference_marker_path_list': [f"{marker_dir}/reference_markers.h5"],
        'n_per_utility': n_per_utility,
        'n_processors': n_cpu,
        'query_path' : query_path
    }
    query_runner = QueryMarkerRunner(
        args=[], input_data=query_config)
    query_runner.run()
    print("======done with query marker downsampling======\n")

    ### TODO: 
    # Save metadata about the map my cells annotation store into a metadtata LIMS store. 
    # Not necessary? 

    
def _mmc_runner(mmc_store_path:str, 
                 BRAIN_REGION:str, 
                 CODEBOOK:str,
                 q_path:ad.AnnData, 
                 output_path:str, 
                 q_norm:str = "log2CPM",                 
                 **kwargs):

    from cell_type_mapper.cli.from_specified_markers import ( # type: ignore
        FromSpecifiedMarkersRunner
    )
    
    mmc_data_dir = f"{mmc_store_path}/{BRAIN_REGION}-{CODEBOOK}" ### An entry in the reference file table
    markers = f"{mmc_data_dir}/marker_dir/calc_markers.json"
    precomp_path = glob.glob(f"{mmc_data_dir}/*_precomputed_stats.h5")[0]

    # MMC params
    n_cpu = kwargs.get('n_cpu', 1)
    bootstrap_factor = kwargs.get('bootstrap_factor', 0.8)
    bootstrap_iteration = kwargs.get('bootstrap_iteration', 100)
    rng = kwargs.get('rng_seed', 13)

    extended_results_path = f"{output_path}/extended_results.json"
    csv_results_path = f"{output_path}/csv_results.csv"

    ### The actual Run
    print('\n====starting mapping====\n')
    config = {
        'precomputed_stats': {
            'path': precomp_path
        },
        'query_markers': {
            'serialized_lookup': markers
        },
        'type_assignment': {
            'bootstrap_factor': bootstrap_factor,
            'bootstrap_iteration': bootstrap_iteration,
            'rng_seed': rng,
            'n_processors': n_cpu,
            'normalization': q_norm
        },
        'query_path': q_path,
        'extended_result_path': extended_results_path,
        'csv_result_path': csv_results_path
    }

    mapping_runner = FromSpecifiedMarkersRunner(
        args=[], input_data=config)

    mapping_runner.run()
    print('\n====done with mapping====\n')

    return (extended_results_path, csv_results_path)  # Return paths for further processing or verification

def _transfer_labels(adata:ad.AnnData, 
                     results_path:str,): 

    # Read in results
    results_df = pd.read_csv(results_path, comment='#', index_col=0)

    # Getting the heirarchy list from the results_df columns
    heirarchy_list = []
    for c in results_df.columns:
        if c.endswith("bootstrapping_probability"):
            heirarchy_list.append(c.split("_bootstrapping_probability")[0])

    # Adding the heirarchy list to adata with preassigned names 
    for h in heirarchy_list:
        adata.obs[f"mmc_{h}_transfer_score"] = results_df[f"{h}_bootstrapping_probability"]
        adata.obs[f"mmc_{h}"] = results_df[f"{h}_name"]

    return adata

def setup_mmc(ref_path:str, 
              BRAIN_REGION:str,
              CODEBOOK:str,
              heirarchy_list:list,
              codebook_path:str = None,
              mmc_store_path:str = None, 
              ref_norm:str="log2CPM", 
              **kwargs
              ):              
    """
    Setup function for MapMyCells integration.
    """

    print("START")

    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
        if not mmc_store_path:
            raise ValueError("Please provide a mmc_store_path or set the MMC_DIR environment variable.")
    
    if not codebook_path: 
        gene_panel_path = os.getenv("GENE_PANEL_PATH")
        if not gene_panel_path:
            raise ValueError("Please provide a codebook path or set the GENE_PANEL_PATH environment variable.")
        codebook_path = glob.glob(f"{gene_panel_path}/*{CODEBOOK}*.csv")[0]

    print("ref_path:", ref_path)
    print("heirarchy_list:",  heirarchy_list)
    print("BRAIN_REGION:", BRAIN_REGION)
    print("CODEBOOK:", CODEBOOK)
    print("codebook_path:", codebook_path)
    print("mmc_store_path:", mmc_store_path)
    print("ref_norm:", ref_norm)

    _setup_mmc(ref_path = ref_path,
               mmc_store_path=mmc_store_path, 
               BRAIN_REGION=BRAIN_REGION,
               CODEBOOK=CODEBOOK,
               codebook_path=codebook_path,
               heirarchy_list=heirarchy_list,
               ref_norm="log2CPM", 
               kwargs=kwargs)

    return 0

def run_mmc(query_adata:ad.AnnData,
            BRAIN_REGION:str,
            CODEBOOK:str,
            mmc_store_path:str=None, 
            anndata_store_path:str=None, 
            annotations_store_path:str=None,
            **kwargs):
     
    ### TODO: 
    # 1. Write the query adata out in a temporary directory (/bican/data/aggregated?) 
    # 2. Determine the output path (/bican/data/annotated?)
    # 2. Run the annotation using _mmc_runner
    # 3. Transfer the labels using _transfer_labels 
    # (written onto the original adata object so it can be written back into spatialdata memory)

    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
        if not mmc_store_path:
            raise ValueError("Please provide a mmc_store_path or set the MMC_DIR environment variable.")
    
    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")
        if not anndata_store_path:
            raise ValueError("Please provide an anndata_store_path or set the ANNDATA_STORE_PATH environment variable.")
        
    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")
        if not annotations_store_path:
            raise ValueError("Please provide an annotations_store_path or set the ANNOTATIONS_STORE_PATH environment variable.")
        
    exp = query_adata.uns.get("experiment")
    seg_name = query_adata.uns.get("segmentation")
    donor = query_adata.uns.get("donor")
    query_path = f"{anndata_store_path}/{exp}/{seg_name}/adata_{donor}.h5ad"
    output_path = f"{annotations_store_path}/{exp}/{seg_name}/mmc/{donor}"
    Path(output_path).mkdir(parents=True, exist_ok=True)
    print("Writing query adata to temporary path: %s" % (query_path))
    Path(query_path).parent.mkdir(parents=True, exist_ok=True)
    query_adata.write_h5ad(query_path)
    
    print("Running MMC: ")
    print("mmc_store_path:", mmc_store_path)
    print("BRAIN_REGION:", BRAIN_REGION)
    print("CODEBOOK:", CODEBOOK)
    print("query path", query_path)
    print("output_path:", output_path)
    

    res_paths = _mmc_runner(mmc_store_path = mmc_store_path, 
                            BRAIN_REGION = BRAIN_REGION, 
                            CODEBOOK = CODEBOOK,
                            q_path = query_path, 
                            output_path = output_path, 
                            q_norm = "log2CPM",
                            kwargs = kwargs)

    print("Transferring labels to the adata object")
    pre_cols = query_adata.obs.columns
    query_adata = _transfer_labels(query_adata, res_paths[1])
    post_cols = query_adata.obs.columns
    print("Added columns to adata.obs: ", set(post_cols) - set(pre_cols))

    query_adata.write_h5ad(query_path)
    
    return 0



#### Runners 
def mmc_setup(ref_path:str,
            heirarchy_list:list,
            BRAIN_REGION:str,
            CODEBOOK:str,
            codebook_path:str=None,
            mmc_store_path:str=None,
            ref_norm:str="log2CPM", 
            **kwargs
            ): 
    """
    Setup function for MapMyCells integration.
    
    Parameters:
    ref_path (str): Path to the reference data.
    heirarchy_list (list): List of hierarchy levels for the annotation.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    codebook_path (str, optional): Path to the codebook file. Defaults to None.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    ref_norm (str, optional): Normalization method for the reference AnnData.X. Defaults to "log2CPM".
    **kwargs: Additional keyword arguments.
    """
    setup_mmc(ref_path=ref_path,
            BRAIN_REGION=BRAIN_REGION,
            CODEBOOK=CODEBOOK,
            codebook_path=codebook_path,
            heirarchy_list=heirarchy_list,
            mmc_store_path=mmc_store_path,
            ref_norm=ref_norm,
            **kwargs)
    """
    Setup function for MapMyCells integration.
    """

    print("DONE")
    return 0


def mmc_annotation_region(
                        exp_name:str, 
                        reg_name:str, 
                        prefix_name:str, 
                        BRAIN_REGION:str, 
                        CODEBOOK:str, 
                        mmc_store_path:str=None, 
                        anndata_store_path:str=None, 
                        annotations_store_path:str=None,
                        **kwargs
                        ):
    """
    Run MapMyCells annotation on a given experiment and region.
    
    Parameters: 
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    print("RUNNING MMC ANNOTATIONS, EXPERIMENT %s, REGION %s, PREFIX %s" %(exp_name, reg_name, prefix_name) )

    # # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")

    ### Need to get the adata_path from this function (it must have been written out beforehand)

    return run_mmc(
            adata,
            BRAIN_REGION,
            CODEBOOK,
            mmc_store_path,
            anndata_store_path, 
            annotations_store_path,
            kwargs=kwargs,
        )


def mmc_annotation_experiment(
        exp_name:str, 
        prefix_name:str, 
        BRAIN_REGION:str, 
        CODEBOOK:str, 
        mmc_store_path:str=None, 
        anndata_store_path:str=None, 
        annotations_store_path:str=None,
        **kwargs
        ):
    """
    Run MapMyCells annotation for an entire experiment.
    
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    
    for reg in regions: 
        reg_name = reg.split("/")[-1]
        mmc_annotation_region(exp_name, reg_name, prefix_name, 
                            BRAIN_REGION=BRAIN_REGION, CODEBOOK=CODEBOOK,
                            mmc_store_path=mmc_store_path,
                            anndata_store_path=anndata_store_path,
                            annotations_store_path=annotations_store_path,
                            **kwargs)