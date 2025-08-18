import os
import glob
from pathlib import Path
import logging
import warnings
from dotenv import load_dotenv

# assuming that all genes in the merscope dataset should overlap with the genes in the scRNA dataset.
import numpy as np
import pandas as pd
import anndata as ad
from ALLCools.clustering import significant_pc_test  # type : ignore
from ALLCools.integration.seurat_class import SeuratIntegration  # type: ignore
from sklearn.decomposition import TruncatedSVD
import time

load_dotenv()
logger = logging.getLogger(__package__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

def _downsample_reference(
    ref_adata : ad.AnnData,
    cluster_col : str,
    max_cluster_size: int = 3000,
    min_cluster_size: int = 0,
):
    """
    Remove clusters from the reference that have less than min_cluster_size cells.
    Downsample larger clusters that have more than max_cluster_size cells.
    """
    from ..utilities._ad_utils import _downsample_ref_clusters, _remove_small_clusters
    if min_cluster_size > 0: 
        ref_adata = _remove_small_clusters(ref_adata, cluster_col, min_cells=min_cluster_size)
    if max_cluster_size > 0:
        ref_adata = _downsample_ref_clusters(ref_adata, cluster_col, max_cells=max_cluster_size)
    return ref_adata

def _get_joined_deg_list(
    ref_adata, 
    ref_col_level,
    qry_adata,
    qry_col_level,
    top_n=20,
    min_cells=25,
    logfc_threshold=0.25,
    pval_threshold=0.05,
    method='wilcoxon',
    correction_method='benjamini-hochberg',
    verbose=1,
): 
    """
    Running DEG analysis on both the reference and query datasets to get a joint list of genes to use for integration.
    """
    from ..utilities._degs import call_degs_by_celltype, summarize_deg_results


    verbosity = False
    inner_verbose = False
    if verbose >= 2: 
        inner_verbose = True
    if verbose >= 1: 
        verbosity = True

    deg_ref = call_degs_by_celltype(
        adata=ref_adata,
        celltype_col=ref_col_level,
        min_cells=min_cells,
        logfc_threshold=logfc_threshold,
        pval_threshold=pval_threshold,
        method=method,
        correction_method=correction_method,
        save_results=False,
        verbose=inner_verbose
    )  

    summary_ref = summarize_deg_results(deg_ref, top_n=top_n)

    up_genes = np.unique([item for sublist in summary_ref['top_upregulated'] for item in sublist])
    down_genes = np.unique([item for sublist in summary_ref['top_downregulated'] for item in sublist])
    ref_genes = np.unique(np.concatenate((up_genes, down_genes)))
    if verbosity: logger.info(f"Total unique genes in top {top_n} reference: {len(ref_genes)}")

    deg_qry = call_degs_by_celltype(
        adata=qry_adata,
        celltype_col=qry_col_level,
        min_cells=min_cells,
        logfc_threshold=logfc_threshold,
        pval_threshold=pval_threshold,
        method=method,
        correction_method=correction_method,
        save_results=False,
        verbose=inner_verbose
    )  

    deg_qry = summarize_deg_results(deg_qry, top_n=top_n)

    up_genes = np.unique([item for sublist in deg_qry['top_upregulated'] for item in sublist])
    down_genes = np.unique([item for sublist in deg_qry['top_downregulated'] for item in sublist])
    qry_genes = np.unique(np.concatenate((up_genes, down_genes)))
    if verbosity: logger.info(f"Total unique genes in top {top_n} query: {len(qry_genes)}")

    total_genes = np.unique(np.concatenate((ref_genes, qry_genes)))
    if verbosity: logger.info(f"Total unique genes in both datasets: {len(total_genes)}")
    return total_genes


def _allcools_pca(
    ref_adata: ad.AnnData,
    qry_adata: ad.AnnData,
    n_train_cell: int = 100000,
    chunk_size: int = 50000,
):
    """
    Perform PCA on the concatenated AnnData object.

    Parameters:
    ref_adata (AnnData): Reference AnnData object.
    qry_adata (AnnData): Query AnnData object.
    **kwargs: Additional keyword arguments for PCA.

    Returns:
    AnnData: PCA transformed AnnData object.
    """

    # run PCA - RNA
    np.random.seed(13)

    logger.info("Downsampling Reference for PCA")

    # select cells to fit the model
    train_cell = np.zeros(ref_adata.shape[0]).astype(bool)
    if ref_adata.shape[0] > n_train_cell:
        train_cell[
            np.random.choice(np.arange(ref_adata.shape[0]), n_train_cell, False)
        ] = True
    else:
        train_cell[:] = True
    ref_adata.obs["Train"] = train_cell.copy()

    ndim = min(100, ref_adata.obs["Train"].sum() - 1, ref_adata.shape[1] - 1)
    model = TruncatedSVD(n_components=ndim, algorithm="arpack", random_state=0)
    model.fit(
        ref_adata.X[ref_adata.obs["Train"].values]
    )  # fit the model on reference only
    sel_dim = model.singular_values_ != 0
    logger.info(f"Selecting number of PCA dimensions: {sel_dim.sum()}")

    logger.info("Transforming Reference for PCA")
    chunks = []
    for chunk_start in range(0, ref_adata.shape[0], chunk_size):
        chunks.append(
            model.transform(ref_adata.X[chunk_start : (chunk_start + chunk_size)])
        )

    ref_adata.obsm["X_pca"] = np.concatenate(chunks, axis=0)[:, sel_dim]
    ref_adata.obsm["X_pca"] /= model.singular_values_[
        sel_dim
    ]  # unit-variance “whitening”

    logger.info("Transforming Query for PCA")
    chunks = []
    for chunk_start in range(0, qry_adata.shape[0], chunk_size):
        chunks.append(
            model.transform(qry_adata.X[chunk_start : (chunk_start + chunk_size)])
        )

    qry_adata.obsm["X_pca"] = np.concatenate(chunks, axis=0)[:, sel_dim]
    qry_adata.obsm["X_pca"] /= model.singular_values_[sel_dim]

    return (ref_adata, qry_adata)


def _allcools_seurat_wrapper(ref_adata, qry_adata, **kwargs):
    """
    Wrapper for ALLCools Seurat integration.

    Parameters:
    ref_adata (AnnData): Reference AnnData object.
    qry_adata (AnnData): Query AnnData object.
    **kwargs: Additional keyword arguments for SeuratIntegration.

    Returns:
    AnnData: Integrated AnnData object.
    """

    logger.info("Concatenating reference and query")
    adata = ref_adata.concatenate(
        qry_adata,
        batch_categories=["ref", "query"],
        batch_key="Modality",
        index_unique=None,
    )

    # ncell = ref_adata.shape[0] + qry_adata.shape[0]
    ncc = significant_pc_test(ref_adata, p_cutoff=0.1, update=False, obsm="X_pca")
    ncc = min(
        50, ncc, ref_adata.shape[0] - 1, qry_adata.shape[0] - 1, ref_adata.shape[1] // 5
    )
    ncc = max(ncc, 5)
    logger.info(
        f"number of pcs used in ref_data pca: {ref_adata.obsm['X_pca'].shape[1]}"
    )
    npc = min([50, ncc + 10, ref_adata.shape[0] - 1, ref_adata.obsm["X_pca"].shape[1]])
    logger.info(
        f"npc: {npc}, ncc: {ncc}, ref_adata.shape[0]: {ref_adata.shape[0]}, qry_adata.shape[0]: {qry_adata.shape[0]}"
    )

    logger.info("Running Seurat Integration on reference and query data")
    integrator = SeuratIntegration()
    adata_list = [ref_adata, qry_adata]

    start_time = time.time()
    integrator.find_anchor(
        adata_list,
        k_local=None,
        key_local="X_pca",
        k_anchor=5,
        key_anchor="X",  # where the data lives
        dim_red="cca",
        max_cc_cells=50000,
        k_score=30,
        # k_filter=min(200, ref_adata.shape[0] // 10),
        k_filter=None,
        scale1=True,
        scale2=True,
        # scale =[False, True]
        n_components=ncc,
        alignments=[[[0], [1]]],
    )

    start_time = time.time()
    corrected = integrator.integrate(
        key_correct="X_pca",
        row_normalize=True,
        n_components=npc,
        k_weight=min(100, integrator.anchor[(0, 1)].shape[0]),
        sd=1,
        alignments=[[[0], [1]]],
    )
    logger.info(f"Integration time: {time.time() - start_time}")

    adata.obsm["X_pca_integrate"] = np.concatenate(corrected)

    return adata, integrator


def _transfer_labels(adata, integrator, rna_cell_type="supercluster_term"):
    """
    Transfer labels from the reference to the query AnnData object using the Seurat integration results.

    Parameters:
    adata (AnnData): Integrated AnnData object containing both reference and query data.
    integrator (SeuratIntegration): The Seurat integration object containing the label transfer results.
    rna_cell_type (str): The column name in adata.obs that contains the RNA cell type labels.
    """

    logger.info("Running Label Transfer")
    key_to_transfer = [rna_cell_type]
    label_transfer = integrator.label_transfer(
        ref=[0],
        qry=[1],
        categorical_key=key_to_transfer,
        continuous_key=None,
        key_dist="X_pca",
    )
    labels = label_transfer[rna_cell_type].columns.tolist()
    adata.obs[f"allcools_{rna_cell_type}"] = adata.obs.index.to_series().map(
        label_transfer[rna_cell_type]
        .apply(lambda x: labels[np.argmax(x)], axis=1)
        .to_dict()
    )
    adata.obs[f"allcools_{rna_cell_type}_transfer_score"] = (
        adata.obs.index.to_series().map(
            label_transfer[rna_cell_type]
            .apply(lambda x: x[np.argmax(x)], axis=1)
            .to_dict()
        )
    )
    adata.obs["all_annot"] = [
        adata.obs.loc[x, rna_cell_type]
        if pd.isnull(adata.obs.loc[x, f"allcools_{rna_cell_type}"])
        else adata.obs.loc[x, f"allcools_{rna_cell_type}"]
        for x in adata.obs.index
    ]
    return adata


def _joint_embeddings(
    adata : ad.AnnData,
    use_rep: str = "X_pca_integrate",
    key_added: str = "integrated_",
    min_dist:float = 0.25,
    leiden_res:int = 1,
    knn:int = 50,
): 
    """
    For calculating the joint embeddings from the integrated PCA space
    """
    from ..P.setup_adata import _calc_embeddings
    return _calc_embeddings(adata, use_rep=use_rep, key_added=key_added, leiden_res=leiden_res, knn=knn, min_dist=min_dist)

def run_allcools_seurat(
    ref_adata: ad.AnnData,
    qry_adata: ad.AnnData,
    anndata_store_path: str,
    annotations_store_path: str,
    rna_cell_type_column: str = "supercluster_name",
    qry_cluster_column:str = "leiden",
    top_deg_genes:int = -1,
    max_cells_per_cluster:int = 3000,
    min_cells_per_cluster:int = 0,
    save_integrator:bool = True, 
    save_adata_comb: bool = True,
    run_joint_embeddings: bool = False,
    **kwargs,
):
    """
    Run ALLCools Seurat integration on two AnnData objects.

    Parameters:
    ref_adata (AnnData): Reference AnnData object.
    qry_adata (AnnData): Query AnnData object.
    **kwargs: Additional keyword arguments for SeuratIntegration.
    """

    logger.info("in the function run_allcools_seurat")

    exp = qry_adata.uns.get("experiment")
    seg_name = qry_adata.uns.get("segmentation")
    donor = qry_adata.uns.get("donor")
    qry_path = Path(f"{anndata_store_path}/{exp}/{seg_name}/adata_{donor}.h5ad")

    # All integration parameters defined here
    n_train_cell = kwargs.get("n_train_cell", 100000)
    chunk_size = kwargs.get("chunk_size", 50000)
    logger.info(f"RNA Cell Type Column: {rna_cell_type_column}")

    ### Downsample reference by clusters:
    ref_adata = _downsample_reference(
        ref_adata,
        cluster_col = rna_cell_type_column,
        max_cluster_size=max_cells_per_cluster,
        min_cluster_size=min_cells_per_cluster
    )

    ### Downsample genes to shared gene space:
    shared_genes = ref_adata.var.index.intersection(qry_adata.var.index)
    logger.info(f"Shared genes between reference and query: {len(shared_genes)}")
    ref_adata = ref_adata[:, shared_genes].copy()
    qry_adata_ds = qry_adata[:, shared_genes].copy()

    # TODO: call DEGs between clusters / cell types in ref and qry
    if top_deg_genes > 0:
        deg_genes =_get_joined_deg_list(
            ref_adata,
            rna_cell_type_column,
            qry_adata_ds,
            qry_cluster_column,
            top_n=top_deg_genes
        )
        logger.info(f"Using {len(deg_genes)} DEGs for integration")
        ref_adata = ref_adata[:, ref_adata.var_names.isin(deg_genes)].copy()
        qry_adata_ds = qry_adata_ds[:, qry_adata_ds.var_names.isin(deg_genes)].copy()

    # Run PCA on the reference and query data
    ref_adata, qry_adata_ds = _allcools_pca(
        ref_adata, qry_adata_ds, n_train_cell=n_train_cell, chunk_size=chunk_size
    )

    # Get the Seurat Integration object
    adata_comb, integrator = _allcools_seurat_wrapper(ref_adata, qry_adata_ds, **kwargs)

    # Transfer labels to adata:
    adata_comb = _transfer_labels(adata_comb, integrator, rna_cell_type=rna_cell_type_column)

    pre_cols = set(qry_adata_ds.obs.columns).union(ref_adata.obs.columns)
    post_cols = adata_comb.obs.columns
    cols_to_add = list(set(post_cols) - set(pre_cols))
    logger.info(f"Added columns to adata.obs: {cols_to_add}")
    qry_adata.obs[cols_to_add] = adata_comb.obs[cols_to_add].copy()

    qry_path.parent.mkdir(parents=True, exist_ok=True)
    qry_adata.write_h5ad(qry_path)

    if save_integrator:
        ### Making sure that there are no category dtypes in the integrator (they have nan's)
        for k, v in integrator.adata_dict.items():
            for col, val in v.obs.items():
                if val.dtype == "O":
                    logger.info(f"O: {col}")
                elif val.dtype == "category":
                    v.obs[col] = v.obs[col].cat.add_categories(["nan"])
                    v.obs[col] = v.obs[col].fillna("nan").astype(str)
        # Save the integrator object
        integrator_path = f"{annotations_store_path}/{exp}/{seg_name}/allcools/{donor}"
        Path(integrator_path).mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving integrator to {integrator_path}/integrator_rna_merfish")
        integrator.save(f"{integrator_path}/integrator_rna_merfish", save_adata=True)

    if save_adata_comb:
        adata_comb_path = f"{annotations_store_path}/{exp}/{seg_name}/allcools/{donor}"
        Path(adata_comb_path).mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving combined adata to {adata_comb_path}/adata_comb_rna_merfish.h5ad")
        adata_comb.write_h5ad(f"{adata_comb_path}/adata_comb_rna_merfish.h5ad")

    if run_joint_embeddings and save_adata_comb: 
        logger.info("Calculating joint embeddings")
        _joint_embeddings(adata_comb, use_rep="X_pca_integrate", key_added="integrated_", **kwargs)
        adata_comb.write_h5ad(f"{adata_comb_path}/adata_comb_rna_merfish.h5ad")

    return qry_adata


# ### Runners
# def allcools_integration_region(
#     exp_name: str,
#     reg_name: str,
#     prefix_name: str,
#     ref_path: str,
#     suffix:str = "_filt",
#     anndata_store_path: str = None,
#     annotations_store_path: str = None,
#     **kwargs,
# ):
#     """
#     Run ALLCools integration on a given experiment and region.
#     Parameters:
#     exp_name (str): Name of the experiment.
#     reg_name (str): Name of the region.
#     prefix_name (str): Prefix for the keys in the spatialdata object.
#     ref_path (str): Path to the reference RNA AnnData object .
#     anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
#     annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
#     to None.
#     **kwargs: Additional keyword arguments for ALLCools integration.
#     """
#     logger.info(
#         "RUNNING ALLCOOLS INTEGRATION, EXPERIMENT %s, REGION %s, PREFIX %s"
#         % (exp_name, reg_name, prefix_name)
#     )
#     if kwargs is not None: 
#         kwargs = kwargs['kwargs']
#     logger.info(f"ALLCOOLS kwargs: {kwargs}")

#     # # Getting the sdata object (right now from a constant zarr store path)
#     zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
#     logger.info(f"reading data, zarr_store {zarr_store}")
#     zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
#     logger.info(f"zarr_path: {zarr_path}")
#     ad_path = f"{zarr_path}/tables/{prefix_name}_table{suffix}"
#     logger.info(f"ad_path: {ad_path}")
#     adata = ad.read_zarr(ad_path)
#     logger.info(f"adata shape: {adata.shape}")

#     ref_adata = ad.read_h5ad(ref_path)
#     logger.info(f"ref_adata shape: {ref_adata.shape}")

#     logger.info("read data and calling the function")
#     adata = run_allcools_seurat(
#         ref_adata, adata, anndata_store_path, annotations_store_path, **kwargs
#     )

#     logger.info("DONE WITH ALLCOOLS INTEGRATION")


# def allcools_integration_experiment(
#     exp_name: str,
#     prefix_name: str,
#     ref_path: str,
#     suffix:str = "_filt",
#     anndata_store_path: str = None,
#     annotations_store_path: str = None,
#     **kwargs,
# ):
#     """
#     Run ALLCools integration for an entire experiment.

#     Parameters:
#     exp_name (str): Name of the experiment.
#     prefix_name (str): Prefix for the keys in the spatialdata object.
#     ref_path (str): Path to the reference RNA AnnData object .
#     anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
#     annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
#     to None.
#     **kwargs: Additional keyword arguments for ALLCools integration.
#     """

#     # Getting the regions for the experiment
#     zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
#     exp_path = Path(f"{zarr_store}/{exp_name}")
#     regions = glob.glob(f"{exp_path}/region_*")

#     for reg in regions:
#         reg_name = reg.split("/")[-1]
#         allcools_integration_region(
#             exp_name,
#             reg_name,
#             prefix_name,
#             ref_path,
#             suffix=suffix,
#             anndata_store_path=anndata_store_path,
#             annotations_store_path=annotations_store_path,
#             **kwargs,
#         )
