import os
import glob
from pathlib import Path
import logging

# assuming that all genes in the merscope dataset should overlap with the genes in the scRNA dataset.
import numpy as np
import pandas as pd
import anndata as ad
from ALLCools.clustering import significant_pc_test  # type : ignore
from ALLCools.integration.seurat_class import SeuratIntegration  # type: ignore
from sklearn.decomposition import TruncatedSVD
import time


def _downsample_reference(ref_adata, max_cluster_size=100000, min_cluster_size=100):
    """
    Remove clusters from the reference that have less than min_cluster_size cells.
    Downsample larger clusters that have more than max_cluster_size cells.
    """
    return 0


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

    logging.info("Downsampling Reference for PCA")

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
    logging.info(f"Selecting number of PCA dimensions: {sel_dim.sum()}")

    logging.info("Transforming Reference for PCA")
    chunks = []
    for chunk_start in range(0, ref_adata.shape[0], chunk_size):
        chunks.append(
            model.transform(ref_adata.X[chunk_start : (chunk_start + chunk_size)])
        )

    ref_adata.obsm["X_pca"] = np.concatenate(chunks, axis=0)[:, sel_dim]
    ref_adata.obsm["X_pca"] /= model.singular_values_[
        sel_dim
    ]  # unit-variance “whitening”

    logging.info("Transforming Query for PCA")
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

    logging.info("Concatenating reference and query")
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
    logging.info(
        f"number of pcs used in ref_data pca: {ref_adata.obsm['X_pca'].shape[1]}"
    )
    npc = min([50, ncc + 10, ref_adata.shape[0] - 1, ref_adata.obsm["X_pca"].shape[1]])
    logging.info(
        f"npc: {npc}, ncc: {ncc}, ref_adata.shape[0]: {ref_adata.shape[0]}, qry_adata.shape[0]: {qry_adata.shape[0]}"
    )

    logging.info("Running Seurat Integration on reference and query data")
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
    logging.info(f"Integration time: {time.time() - start_time}")

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

    logging.info("Running Label Transfer")
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


def run_allcools_seurat(
    ref_adata: ad.AnnData,
    qry_adata: ad.AnnData,
    anndata_store_path: str,
    annotations_store_path: str,
    **kwargs,
):
    """
    Run ALLCools Seurat integration on two AnnData objects.

    Parameters:
    ref_adata (AnnData): Reference AnnData object.
    qry_adata (AnnData): Query AnnData object.
    **kwargs: Additional keyword arguments for SeuratIntegration.
    """

    # checking if this works
    pre_cols = set(qry_adata.obs.columns).union(ref_adata.obs.columns)

    # Loading global parameters
    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")
        if not anndata_store_path:
            raise ValueError(
                "Please provide an anndata_store_path or set the ANNDATA_STORE_PATH environment variable."
            )

    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")
        if not annotations_store_path:
            raise ValueError(
                "Please provide an annotations_store_path or set the ANNOTATIONS_STORE_PATH environment variable."
            )

    exp = qry_adata.uns.get("experiment")
    seg_name = qry_adata.uns.get("segmentation")
    donor = qry_adata.uns.get("donor")
    qry_path = f"{anndata_store_path}/{exp}/{seg_name}/adata_{donor}.h5ad"

    # All integration parameters defined here
    # max_cluster_size = kwargs.get('max_cluster_size', 100000)
    # min_cluster_size = kwargs.get('min_cluster_size', 100)
    n_train_cell = kwargs.get("n_train_cell", 100000)
    chunk_size = kwargs.get("chunk_size", 50000)
    rna_cell_type = kwargs.get("rna_cell_type_column", "supercluster_name")
    logging.info(f"RNA Cell Type Column: {rna_cell_type}")
    save_integrator = kwargs.get("save_integrator", True)
    save_adata_comb = kwargs.get("save_adata_comb", True)

    ### TODO: Do I want to implement a downsampling on the reference data?
    # adata_ref = _downsample_reference(ref_adata, max_cluster_size=max_cluster_size, min_cluster_size=min_cluster_size)
    ### Downsample genes to shared gene space:
    shared_genes = ref_adata.var.index.intersection(qry_adata.var.index)
    logging.info(f"Shared genes between reference and query: {len(shared_genes)}")
    ref_adata = ref_adata[:, shared_genes].copy()
    qry_adata = qry_adata[:, shared_genes].copy()

    # Run PCA on the reference and query data
    ref_adata, qry_adata = _allcools_pca(
        ref_adata, qry_adata, n_train_cell=n_train_cell, chunk_size=chunk_size
    )

    # Get the Seurat Integration object
    adata_comb, integrator = _allcools_seurat_wrapper(ref_adata, qry_adata, **kwargs)

    # Transfer labels to adata:
    adata_comb = _transfer_labels(adata_comb, integrator, rna_cell_type=rna_cell_type)

    pre_cols = set(qry_adata.obs.columns).union(ref_adata.obs.columns)
    post_cols = adata_comb.obs.columns
    cols_to_add = list(set(post_cols) - set(pre_cols))
    logging.info(f"Added columns to adata.obs: {cols_to_add}")
    qry_adata.obs[cols_to_add] = adata_comb.obs[cols_to_add].copy()

    qry_adata.write_h5ad(qry_path)

    if save_integrator:
        ### Making sure that there are no category dtypes in the integrator (they have nan's)
        for k, v in integrator.adata_dict.items():
            for col, val in v.obs.items():
                if val.dtype == "O":
                    logging.info(f"O: {col}")
                elif val.dtype == "category":
                    v.obs[col] = v.obs[col].cat.add_categories(["nan"])
                    v.obs[col] = v.obs[col].fillna("nan").astype(str)
        # Save the integrator object
        integrator_path = f"{annotations_store_path}/{exp}/{seg_name}/allcools/{donor}"
        Path(integrator_path).mkdir(parents=True, exist_ok=True)
        integrator.save(f"{integrator_path}/integrator_rna_merfish", save_adata=True)

    if save_adata_comb:
        adata_comb_path = f"{annotations_store_path}/{exp}/{seg_name}/allcools/{donor}"
        Path(adata_comb_path).mkdir(parents=True, exist_ok=True)
        adata_comb.write_h5ad(f"{adata_comb_path}/adata_comb_rna_merfish.h5ad")

    return qry_adata


### Runners
def allcools_integration_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    ref_path: str,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Run ALLCools integration on a given experiment and region.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    ref_path (str): Path to the reference RNA AnnData object .
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
    to None.
    **kwargs: Additional keyword arguments for ALLCools integration.
    """
    logging.info(
        "RUNNING ALLCOOLS INTEGRATION, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    # # Getting the sdata object (right now from a constant zarr store path)
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")

    ref_adata = ad.read_h5ad(ref_path)

    adata = run_allcools_seurat(
        ref_adata, adata, anndata_store_path, annotations_store_path, **kwargs
    )

    logging.info("DONE WITH ALLCOOLS INTEGRATION")


def allcools_integration_experiment(
    exp_name: str,
    prefix_name: str,
    ref_path: str,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Run ALLCools integration for an entire experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    ref_path (str): Path to the reference RNA AnnData object .
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
    to None.
    **kwargs: Additional keyword arguments for ALLCools integration.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    for reg in regions:
        reg_name = reg.split("/")[-1]
        allcools_integration_region(
            exp_name,
            reg_name,
            prefix_name,
            ref_path,
            anndata_store_path=anndata_store_path,
            annotations_store_path=annotations_store_path,
            **kwargs,
        )
