import os
import glob
from pathlib import Path
import tempfile
import logging
import warnings
from dotenv import load_dotenv

# assuming that all genes in the merscope dataset should overlap with the genes in the scRNA dataset.
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from ALLCools.clustering import significant_pc_test  # type: ignore
from ALLCools.integration import confusion_matrix_clustering # type: ignore
from ALLCools.integration.seurat_class import SeuratIntegration  # type: ignore
from sklearn.decomposition import TruncatedSVD

from ..utilities.ad_utils import normalize_adata
import time

load_dotenv()
logger = logging.getLogger(__name__)
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
    from ..utilities.ad_utils import _downsample_ref_clusters, _remove_small_clusters
    if min_cluster_size > 0: 
        ref_adata = _remove_small_clusters(ref_adata, cluster_col, min_cells=min_cluster_size)
    if max_cluster_size > 0:
        ref_adata = _downsample_ref_clusters(ref_adata, cluster_col, max_cells=max_cluster_size)
    return ref_adata

def _allcools_cef(
    ref_adata,
    qry_adata,
    ref_col_level: str = 'Subclass_label',
    qry_col_level: str = None,
    top_n: int = 200,
    alpha: float = 0.05,
):
    from ALLCools.clustering import cluster_enriched_features # type: ignore
    cluster_enriched_features(ref_adata,
                            cluster_col=ref_col_level,
                            top_n=top_n,
                            alpha=alpha,
                            stat_plot=False,
                            method="rna")
    ref_cef = ref_adata.var[f'{ref_col_level}_enriched_features']
    ref_cef = ref_cef[ref_cef].index

    if qry_col_level is not None and qry_col_level in qry_adata.obs.columns:
        cluster_enriched_features(qry_adata,
                                cluster_col=qry_col_level,
                                top_n=top_n,
                                alpha=alpha,
                                stat_plot=False,
                                method="rna")
        qry_cef = qry_adata.var[f'{qry_col_level}_enriched_features']
        qry_cef = qry_cef[qry_cef].index
        return ref_cef.union(qry_cef)

    return ref_cef

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
    from ..utilities.degs import (
        call_degs_scran, summarize_scran_results,
        call_degs_scanpy, summarize_scanpy_results,
    )

    verbosity = False
    inner_verbose = False
    if verbose >= 2:
        inner_verbose = True
    if verbose >= 1:
        verbosity = True

    def _collect_genes(adata, celltype_col, n):
        try:
            raw = call_degs_scran(
                adata, celltype_col=celltype_col,
                min_cells=min_cells, verbose=inner_verbose,
            )
            summary = summarize_scran_results(raw, n_genes=n)
            top_col = next(c for c in summary.columns if c.startswith('top_') and c.endswith('_markers'))
            return np.unique([g for genes in summary[top_col] for g in genes])
        except Exception as e:
            logger.warning(f"scran DEGs failed ({e}), falling back to scanpy")
            raw = call_degs_scanpy(
                adata, celltype_col=celltype_col,
                min_cells=min_cells, logfc_threshold=logfc_threshold,
                pval_threshold=pval_threshold, method=method,
                correction_method=correction_method,
                save_results=False, verbose=inner_verbose,
            )
            summary = summarize_scanpy_results(raw, top_n=n)
            up = np.unique([g for genes in summary['top_upregulated'] for g in genes])
            dn = np.unique([g for genes in summary['top_downregulated'] for g in genes])
            return np.unique(np.concatenate((up, dn)))

    try:
        ref_genes = _collect_genes(ref_adata, ref_col_level, top_n)
    except Exception as e:
        logger.warning(f"Reference DEGs failed entirely: {e}")
        ref_genes = np.array([])
    if verbosity: logger.info(f"Total unique genes in top {top_n} reference: {len(ref_genes)}")

    try:
        qry_genes = _collect_genes(qry_adata, qry_col_level, top_n)
    except Exception as e:
        logger.warning(f"Query DEGs failed entirely: {e}")
        qry_genes = np.array([])
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
    ref_adata.obs["Train"] = train_cell

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
    adata = ad.concat(
        [ref_adata, qry_adata],
        label="Modality",
        keys=["ref", "query"],
        join="outer",
        index_unique=None,
        uns_merge="same",
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
        scale1=False,
        scale2=False,
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


def _transfer_labels(
    adata,
    integrator,
    rna_cell_type="supercluster_term",
    label_transfer_k=100,
    save_transfer: str | Path | None = None,
):
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
        k_weight=label_transfer_k,
    )
    # IF WANTING TO SAVE THE TRANSFER RESULTS
    if save_transfer is not None:
        logger.info(f"Saving label transfer results to {save_transfer}")
        label_transfer[rna_cell_type].to_csv(save_transfer, sep="\t")
    score_df = label_transfer[rna_cell_type]   # (n_query × n_labels) DataFrame
    labels = score_df.columns.tolist()
    adata.obs[f"allcools_{rna_cell_type}"] = adata.obs.index.to_series().map(
        score_df.apply(lambda x: labels[np.argmax(x)], axis=1).to_dict()
    )
    adata.obs[f"allcools_{rna_cell_type}_transfer_score"] = (
        adata.obs.index.to_series().map(
            score_df.apply(lambda x: x[np.argmax(x)], axis=1).to_dict()
        )
    )
    allcools_col = adata.obs[f"allcools_{rna_cell_type}"]
    ref_col = adata.obs.get(rna_cell_type)
    if ref_col is not None:
        adata.obs["all_annot"] = np.where(allcools_col.isna(), ref_col, allcools_col)
    else:
        adata.obs["all_annot"] = allcools_col
    return adata, score_df

def filt_adata_allcools(
    adata: ad.AnnData,
    label_col: str,
    transfer_score_threshold: float = 0.5,
    score_col: str = None,
    filt_col: str = None,
) -> None:
    """Filter allcools annotation by transfer score, modifying adata.obs in place.

    Cells whose transfer score is below transfer_score_threshold are relabeled
    "unknown" in the output filt_col.  Cells with no score (NaN — typically ref
    cells in a combined object) are left unchanged.

    Parameters
    ----------
    adata : AnnData
        Must contain label_col and the corresponding score column.
    label_col : str
        Full label column name, e.g. "allcools_Subclass".
    transfer_score_threshold : float
        Cells with score < threshold are marked "unknown". Default 0.5.
    score_col : str, optional
        Score column; defaults to f"{label_col}_transfer_score".
    filt_col : str, optional
        Output column; defaults to f"{label_col}_filt".
    """
    if score_col is None:
        score_col = f"{label_col}_transfer_score"
    if filt_col is None:
        filt_col = f"{label_col}_filt"

    adata.obs[filt_col] = adata.obs[label_col].astype("category")
    if "unknown" not in adata.obs[filt_col].cat.categories:
        adata.obs[filt_col] = adata.obs[filt_col].cat.add_categories("unknown")

    low_quality = adata.obs.index[
        adata.obs[score_col].notna() & (adata.obs[score_col] < transfer_score_threshold)
    ]
    adata.obs.loc[low_quality, filt_col] = "unknown"
    n_scored = adata.obs[score_col].notna().sum()
    logger.info(
        f"filt_adata_allcools: {len(low_quality)} / {n_scored} scored cells marked "
        f"unknown (score < {transfer_score_threshold})"
    )


def _joint_embeddings(
    adata : ad.AnnData,
    use_rep: str = "X_pca_integrate",
    key_added: str = "integrated_",
    min_dist:float = 0.25,
    leiden_res:int = 1,
    knn:int = 50,
    target_accuracy=0.85,
    **kwargs
): 
    """
    For calculating the joint embeddings from the integrated PCA space
    """
    from ..utilities.ad_utils import _calc_embeddings
    return _calc_embeddings(
        adata,
        use_rep=use_rep,
        key_added=key_added,
        leiden_res=leiden_res,
        knn=knn,
        min_dist=min_dist,
        target_accuracy=target_accuracy,
        **kwargs
    )

def _clust2clust_transfer(
    adata_comb,
    ref_cell_type_column="Subclass",
    joint_cluster_column="integrated_leiden",
    ref_proportion=0.8,
    ref_cell_count=5,
    qry_cell_count=5,
    ref_batch_modality_key="ref",
    qry_batch_modality_key="query",
):
    """Proportion-threshold cluster-to-cluster label transfer.

    For each integrated Leiden cluster, assigns the dominant reference cell type
    if it clears ref_proportion and minimum cell count thresholds; otherwise marks
    as "unknown".  Mirrors the notebook inline c2c logic with explicit parameters.
    """
    ref_obs = adata_comb.obs.loc[adata_comb.obs["Modality"] == ref_batch_modality_key]
    qry_obs = adata_comb.obs.loc[adata_comb.obs["Modality"] == qry_batch_modality_key]

    ref_vc = (
        ref_obs.groupby(joint_cluster_column, observed=True)[ref_cell_type_column]
        .value_counts(normalize=True)
        .sort_values(ascending=False)
        .reset_index()
    )
    ref_vc.drop_duplicates(joint_cluster_column, keep="first", inplace=True)
    ref_vc.rename(columns={"proportion": f"{ref_cell_type_column}_proportion"}, inplace=True)
    df_map = ref_vc.set_index(joint_cluster_column)

    df_map[f"{ref_cell_type_column}_cell_count"] = (
        df_map.index.to_series()
        .map(ref_obs.groupby(joint_cluster_column, observed=True)[ref_cell_type_column].count())
        .astype(float)   # cast away Categorical before fillna so 0 is a valid fill value
        .fillna(0)
        .astype(int)
    )
    df_map["qry_cell_count"] = (
        df_map.index.to_series()
        .map(qry_obs.groupby(joint_cluster_column, observed=True).size())
        .astype(float)   # same fix
        .fillna(0)
        .astype(int)
    )

    added_col = f"c2c_allcools_label_{ref_cell_type_column}"
    adata_comb.uns[f"{joint_cluster_column}_{ref_cell_type_column}_map"] = df_map.copy()

    df_filtered = df_map.query(
        f"{ref_cell_type_column}_cell_count >= {ref_cell_count} & qry_cell_count >= {qry_cell_count}"
    )
    df_filtered = df_filtered.query(f"{ref_cell_type_column}_proportion >= {ref_proportion}")

    cluster_to_label = df_filtered[ref_cell_type_column].to_dict()
    adata_comb.obs[added_col] = (
        adata_comb.obs[joint_cluster_column]
        .map(cluster_to_label)
        .fillna("unknown")
        .astype("category")
    )

    logger.info(
        f"c2c proportion transfer: {len(cluster_to_label)} / {len(ref_vc)} clusters assigned "
        f"(ref_proportion>={ref_proportion}, ref_cell_count>={ref_cell_count}, qry_cell_count>={qry_cell_count})"
    )
    return adata_comb


def _allcools_clust2clust_transfer(
    adata_comb,
    ref_cell_type_column = "Subclass",
    joint_cluster_column = "integrated_leiden",
    confusion_matrix_cluster_min_value = 0.25,
    confusion_matrix_cluster_max_value = 0.9,
    confusion_matrix_cluster_resolution = 1.5,
    qry_only_cluster_threshold = 50,
    ref_only_cluster_threshold = 20,
    integration_round: int = 1,
    ref_batch_modality_key: str = "ref",
    qry_batch_modality_key: str = "query",
):
    data = adata_comb.obs.loc[~adata_comb.obs[joint_cluster_column].isna(), [joint_cluster_column, ref_cell_type_column]].value_counts().unstack(fill_value=0)
    data = data.drop(columns=['nan']) if 'nan' in data.columns else data

    merfish_only_cluster = data.index[data.sum(axis=1) < qry_only_cluster_threshold]
    rna_only_cluster = data.columns[data.sum(axis=0) < ref_only_cluster_threshold]

    data = data.drop(merfish_only_cluster, axis=0)
    datac = data / data.sum(axis=0)
    datar = data / data.sum(axis=1).values[:, None]
    confusion_matrix = datar.where(datar > datac, datac)

    logger.info(f"RNA only clusters: {rna_only_cluster.tolist()}, Merfish only clusters: {merfish_only_cluster.tolist()}")

    # confusion matrix contains RNA only clusters, but not merfish only clusters
    if data.values.max()>0.1:
        (query_group, ref_group, confusion_matrix, g, modularity_score1) = confusion_matrix_clustering(
            confusion_matrix=confusion_matrix.T,
            min_value=confusion_matrix_cluster_min_value,
            max_value=confusion_matrix_cluster_max_value,
            partition_type=None,
            resolution=confusion_matrix_cluster_resolution,
            seed=0,
        )
    else:
        ref_group = pd.Series(-1, index=data.index)
        query_group = pd.Series(-1, index=data.columns)
        confusion_matrix = confusion_matrix.T
        modularity_score1 = 0
    
    logger.info(f"Modularity score: {modularity_score1}, Ref group unique shapes: {ref_group.unique().shape}, Query group unique shapes: {query_group.unique().shape}")

    integration_group_cluster = {}
    integration_group_cell = {}
    for group,xx in ref_group.groupby(ref_group):
        qry_cluster = xx.index
        ref_cluster = (datac.loc[qry_cluster].sum(axis=0) > 0.3) | (datar.loc[qry_cluster].max(axis=0) > 0.3)  #parameterize?
        ref_cluster = ref_cluster.index[ref_cluster].astype(str)
        ref_cell = adata_comb.obs.index[(adata_comb.obs['Modality']==ref_batch_modality_key) & (adata_comb.obs[ref_cell_type_column].isin(ref_cluster))]
        qry_cell = adata_comb.obs.index[(adata_comb.obs['Modality']==qry_batch_modality_key) & (adata_comb.obs[joint_cluster_column].isin(qry_cluster))]
        if (len(qry_cell)>0) and (len(ref_cell)>0):
            ## integration_group is used for split coclusters into next iteration, so should have cells from both modalities
            ## qry_cell must already be positive since len(qry_cluster)>0 in the last condition, so this if requires the group to have MERFISH cells
            integration_group_cluster[f'IG{group}'] = {
                "ref": ref_cluster.tolist(),
                "qry": qry_cluster.tolist()
            }
            integration_group_cell[f'IG{group}'] = {
                "ref": ref_cell.tolist(),
                "qry": qry_cell.tolist()
            }
        else:
            logger.info(f"non_int_group: {group}, {ref_cluster}, {qry_cluster}")

    logger.info(f"Number of integration groups: {len(integration_group_cluster)}, {len(integration_group_cell)}")
    
    added_col = f"c2c_allcools_label_{ref_cell_type_column}"
    adata_comb.obs.loc[adata_comb.obs['Modality']==qry_batch_modality_key, added_col] = "unknown"
    for group in integration_group_cell: 
        ref_cells = integration_group_cell[group]['ref']
        qry_cells = integration_group_cell[group]['qry']
        ref_annot = integration_group_cluster[group]['ref']
        logger.info(f"{group}, {len(ref_cells)}, {len(qry_cells)}, {ref_annot}")
        if len(ref_annot) > 1: 
            if integration_round == 1: # Doing a second round of integration with only the ambiguous clusters
                anndata_dir = tempfile.TemporaryDirectory()
                _, adata_comb_ret = run_allcools_seurat(
                    ref_adata = adata_comb[ref_cells,:].copy(),
                    qry_adata = adata_comb[qry_cells,:].copy(),
                    anndata_store_path=anndata_dir,
                    annotations_store_path=None,
                    rna_cell_type_column = ref_cell_type_column,
                    qry_cluster_column = joint_cluster_column,
                    top_deg_genes = 20,
                    max_cells_per_cluster=2000,
                    min_cells_per_cluster=0,
                    label_transfer_k=20,
                    save_adata_comb=False,
                    save_integrator=False,
                    save_query=False,
                    run_joint_embeddings=True,
                    joint_embedding_leiden_res=0.5,
                    run_clust_label_transfer=True,
                    c2c_method="allcools",
                    integration_round=2,
                    confusion_matrix_cluster_resolution=1,
                    confusion_matrix_cluster_min_value=0.1,
                    qry_only_cluster_threshold=10,
                    ref_only_cluster_threshold=10,
                )
                adata_comb.obs[added_col] = adata_comb.obs[added_col].astype("str")
                adata_comb.obs.loc[qry_cells,added_col] = adata_comb_ret.obs.loc[qry_cells, added_col].astype("str").values
                adata_comb.obs[added_col] = adata_comb.obs[added_col].astype("category")
            else: #TODO: Decide on an R2 strategy
                # On the second round just take the majority of annotations if there is still an ambiguous cluster type
                annot = adata_comb.obs.loc[ref_cells, ref_cell_type_column].mode()[0]
                if added_col in adata_comb.obs.columns:
                    if isinstance(adata_comb.obs[added_col].dtype, pd.CategoricalDtype):
                        if annot not in adata_comb.obs[added_col].cat.categories:
                            adata_comb.obs[added_col] = adata_comb.obs[added_col].cat.add_categories([annot])
                adata_comb.obs.loc[qry_cells, added_col] = annot
                # If there is still ambiguity on the second round just call it unknown
                # adata_comb.obs.loc[qry_cells, added_col] = "unknown"
        else: 
            annot = ref_annot[0]
            if added_col in adata_comb.obs.columns:
                if isinstance(adata_comb.obs[added_col].dtype, pd.CategoricalDtype):
                    if annot not in adata_comb.obs[added_col].cat.categories:
                        adata_comb.obs[added_col] = adata_comb.obs[added_col].cat.add_categories([annot])
            adata_comb.obs.loc[qry_cells, added_col] = annot

    adata_comb.uns['c2c_allcools_integration_results'] = {
        "ref_group" : pd.DataFrame(ref_group , columns=['ref_group']),
        "qry_group": pd.DataFrame(query_group, columns=['qry_group']),
        "confusion_matrix": confusion_matrix,
        "integration_group_cluster": integration_group_cluster
    }
    return adata_comb

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
    label_transfer_k:int = 100,
    save_label_transfer_path: str | Path | None = None,
    save_query: bool = True,
    save_integrator:bool = True, 
    save_adata_comb: bool = True,
    run_joint_embeddings: bool = False,
    joint_embedding_leiden_res: float = 2.5,
    run_harmony_joint_embeddings: bool = False,
    harmony_nclust_joint_embeddings: int = 30,
    max_iter_harmony_joint_embeddings: int = 20,
    harmony_batch_keys: str | list[str] = "Modality",
    filter_on_score: bool = False,
    transfer_score_threshold: float = 0.5,
    run_clust_label_transfer: bool = True,
    c2c_method: str = "proportion",
    ref_proportion: float = 0.8,
    ref_cell_count: int = 5,
    qry_cell_count: int = 5,
    confusion_matrix_cluster_min_value = 0.25,
    confusion_matrix_cluster_max_value = 0.9,
    confusion_matrix_cluster_resolution = 1.5,
    qry_only_cluster_threshold = 50,
    ref_only_cluster_threshold = 20,
    integration_round: int = 1,
    normalize_ref: bool = False,
    normalize_qry: bool = False,
    ref_raw_key: str = "counts",
    qry_raw_key: str = "counts",
    deg_type: str = "deg", # "cef" or "deg" - whether to use cluster enriched features (cef) or differentially expressed genes (deg) for integration
    cef_column : str = None,  # if using cef, the column in ref_adata.var to use for cefs
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

    if top_deg_genes > 0:
        if deg_type == "cef":
            logger.info("Using Cluster Enriched Features (CEF) for integration")
            cef_col = cef_column if cef_column is not None else rna_cell_type_column
            deg_genes = _allcools_cef(
                ref_adata,
                qry_adata_ds,
                ref_col_level=cef_col,
                qry_col_level=qry_cluster_column,
                top_n=top_deg_genes,
                alpha=0.05,
            )
        else: 
            logger.info("Using Differentially Expressed Genes (DEG) for integration")
            deg_genes = _get_joined_deg_list(
                ref_adata,
                rna_cell_type_column,
                qry_adata_ds,
                qry_cluster_column,
                top_n=top_deg_genes
            )
        logger.info(f"Using {len(deg_genes)} DEGs for integration")
        ref_adata = ref_adata[:, ref_adata.var_names.isin(deg_genes)].copy()
        qry_adata_ds = qry_adata_ds[:, qry_adata_ds.var_names.isin(deg_genes)].copy()

    # Normalize the data if needed
    if normalize_ref:
        logger.info("Normalizing reference data")
        normalize_adata(ref_adata, layer=ref_raw_key)
    if normalize_qry:
        logger.info("Normalizing query data")
        normalize_adata(qry_adata_ds, layer=qry_raw_key)

    # Run PCA on the reference and query data
    ref_adata, qry_adata_ds = _allcools_pca(
        ref_adata, qry_adata_ds, n_train_cell=n_train_cell, chunk_size=chunk_size
    )

    # Get the Seurat Integration object
    adata_comb, integrator = _allcools_seurat_wrapper(ref_adata, qry_adata_ds, **kwargs)
    pre_cols = set(qry_adata_ds.obs.columns).union(ref_adata.obs.columns)
    del ref_adata, qry_adata_ds


    # Transfer labels to adata:
    adata_comb, _score_df = _transfer_labels(
        adata_comb,
        integrator,
        rna_cell_type=rna_cell_type_column,
        label_transfer_k=label_transfer_k,
        save_transfer=save_label_transfer_path
    )

    # Filter by transfer score if requested
    if filter_on_score:
        filt_adata_allcools(
            adata_comb,
            label_col=f"allcools_{rna_cell_type_column}",
            transfer_score_threshold=transfer_score_threshold,
        )

    # Run joint embeddings if needed
    if run_joint_embeddings:
        logger.info("Calculating joint embeddings")
        _joint_embeddings(
            adata_comb,
            use_rep="X_pca_integrate",
            key_added="integrated_",
            leiden_res=joint_embedding_leiden_res,
            run_harmony=run_harmony_joint_embeddings,
            batch_key=harmony_batch_keys,
            harmony_nclust=harmony_nclust_joint_embeddings,
            max_iter_harmony=max_iter_harmony_joint_embeddings,
            **kwargs)
        
    # Running the cluster2cluster label transfer
    if run_clust_label_transfer:
        if c2c_method == "allcools":
            adata_comb = _allcools_clust2clust_transfer(
                adata_comb,
                ref_cell_type_column=rna_cell_type_column,
                joint_cluster_column="integrated_leiden",
                confusion_matrix_cluster_min_value=confusion_matrix_cluster_min_value,
                confusion_matrix_cluster_max_value=confusion_matrix_cluster_max_value,
                confusion_matrix_cluster_resolution=confusion_matrix_cluster_resolution,
                qry_only_cluster_threshold=qry_only_cluster_threshold,
                ref_only_cluster_threshold=ref_only_cluster_threshold,
                integration_round=integration_round,
            )
        else:
            adata_comb = _clust2clust_transfer(
                adata_comb,
                ref_cell_type_column=rna_cell_type_column,
                joint_cluster_column="integrated_leiden",
                ref_proportion=ref_proportion,
                ref_cell_count=ref_cell_count,
                qry_cell_count=qry_cell_count,
            )


    post_cols = adata_comb.obs.columns
    cols_to_add = list(set(post_cols) - set(pre_cols))
    qry_adata = _add_obs_columns(adata_comb, qry_adata, cols_to_add)

    # Store full per-label score matrix so downstream QC can inspect non-max scores.
    # obsm key: allcools_{col}_score_matrix  (n_cells × n_labels numpy array)
    # uns  key: allcools_{col}_score_labels  (list of label names)
    _score_key  = f"allcools_{rna_cell_type_column}_score_matrix"
    _labels_key = f"allcools_{rna_cell_type_column}_score_labels"
    qry_adata.obsm[_score_key] = _score_df.loc[qry_adata.obs_names].values
    qry_adata.uns[_labels_key] = _score_df.columns.tolist()

    if save_query:
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

    return qry_adata, adata_comb



def _add_obs_columns(adata_source, adata_target, cols_to_add): 
    """
    Add new columns from adata_comb to qry_adata.obs, ensuring categorical columns have matching categories.
    """
    logger.info(f"Adding columns to qry_adata.obs: {cols_to_add}")
    # Ensure categorical columns have matching categories before assignment
    for col in cols_to_add:
        if col in adata_source.obs.columns:
            source_series = adata_source.obs[col].copy()
            if source_series.dtype.name == 'category':
                # If the column already exists in qry_adata.obs and is categorical
                if col in adata_target.obs.columns and adata_target.obs[col].dtype.name == 'category':
                    # Union the categories from both series
                    combined_categories = source_series.cat.categories.union(adata_target.obs[col].cat.categories)
                    source_series = source_series.cat.add_categories(
                        combined_categories.difference(source_series.cat.categories)
                    )
                    adata_target.obs[col] = adata_target.obs[col].cat.add_categories(
                        combined_categories.difference(adata_target.obs[col].cat.categories)
                    )
            adata_target.obs[col] = source_series
        
    # For any remaining columns, use the original assignment
    remaining_cols = [col for col in cols_to_add if col not in adata_target.obs.columns]
    if remaining_cols:
        adata_target.obs[remaining_cols] = adata_source.obs[remaining_cols]
    return adata_target