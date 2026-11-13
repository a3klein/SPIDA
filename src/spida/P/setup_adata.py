import scanpy as sc

from spida.io.layout import RegionMeta
import anndata as ad
import scipy.sparse as scp
import numpy as np

from ..utilities.ad_utils import (
    normalize_adata, _calc_embeddings
)

# Individual adatas
def run_setup(
    adata: ad.AnnData,
    meta: RegionMeta,
    do_tsne: bool = False,
    consensus_cluster: bool = True,
    n_jobs: int | None = None,
):
    """
    Setup the AnnData object for further analysis.

    Parameters:
    adata (AnnData): The AnnData object to be set up.
    meta (RegionMeta): Experiment / region / segmentation / donor names (stored in uns).
    do_tsne (bool, optional): Whether to also compute the tSNE embedding (in addition to UMAP). Defaults to False.
    consensus_cluster (bool, optional): Use ALLCools consensus clustering (200 Leiden
        repeats across `n_jobs` workers) rather than a single Leiden run. Defaults to True.
    n_jobs (int | None, optional): Consensus clustering workers. Defaults to None: the job's cores, max 24.
    Returns:
    AnnData: The modified AnnData object.
    """

    # filtering for QC'd cells + doublet cells
    adata = adata[adata.obs["pass_qc"]].copy()
    if "doublet_bool" in adata.obs.columns:
        adata = adata[~adata.obs["doublet_bool"]].copy()
    # Set metadata
    adata.uns["experiment"] = meta.experiment
    adata.uns["region"] = meta.region
    adata.uns["segmentation"] = meta.segmentation
    adata.uns["donor"] = meta.donor

    # backup raw
    adata.layers["raw"] = adata.X.copy()
    # create the spatial coordinates in obsm
    adata.obsm["spatial"] = adata.obs[["CENTER_X", "CENTER_Y"]].values

    # Normalize gene count by volume
    adata.layers['volume_norm'] = scp.csr_matrix(adata.X / adata.obs['volume'].values[:, np.newaxis])
    adata.X = adata.layers['raw'].copy()
    normalize_adata(adata, layer='raw', log1p=True)
    _calc_embeddings(adata, layer=None, key_added="base_", leiden_res=1, knn=35,
                     consensus_cluster=consensus_cluster, do_tsne=do_tsne, n_jobs=n_jobs)

    try:
        from ..utilities.degs import call_degs_scran
        call_degs_scran(adata, celltype_col="base_leiden", uns_key="scran_markers", num_threads=4)
    except Exception:
        pass
    sc.tl.rank_genes_groups(adata, groupby="base_leiden", method="t-test_overestim_var")

    return adata


def combined_setup(
    adata: ad.AnnData,
    batch_key: str | list[str] | None = None,
    key_added: str = "",
    leiden_res: float = 1,
    knn: int = 35,
    min_dist: float = 0.25,
    do_tsne: bool = False,
    consensus_cluster: bool = True,
    n_jobs: int | None = None,
):
    """
    Joint embedding + clustering of a project's combined table.

    Parameters:
    adata (AnnData): The combined table (layers["raw"] = counts).
    batch_key (str | list[str], optional): obs column(s) to Harmony-correct on (X_<key_added>pca_harmony). None: no correction.
    key_added (str, optional): Prefix for the embedding / cluster keys. Defaults to "".
    leiden_res (float, optional): Leiden resolution. Defaults to 1.
    knn (int, optional): Neighbours for the kNN graph. Defaults to 35.
    min_dist (float, optional): UMAP min_dist. Defaults to 0.25.
    do_tsne (bool, optional): Also compute the tSNE embedding. Defaults to False.
    consensus_cluster (bool, optional): ALLCools consensus clustering rather than a single Leiden run. Defaults to True.
    n_jobs (int | None, optional): Consensus clustering workers. Defaults to None: the job's cores, max 24.
    Returns:
    AnnData: The modified AnnData object.
    """
    normalize_adata(adata, layer="raw", log1p=True)
    _calc_embeddings(adata, layer=None, key_added=key_added, leiden_res=leiden_res, knn=knn, min_dist=min_dist,
                     run_harmony=batch_key is not None, batch_key=batch_key or "dataset_id",
                     consensus_cluster=consensus_cluster, do_tsne=do_tsne, n_jobs=n_jobs)
    
    for name in (f"{key_added}pca", f"{key_added}pca_harmony"):
        if name in adata.obsm:
            adata.obsm[f"X_{name}"] = adata.obsm.pop(name)
    # TODO: NSForest markers (replacing scran DEGs) for plotting
    return adata


# combine adatas (into a single object across donors?)
# --> This will not be associated with a single spatialdata object?
# _scvi_merge adatas (?)
# --> Need to make one spatialdata object per experiment as opposed to splitting them up by region!
# This is for segmentation testing, so each region is individually done, once there is
# a consensus segmentation the merged spatialdata objects would make more sense.
