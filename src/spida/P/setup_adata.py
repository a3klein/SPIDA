import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scipy.sparse as scp
import numpy as np
import pandas as pd

from ALLCools.clustering import tsne, significant_pc_test  # type: ignore

from ..utilities.ad_utils import (
    normalize_adata, dump_embedding, _calc_embeddings, multi_round_clustering
)

# Individual adatas
def run_setup(
    adata: ad.AnnData,
    exp_name: str,
    reg_name: str,
    seg_name: str,
    donor_name: str = None,
):
    """
    Setup the AnnData object for further analysis.

    Parameters:
    adata (AnnData): The AnnData object to be set up.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    seg_name (str): Name of the segmentation.
    donor_name (str, optional): Name of the donor. Defaults to None.

    Returns:
    AnnData: The modified AnnData object.
    """

    # filtering for QC'd cells + doublet cells
    adata = adata[adata.obs["pass_qc"]].copy()
    if "doublet_bool" in adata.obs.columns:
        adata = adata[~adata.obs["doublet_bool"]].copy()
    # Set metadata
    adata.uns["experiment"] = exp_name
    adata.uns["region"] = reg_name
    adata.uns["segmentation"] = seg_name
    adata.uns["donor"] = donor_name

    # backup raw
    adata.layers["raw"] = adata.X.copy()
    # create the spatial coordinates in obsm
    adata.obsm["spatial"] = adata.obs[["CENTER_X", "CENTER_Y"]].values

    # Normalize gene count by volume
    adata.layers['volume_norm'] = scp.csr_matrix(adata.X / adata.obs['volume'].values[:, np.newaxis])
    adata.X = adata.layers['volume_norm'].copy()
    normalize_adata(adata, layer='volume_norm', log1p=True)
    _calc_embeddings(adata, layer=None, key_added="base_", leiden_res=1, knn=35, consensus_cluster=True)

    try:
        from ..utilities.degs import call_degs_scran
        call_degs_scran(adata, celltype_col="base_leiden", uns_key="scran_markers", num_threads=4)
    except Exception:
        pass
    sc.tl.rank_genes_groups(adata, groupby="base_leiden", method="t-test_overestim_var")

    return adata


def combined_setup(
    adata : ad.AnnData, 
    scale: bool = False
): 
    """
    Setup the AnnData object for the combined datasets. 
    """
    import scipy.sparse as sp
    
    if "raw" in adata.layers:
        adata.X = adata.layers['raw'].copy()
        del adata.layers['raw']  # remove raw layer to save memory

        adata.layers["counts"] = sp.csr_matrix(adata.X.copy())
        adata.X = adata.layers['counts'].copy()
    else: 
        adata.X = adata.layers['counts'].copy()

    # TODO: Fix this stuff (right now this is just totally wrong so commenting out)
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    # adata.layers["normalized"] = adata.X.copy()

    # if scale: 
    #     sc.pp.scale(adata)
    #     adata.layers["scaled"] = adata.X.copy()

    # sc.pp.pca(adata, n_comps=50, chunked=True)
    # sc.pp.neighbors(adata)
    # sc.tl.umap(adata, random_state=0, min_dist=0.25, spread=1)
    # sc.tl.leiden(adata, resolution=0.8, random_state=13, flavor="igraph", n_iterations=2)

    # # manifold projections:
    # dump_embedding(adata, "umap", "base_umap")
    # adata.obsm["X_base_umap"] = adata.obsm["X_umap"].copy()

    # tsne(
    #     adata,
    #     obsm="X_pca",
    #     metric="euclidean",
    #     exaggeration=-1,
    #     perplexity=50,
    #     n_jobs=-1,
    # )
    # dump_embedding(adata, "tsne", "base_tsne")
    # adata.obsm["X_base_tsne"] = adata.obsm["X_tsne"].copy()

    return adata


# combine adatas (into a single object across donors?)
# --> This will not be associated with a single spatialdata object?
# _scvi_merge adatas (?)
# --> Need to make one spatialdata object per experiment as opposed to splitting them up by region!
# This is for segmentation testing, so each region is individually done, once there is
# a consensus segmentation the merged spatialdata objects would make more sense.
