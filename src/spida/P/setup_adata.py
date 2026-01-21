import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scipy.sparse as scp
import numpy as np
import pandas as pd

from ALLCools.clustering import tsne, significant_pc_test  # type: ignore

from ..utilities._ad_utils import normalize_adata


# A function that moves the manifold coordinate
def dump_embedding(adata, from_name, to_name=None, n_dim=2):
    if to_name is None:
        to_name = from_name
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f"{to_name}_{i}"] = adata.obsm[f"X_{from_name}"][:, i]
    return adata


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
    _calc_embeddings(adata, layer=None, key_added="base_", leiden_res=1, knn=35)
    
    # sc.pp.normalize_total(adata)
    # sc.tl.pca(adata)
    # sc.pp.neighbors(adata)

    # # manifold projections:
    # sc.tl.umap(adata)
    # dump_embedding(adata, "umap")

    # tsne(
    #     adata,
    #     obsm="X_pca",
    #     metric="euclidean",
    #     exaggeration=-1,
    #     perplexity=50,
    #     n_jobs=-1,
    # )
    # dump_embedding(adata, "tsne")

    # # Clustering using leiden
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    return adata

from ALLCools.clustering import tsne, significant_pc_test  # type: ignore
from spida.P.setup_adata import dump_embedding
def _calc_embeddings(
    adata : ad.AnnData, 
    layer : str | None = None,
    use_rep : str | None = None, 
    key_added : str = "X_",
    leiden_res : int = 0.7,
    min_dist=0.25,
    knn:int=25,
    p_cutoff: float = 0.1,
    run_embedding: bool = True,
    run_harmony: bool = False,
    batch_key: str | list[str] = "dataset_id",
    harmony_nclust: int = 20,
    max_iter_harmony: int = 20,
): 
    """ Calculate the UMAP and tSNE embeddings for the AnnData object and given layer"""
    if (use_rep is not None) and (layer is not None): 
        raise ValueError("Either use_rep or layer can be specified, not both.")
    
    chunked = (adata.shape[0] > 50000) and (layer is None)
    
    if use_rep is None: 
        sc.pp.pca(
            adata,
            n_comps=min(min(adata.shape) - 1, 100),
            chunked=chunked,
            layer=layer,
            key_added=f"{key_added}pca"
        )
        significant_pc_test(adata, p_cutoff=p_cutoff, update=True, obsm=f"{key_added}pca", downsample=100000)
        use_rep = f"{key_added}pca"

    
    if run_harmony:
        if isinstance(batch_key, str):
            if batch_key not in adata.obs.columns:
                raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
            batch_key = [batch_key]
        elif isinstance(batch_key, list):
            for key in batch_key:
                if key not in adata.obs.columns:
                    raise ValueError(f"batch_key '{key}' not found in adata.obs")
        else:
            raise ValueError("batch_key must be a string or a list of strings")

        for i, _key in enumerate(batch_key):
            adata.obs[_key] = adata.obs[_key].astype("category")
            basis = f"{key_added}pca" if use_rep is None else use_rep
            sce.pp.harmony_integrate(adata,
                                    key=_key,
                                    basis=basis if i == 0 else f"{key_added}pca_harmony",
                                    adjusted_basis=f"{key_added}pca_harmony",
                                    nclust=harmony_nclust,
                                    max_iter_harmony=max_iter_harmony) #X_pca_harmony
        use_rep = f"{key_added}pca_harmony"
    
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=knn)
    
    # If running embedding (usually in R1 = True, R2... = False)
    if run_embedding:
        sc.tl.umap(adata, min_dist=min_dist, random_state=13, key_added=f"X_{key_added}umap")
        dump_embedding(adata, from_name=f"{key_added}umap", to_name=f"{key_added}umap")

        
        temp_tsne = adata.obsm['X_tsne'].copy() if "X_tsne" in adata.obsm.keys() else None
        tsne(
            adata,
            obsm=use_rep,
            metric="euclidean",
            exaggeration=-1,
            perplexity=50,
            n_jobs=-1,
        )
        adata.obsm[f"X_{key_added}tsne"] = adata.obsm['X_tsne'].copy()
        if temp_tsne is not None:
            adata.obsm['X_tsne'] = temp_tsne
        dump_embedding(adata, from_name=f"{key_added}tsne", to_name=f"{key_added}tsne")
    
    # Clustering using leiden
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, key_added=f"{key_added}leiden",
                 resolution=leiden_res)

    return adata
    
def multi_round_clustering(
    adata : ad.AnnData,
    layer : str | None = None,
    use_rep : str | None = None,
    key_added : str = "X_",
    num_rounds : int = 2,
    leiden_res : int | list[float] = 1,
    min_dist=0.25,
    knn:int=25,
    p_cutoff: float = 0.1,
    min_group_size:int=50,
    run_harmony: bool = False,
    batch_key: str | list[str] = "dataset_id",
    harmony_nclust: int = 20,
    max_iter_harmony: int = 20,
):
    """ Perform multiple rounds of clustering on the AnnData object """
    if isinstance(leiden_res, (int, float)):
        leiden_res = [leiden_res] * num_rounds
    elif len(leiden_res) != num_rounds:
        raise ValueError("leiden_res must be a single value or a list of length num_rounds")
    
    adata = _calc_embeddings(
        adata,
        layer=layer,
        use_rep=use_rep,
        key_added=f"{key_added}round1_",
        leiden_res=leiden_res[0],
        min_dist=min_dist,
        knn=knn,
        p_cutoff=p_cutoff,
        run_embedding=True,
        run_harmony=run_harmony,
        batch_key=batch_key,
        harmony_nclust=harmony_nclust,
        max_iter_harmony=max_iter_harmony,
    )
    if run_harmony: 
        use_rep = f"{key_added}round1_pca_harmony"
        layer = None

    # TODO: Decide whether to recalculate the pca for each subgroup or use the same pca (harmony on the sub_pcas as well?)
    for i in range(1, num_rounds):
        df_to_add = pd.DataFrame(index=adata.obs.index)
        df_to_add.loc[:, f"{key_added}round{i+1}_leiden"] = "unassigned"
        for _group in adata.obs[f"{key_added}round{i}_leiden"].cat.categories:
            idx = adata.obs[f"{key_added}round{i}_leiden"] == _group
            if np.sum(idx) < min_group_size:
                df_to_add.loc[idx, f"{key_added}round{i+1}_leiden"] = _adata_sub.obs[f"{key_added}round{i}_leiden"].astype(str)
                continue
            _adata_sub = adata[idx].copy()
            _adata_sub = _calc_embeddings(
                _adata_sub,
                layer=layer,
                use_rep=use_rep,
                key_added=f"{key_added}round{i+1}_sub_",
                leiden_res=0.5,
                min_dist=min_dist,
                knn=knn,
                p_cutoff=p_cutoff,
                run_embedding=False,
                run_harmony=False,
            )
            # transfer back to original adata
            df_add = _adata_sub.obs[f"{key_added}round{i}_leiden"].astype(str) + "_" + _adata_sub.obs[f"{key_added}round{i+1}_sub_leiden"].astype(str)
            df_to_add.loc[idx, f"{key_added}round{i+1}_leiden"] = df_add
        adata.obs[f"{key_added}round{i+1}_leiden"] = df_to_add[f"{key_added}round{i+1}_leiden"].astype("category")
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
