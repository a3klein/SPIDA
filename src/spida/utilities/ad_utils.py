import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp


def dump_embedding(
    adata : ad.AnnData,
    from_name : str,
    to_name : str = None,
    n_dim : int = 2
):
    """
    Move coordinates from adata.obsm to adata.obs
    """
    if to_name is None:
        to_name = from_name
    for i in range(n_dim):
        adata.obs[f"{to_name}_{i}"] = adata.obsm[f"X_{from_name}"][:, i]
    return adata

def _downsample_ref_clusters(
    adata : ad.AnnData, 
    col : str,
    max_cells = 3000,
    random_state = 13
): 
    """Downsample reference AnnData object by clusters to a maximum number of cells per cluster."""
    if col not in adata.obs.columns: 
        raise ValueError(f"Column {col} not found in adata.obs")
    vc = adata.obs[col].value_counts()
    keep_cells = set()
    for cat, n in vc.items(): 
        if n > max_cells: 
            sampled = adata.obs_names[adata.obs[col] == cat].to_series().sample(max_cells, random_state=random_state).tolist()
            keep_cells.update(sampled)
        else: 
            sampled = adata.obs_names[adata.obs[col] == cat].to_list()
            keep_cells.update(sampled)
    return adata[adata.obs.index.isin(keep_cells)].copy()


def _remove_small_clusters(
    adata : ad.AnnData, 
    col : str, 
    min_cells : int = 20
): 
    """Remove clusters from the AnnData object that have fewer than min_cells."""
    if col not in adata.obs.columns: 
        raise ValueError(f"Column {col} not found in adata.obs")

    vc = adata.obs[col].value_counts()
    keep_cells = vc[vc >= min_cells].index
    return adata[adata.obs[col].isin(keep_cells)].copy()


# Apply normalization to the .X in either AnnData Object
def normalize_adata(
    adata: ad.AnnData,
    layer: str = None,
    target_sum: float = None,
    log1p: bool = True,
): 
    if layer is not None: 
        adata.X = adata.layers[layer].copy() # X needs to be a csr_matrix 
    if not isinstance(adata.X, sp.csr_matrix): 
        adata.X = sp.csr_matrix(adata.X)
    n_counts = np.ravel(adata.X.sum(axis=1))
    target = np.median(n_counts) if target_sum is None else target_sum
    adata.X.data = adata.X.data/np.repeat(n_counts, adata.X.getnnz(axis=1)) * target
    if log1p: 
        sc.pp.log1p(adata)

#TODO: @a3klein Look to ALLCools consensus clustering for ideas on how to improve clustering 
# Those functions are backend, I should expose them as importable functions from spida.P in general.
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
    from ALLCools.clustering import tsne, significant_pc_test  # type: ignore
    
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
    # TODO: @a3klein: Replace leiden clustering with ALLCools Consensus Clustering (or at least add that as an option) for more robust clustering results.

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

