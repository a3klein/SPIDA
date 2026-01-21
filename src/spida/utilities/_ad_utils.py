import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

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

#TODO: @a3klein move _calc_embeddings and multi_round_clustering to _ad_utils. Perhaps look to ALLCools consensus clustering for ideas on how to improve clustering 
# Those functions are backend, I should expose them as importable functions from spida.P in general.