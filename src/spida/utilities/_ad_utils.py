import numpy as np
import pandas as pd
import anndata as ad

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
