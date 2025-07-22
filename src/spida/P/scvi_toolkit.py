import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi


# For doublet detection and removal using SOLO
def identify_doublets(adata, threshold=0.5): 
    scvi.model.SCVI.setup_anndata(adata, layer='raw')
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    predictions = solo.predict()
    adata.obs['doublet_score'] = predictions['doublet']
    adata.obs['doublet_bool'] = adata.obs['doublet_score'] > threshold
    return adata

def remove_doublets(adata): 
    adata = adata[~adata.obs['doublet_bool'], :].copy()
    adata.obs.drop(columns=['doublet_score', 'doublet_bool'], inplace=True)
    return adata

def resolvi_cluster(adata): 
    # SCVI stores spatial coordinates in 'X_spatial'
    adata.obsm['X_spatial'] = adata.obsm['spatial']
    
    # Making a copy of adata to avoid modifying the original
    adata_train = adata.copy()
    scvi.external.RESOLVI.setup_anndata(adata_train, layer='raw', batch_key='region')
    
    resolvi = scvi.external.RESOLVI(adata_train, downsample_counts=False) 
    resolvi.train(max_epochs=100)
    adata_train.obsm['X_resolvi'] = resolvi.get_latent_representation(adata_train)
    sc.pp.neighbors(adata_train, use_rep='X_resolvi')
    sc.tl.umap(adata_train)
    sc.tl.leiden(adata, resolution=0.5, random_state=13, flavor="igraph", n_iterations=2)

    sample_mixtures = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=[
            'mixture_proportions'
            ],
        summary_fun={"post_sample_means": np.mean},
        num_samples=3)
    sample_mixtures = pd.DataFrame(sample_mixtures).T
    # adata_train.obs[['true_proportion', 'diffusion_proportion', 'background_proportion']] = sample_mixtures.loc['post_sample_means', 'mixture_proportions']

    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={"post_sample_q50": np.median},
        num_samples=3,
        summary_frequency=30,
    )
    samples_corr = pd.DataFrame(samples_corr).T

    adata.layers["generated_expression"] = samples_corr.loc["post_sample_q50", "px_rate"]


