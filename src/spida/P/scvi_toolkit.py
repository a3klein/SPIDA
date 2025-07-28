import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi

from ALLCools.clustering import tsne # type: ignore


# A function that moves the manifold coordinate
def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f"{name}_{i}"] = adata.obsm[f"X_{name}"][:, i]
    return adata

# For doublet detection and removal using SOLO
def identify_doublets(
        adata:ad.AnnData,
        threshold:float=0.5, 
        **model_kwargs): 
    
    # For filtering out the already qc'd cells (so assignment does not fail downstream)
    adata.obs['doublet_score'] = np.nan
    adata.obs['doublet_bool'] = False
    # Making a copy of adata to avoid modifying the original adata
    adata_train = adata.copy()
    if "pass_qc" in adata_train.obs.columns:
        adata_train = adata_train[adata_train.obs['pass_qc'], :].copy()

    if 'raw' not in adata_train.layers:
        adata_train.layers['raw'] = adata_train.X.copy()

    scvi.model.SCVI.setup_anndata(adata_train, layer='raw')
    vae = scvi.model.SCVI(adata_train)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    predictions = solo.predict()
    adata.obs['doublet_score'] = predictions['doublet']
    adata.obs['doublet_bool'] = adata.obs['doublet_score'] > threshold
    
    del adata_train
    return adata

def remove_doublets(adata): 
    adata = adata[~adata.obs['doublet_bool'], :].copy()
    adata.obs.drop(columns=['doublet_score', 'doublet_bool'], inplace=True)
    return adata

def resolvi_cluster(adata,
                    max_epochs:int=200,
                    **model_kwargs,
                    ): 
    
    # SCVI stores spatial coordinates in 'X_spatial'
    adata.obsm['X_spatial'] = adata.obsm['spatial']
    
    # Making a copy of adata to avoid modifying the original
    adata_train = adata.copy()
    scvi.external.RESOLVI.setup_anndata(adata_train, layer='raw', batch_key='region')
    
    resolvi = scvi.external.RESOLVI(adata_train, downsample_counts=False) 
    resolvi.train(max_epochs=max_epochs)
    
    adata_train.obsm['X_resolvi'] = resolvi.get_latent_representation(adata_train)
    sc.pp.neighbors(adata_train, use_rep='X_resolvi')
    sc.tl.umap(adata_train)
    dump_embedding(adata_train, "umap")
    
    tsne(adata_train, obsm="X_resolvi", metric="euclidean",
         exaggeration=-1, perplexity=50, n_jobs=-1,)
    dump_embedding(adata_train, "tsne")

    # Clustering using leiden 
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    sc.tl.leiden(adata_train, resolution=0.5, random_state=13, flavor="igraph", n_iterations=2, key_added='leiden_resolvi')

    sample_mixtures = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=[
            'mixture_proportions'
            ],
        summary_fun={"post_sample_means": np.mean},
        num_samples=3)
    sample_mixtures = pd.DataFrame(sample_mixtures).T
    adata_train.obs[['true_proportion', 'diffusion_proportion', 'background_proportion']] = sample_mixtures.loc['post_sample_means', 'mixture_proportions']

    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={"post_sample_q50": np.median},
        num_samples=3,
        summary_frequency=30,
    )
    samples_corr = pd.DataFrame(samples_corr).T
    adata_train.layers["generated_expression"] = samples_corr.loc["post_sample_q50", "px_rate"]

    ### move from adata_train to adata
    adata.obsm['X_resolvi'] = adata_train.obsm['X_resolvi']
    adata.obsm['X_umap'] = adata_train.obsm['X_umap']
    adata.obsm['X_tsne'] = adata_train.obsm['X_tsne']
    adata.obs['leiden_resolvi'] = adata_train.obs['leiden_resolvi']
    adata.obs['true_proportion'] = adata_train.obs['true_proportion']
    adata.obs['diffusion_proportion'] = adata_train.obs['diffusion_proportion']
    adata.obs['background_proportion'] = adata_train.obs['background_proportion']
    adata.obs['umap_0'] = adata_train.obs['umap_0']
    adata.obs['umap_1'] = adata_train.obs['umap_1']
    adata.obs['tsne_0'] = adata_train.obs['tsne_0']
    adata.obs['tsne_1'] = adata_train.obs['tsne_1']
    adata.layers['generated_expression'] = adata_train.layers['generated_expression']

    return adata


