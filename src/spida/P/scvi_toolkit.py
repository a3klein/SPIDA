import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi

from ALLCools.clustering import tsne  # type: ignore


# A function that moves the manifold coordinate
def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f"{name}_{i}"] = adata.obsm[f"X_{name}"][:, i]
    return adata


# For doublet detection and removal using SOLO
def identify_doublets(adata: ad.AnnData, threshold: float = 0.5, **model_kwargs):
    # For filtering out the already qc'd cells (so assignment does not fail downstream)
    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_bool"] = False
    # Making a copy of adata to avoid modifying the original adata
    adata_train = adata.copy()
    if "pass_qc" in adata_train.obs.columns:
        adata_train = adata_train[adata_train.obs["pass_qc"], :].copy()

    if "raw" not in adata_train.layers:
        adata_train.layers["raw"] = adata_train.X.copy()

    scvi.model.SCVI.setup_anndata(adata_train, layer="raw")
    vae = scvi.model.SCVI(adata_train)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    predictions = solo.predict()
    adata.obs["doublet_score"] = predictions["doublet"]
    adata.obs["doublet_bool"] = adata.obs["doublet_score"] > threshold

    del adata_train
    return adata


def remove_doublets(adata):
    adata = adata[~adata.obs["doublet_bool"], :].copy()
    adata.obs.drop(columns=["doublet_score", "doublet_bool"], inplace=True)
    return adata


def resolvi_cluster(
    adata,
    max_epochs: int = 200,
    layer : str = "raw",
    batch_key: str = 'dataset_id',
    categorical_covariates: list = None,
    sample_quantile : str = 'post_sample_q50',
    model_save_path : str = None, 
    model_save_ext : str = None,
    trained : bool = False,
    **model_kwargs,
):
    """ Running RESOLVI clustering on the given AnnData object.
    Parameters
    ----------
    adata : anndata.AnnData
        The annotated data matrix.
    max_epochs : int, optional
        The maximum number of epochs to train the model, by default 200.
    layer : str, optional
        The layer of the AnnData object to use for training, by default "raw".
    batch_key : str, optional
        The key in `adata.obs` that contains batch information, by default None.
    categorical_covariates : list, optional
        List of categorical covariate keys in `adata.obs`, by default None.
    sample_quantile : str, optional
        The quantile to use for sampling posterior predictive, by default 'post_sample_q50'.
    model_save_path : str, optional
        Path to save the trained model, by default None.
    model_kwargs : dict, optional
        Additional keyword arguments for the RESOLVI model.
    """
    # import scipy.sparse as scp

    import logging
    logger = logging.getLogger(__package__)

    logger.info("Arguments received for RESOLVI clustering:")
    logger.info(f"max_epochs: {max_epochs}")
    logger.info(f"layer: {layer}")
    logger.info(f"batch_key: {batch_key}")
    logger.info(f"categorical_covariates: {categorical_covariates}")
    logger.info(f"model_save_path: {model_save_path}")
    logger.info(f"model_save_ext: {model_save_ext}")
    logger.info(f"trained: {trained}")
    

    # SCVI stores spatial coordinates in 'X_spatial'
    adata.obsm["X_spatial"] = adata.obsm["spatial"]

    # Making a copy of adata to avoid modifying the original
    adata_train = adata.copy()
    
    if trained and (model_save_path is not None): 
        resolvi = scvi.external.RESOLVI.load(
            model_save_path / f'{model_save_ext}_resolvi/',
            adata=adata_train,
        )
    else: 
        # Setup the model
        scvi.external.RESOLVI.setup_anndata(
            adata_train,
            layer=layer,
            batch_key=batch_key,
            categorical_covariate_keys=categorical_covariates,
        )

        resolvi = scvi.external.RESOLVI(adata_train, downsample_counts=False)
        resolvi.train(max_epochs=max_epochs)
        # resolvae.train(max_epochs=100, lr=1e-3, lr_extra=5e-2, enable_progress_bar=False)
        if model_save_path is not None:
            temp_out_path = model_save_path / f'{model_save_ext}_resolvi/'
            logger.info(f"Saving RESOLVI model to {temp_out_path}")
            resolvi.save(temp_out_path, save_anndata=False, overwrite=True)

    adata_train.obsm["X_resolvi"] = resolvi.get_latent_representation(adata_train)

    sample_mixtures = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=["mixture_proportions"],
        summary_fun={"post_sample_means": np.mean},
        num_samples=30, return_samples=False, batch_size=1000
    )
    sample_mixtures = pd.DataFrame(sample_mixtures).T
    adata_train.obs[
        ["true_proportion", "diffusion_proportion", "background_proportion"]
    ] = sample_mixtures.loc["post_sample_means", "mixture_proportions"]
    
    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={sample_quantile: np.median},
        summary_frequency=30,
        num_samples=30,
        return_samples=False,
        batch_size=1000,        
    )
    samples_corr = pd.DataFrame(samples_corr).T
    adata_train.layers["generated_expression"] = samples_corr.loc[sample_quantile, "px_rate"]

    ### move from adata_train to adata
    adata.obsm["X_resolvi"] = adata_train.obsm["X_resolvi"]
    adata.obs["leiden_resolvi"] = adata_train.obs["leiden_resolvi"]
    adata.obs["true_proportion"] = adata_train.obs["true_proportion"]
    adata.obs["diffusion_proportion"] = adata_train.obs["diffusion_proportion"]
    adata.obs["background_proportion"] = adata_train.obs["background_proportion"]
    adata.layers["generated_expression"] = adata_train.layers["generated_expression"]

    return adata

# CODE FROM THE RESOLVI PAPER 
# samples_corr = resolvae.sample_posterior_predictive(
#     model=resolvae.module.model_corrected,
#     return_sites=['px_rate', 'obs'],
#     num_samples=30, return_samples=False, batch_size=1000, macro_batch_size=50000)
# samples_corr = pd.DataFrame(samples_corr).T

# samples = resolvae.sample_posterior_predictive(
#     model=resolvae.module.model_residuals,
#     return_sites=[
#         'mixture_proportions', 'mean_poisson', 'per_gene_background', 
#         'diffusion_mixture_proportion', 'per_neighbor_diffusion', 'px_r_inv'
#         ],
#     num_samples=30, return_samples=False, batch_size=1000, macro_batch_size=50000)
# samples = pd.DataFrame(samples).T

# adata.obs['true_proportion'] = samples.loc['post_sample_means', 'mixture_proportions'][:, 0]
# adata.obs['diffusion_proportion'] = samples.loc['post_sample_means', 'mixture_proportions'][:, 1]
# adata.obs['background_proportion'] = samples.loc['post_sample_means', 'mixture_proportions'][:, 2]
# adata.varm['background'] = samples.loc['post_sample_means', 'per_gene_background'].squeeze().T
# adata.var['px_r'] = 1/(1e-6 + samples.loc['post_sample_means', 'px_r_inv'][0, :])

# _ = plt.hist(adata.obs['true_proportion'], bins=30, range=(0,1))
# _ = plt.hist(adata.obs['diffusion_proportion'], bins=30, range=(0,1))
# _ = plt.hist(adata.obs['background_proportion'], bins=30, range=(0,1))
# plt.legend(['True_proportions', 'Diffusion_Proportion', 'Background_Proportion'])
# plt.savefig(f'{output_dir}/histogram_proportions.pdf')

# adata.obsm["X_resolVI"] = resolvae.get_latent_representation()
# adata.layers["generated_expression"] = scp.csr_matrix(samples_corr.loc[sample_quantile, 'obs'])
# adata.layers["corrected_counts"] = adata.layers['counts'].multiply((samples_corr.loc[sample_quantile, 'px_rate'] / (
#     1.0 + samples_corr.loc[sample_quantile, 'px_rate'] + samples.loc['post_sample_means', 'mean_poisson']))).tocsr()


# sc.pp.neighbors(adata_train, use_rep="X_resolvi")
    # sc.tl.umap(adata_train)
    # dump_embedding(adata_train, "umap")

    # tsne(
    #     adata_train,
    #     obsm="X_resolvi",
    #     metric="euclidean",
    #     exaggeration=-1,
    #     perplexity=50,
    #     n_jobs=-1,
    # )
    # dump_embedding(adata_train, "tsne")

    # # Clustering using leiden
    # sc.tl.leiden(
    #     adata_train,
    #     resolution=0.5,
    #     random_state=13,
    #     flavor="igraph",
    #     n_iterations=2,
    #     key_added="leiden_resolvi",
    # )


    # adata.obs["umap_0"] = adata_train.obs["umap_0"]
    # adata.obs["umap_1"] = adata_train.obs["umap_1"]
    # adata.obs["tsne_0"] = adata_train.obs["tsne_0"]
    # adata.obs["tsne_1"] = adata_train.obs["tsne_1"]
    # adata.obsm["X_umap"] = adata_train.obsm["X_umap"]
    # adata.obsm["X_tsne"] = adata_train.obsm["X_tsne"]