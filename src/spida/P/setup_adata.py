import scanpy as sc
import anndata as ad

from ALLCools.clustering import tsne, significant_pc_test  # type: ignore


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

    # standard scanpy analysis pipeline
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)

    # manifold projections:
    sc.tl.umap(adata)
    dump_embedding(adata, "umap")

    tsne(
        adata,
        obsm="X_pca",
        metric="euclidean",
        exaggeration=-1,
        perplexity=50,
        n_jobs=-1,
    )
    dump_embedding(adata, "tsne")

    # Clustering using leiden
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    return adata

def _calc_embeddings(
    adata : ad.AnnData, 
    layer : str | None = None,
    use_rep : str | None = None, 
    key_added : str = "X_",
    leiden_res : int = 0.7,
    min_dist=0.25,
    knn:int=25,
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
        significant_pc_test(adata, p_cutoff=0.1, update=True, obsm=f"{key_added}pca", downsample=100000)
        use_rep = f"{key_added}pca"

    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=knn)
    sc.tl.umap(adata, min_dist=min_dist, random_state=13, key_added=f"X_{key_added}umap")
    dump_embedding(adata, from_name=f"{key_added}umap", to_name=f"{key_added}umap")

    temp_tsne = adata.obsm['X_tsne'].copy()
    tsne(
        adata,
        obsm=use_rep,
        metric="euclidean",
        exaggeration=-1,
        perplexity=50,
        n_jobs=-1,
    )
    adata.obsm[f"X_{key_added}tsne"] = adata.obsm['X_tsne'].copy()
    adata.obsm['X_tsne'] = temp_tsne  # restore original tsne
    dump_embedding(adata, from_name=f"{key_added}tsne", to_name=f"{key_added}tsne")
    
    # Clustering using leiden
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, key_added=f"{key_added}leiden",
                 resolution=leiden_res)

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

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["normalized"] = adata.X.copy()

    if scale: 
        sc.pp.scale(adata)
        adata.layers["scaled"] = adata.X.copy()

    sc.pp.pca(adata, n_comps=50, chunked=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, random_state=0, min_dist=0.25, spread=1)
    sc.tl.leiden(adata, resolution=0.8, random_state=13, flavor="igraph", n_iterations=2)

    # manifold projections:
    dump_embedding(adata, "umap", "base_umap")
    adata.obsm["X_base_umap"] = adata.obsm["X_umap"].copy()

    tsne(
        adata,
        obsm="X_pca",
        metric="euclidean",
        exaggeration=-1,
        perplexity=50,
        n_jobs=-1,
    )
    dump_embedding(adata, "tsne", "base_tsne")
    adata.obsm["X_base_tsne"] = adata.obsm["X_tsne"].copy()

    return adata


# combine adatas (into a single object across donors?)
# --> This will not be associated with a single spatialdata object?
# _scvi_merge adatas (?)
# --> Need to make one spatialdata object per experiment as opposed to splitting them up by region!
# This is for segmentation testing, so each region is individually done, once there is
# a consensus segmentation the merged spatialdata objects would make more sense.
