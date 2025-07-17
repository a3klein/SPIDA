import scanpy as sc
import anndata as ad

from ALLCools.clustering import tsne # type: ignore

# A function that moves the manifold coordinate
def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f"{name}_{i}"] = adata.obsm[f"X_{name}"][:, i]
    return adata


# Individual adatas
def run_setup(adata:ad.AnnData,
              exp_name:str, 
              reg_name:str,
              seg_name:str,
              donor_name:str=None,):  
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
        
    # filtering for QC'd cells
    adata = adata[adata.obs['pass_qc']].copy()
    # Set metadata
    adata.uns['experiment'] = exp_name
    adata.uns['region'] = reg_name
    adata.uns['segmentation'] = seg_name
    adata.uns['donor'] = donor_name

    # backup raw
    adata.layers['raw'] = adata.X.copy()
    # create the spatial coordinates in obsm
    adata.obsm['spatial'] = adata.obs[['CENTER_X', 'CENTER_Y']].values

    # standard scanpy analysis pipeline
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)

    # manifold projections: 
    sc.tl.umap(adata)
    dump_embedding(adata, "umap")
    
    tsne(adata, obsm="X_pca", metric="euclidean",
         exaggeration=-1, perplexity=50, n_jobs=-1,)
    dump_embedding(adata, "tsne")

    # Clustering using leiden 
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    return adata

# combine adatas (into a single object across donors?)
# --> This will not be associated with a single spatialdata object? 
# _scvi_merge adatas (?)
# --> Need to make one spatialdata object per experiment as opposed to splitting them up by region!
# This is for segmentation testing, so each region is individually done, once there is 
# a consensus segmentation the merged spatialdata objects would make more sense. 
