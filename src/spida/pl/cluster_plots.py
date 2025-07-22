import matplotlib.pyplot as plt
import seaborn as sns


def plot_scatter(adata, title="", colors=['volume', 'nCount_RNA'], ncols=2, 
              coord_base="umap", pdf_file=None):
    
    # Set rcparam
    plt.rcParams['axes.facecolor'] = 'white'

    ds = min(200000/adata.shape[0], 4) # the point size to use for plotting
    nrows = 1

    if coord_base in adata.obsm:
        arr = adata.obsm[coord_base]
    elif 'X_' + coord_base in adata.obsm:
        arr = adata.obsm['X_' + coord_base]
    else: 
        raise ValueError(f"Coordinate base '{coord_base}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")

    fig, axes = plt.subplots(nrows, ncols, figsize = (ncols*5, nrows*3), # gridspec_kw={"width_ratios": [4, 1, 4, 1, 4, 1]},
                             dpi=300, constrained_layout=True,)


    for c, ax in zip(colors, axes):
        if adata.obs[c].dtype.name == 'category':
            sns.scatterplot(x=arr[:,0], y=arr[:,1], s=ds, hue=adata.obs[c], marker='.', ax=ax,
                            palette='tab20', linewidth=0, alpha=1, legend=False, rasterized=True)
        else: 
            sns.scatterplot(x=arr[:,0], y=arr[:,1], s=ds, hue=adata.obs[c], marker='.', ax=ax,
                            palette='viridis', linewidth=0, alpha=1, legend=False, rasterized=True)
    plt.suptitle(title)
    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    else: 
        plt.show()
        plt.close()