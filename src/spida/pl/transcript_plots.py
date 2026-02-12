import numpy as np
import matplotlib.pyplot as plt
from spida.pl import plot_categorical, plot_continuous

def plot_hex_qc(grid, cmap="RdYlGn_r"):
    grid_filt = grid[~grid["filtered"]] if "filtered" in grid.columns else grid
    grid_filt['log_tz_count'] = np.log1p(grid_filt['tz_count'])

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=150)
    ax = axes[0]
    grid_filt.plot(ax=ax, column="tz_count", cmap=cmap, legend=True).axis("off")
    ax.set_title("total transcripts")

    ax = axes[1]
    grid_filt.plot(ax=ax, column="log_tz_count", cmap=cmap, legend=True).axis("off")
    ax.set_title("log(total transcripts+1)")

    ax = axes[2]
    grid_filt.plot(ax=ax, column="n_genes", cmap=cmap, legend=True).axis("off")
    ax.set_title("n_genes")

    plt.tight_layout()
    return fig

def plot_hex_clusters(grid, adata_hex, cmap="RdYlGn_r"):
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=150)

    if "X_umap" in adata_hex.obsm:
        ax = axes[0]
        plot_categorical(adata_hex, cluster_col="leiden", coord_base="umap", show=False, text_anno=True, coding=True, ax=ax)
        ax.set_title("UMAP - leiden")

        ax = axes[1]
        plot_continuous(adata_hex, color_by="tz_count", coord_base="umap", show=False, ax=ax)
        ax.set_title("UMAP - tz_count")

    grid_filt = grid[~grid["filtered"]] if "filtered" in grid.columns else grid
    if "leiden" in adata_hex.obs.columns:
        grid_filt = grid_filt.copy()
        grid_filt["leiden"] = adata_hex.obs["leiden"]

    fig2, axes2 = plt.subplots(1, 3, figsize=(12, 4), dpi=150)
    ax = axes2[0]
    grid_filt.plot(ax=ax, column="tz_count", cmap=cmap, legend=True).axis("off")
    ax.set_title("total transcripts")

    ax = axes2[1]
    if "leiden" in grid_filt.columns:
        grid_filt.plot(ax=ax, column="leiden", categorical=True, legend=True).axis("off")
        ax.set_title("leiden")
    else:
        ax.axis("off")

    ax = axes2[2]
    grid_filt.plot(ax=ax, column="n_genes", cmap=cmap, legend=True).axis("off")
    ax.set_title("n_genes")

    plt.tight_layout()
    return fig, fig2