import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from spida.pl import plot_categorical, plot_continuous

def _get_leiden_color_dict(adata_hex):
    if "leiden" not in adata_hex.obs.columns:
        return {}
    cats = adata_hex.obs["leiden"].astype("category").cat.categories.tolist()
    colors = adata_hex.uns.get("leiden_colors", None)
    if colors is None or len(colors) < len(cats):
        colors = [plt.cm.tab20(i / max(1, len(cats) - 1)) for i in range(len(cats))]
    return {c: colors[i] for i, c in enumerate(cats)}


# TODO: in adataviz move this away from here. 

# Helper Scripts for custom dotplot based on scanpy's rank_genes_groups results.
# This avoids importing scnapy in the plotting module 
def _top_genes_from_scran_uns(adata, uns_key: str = 'scran_markers', n_genes: int = 5) -> list[str]:
    """Return a flat, deduplicated gene list from scran results stored in adata.uns."""
    from spida.utilities.degs import summarize_scran_results_from_uns
    summary = summarize_scran_results_from_uns(adata, uns_key=uns_key, top_n=n_genes)
    top_col = next((c for c in summary.columns if c.startswith('top_') and c.endswith('_markers')), None)
    if top_col is None:
        return []
    seen, genes = set(), []
    for gene_list in summary[top_col]:
        for g in (gene_list or []):
            if g not in seen:
                genes.append(g)
                seen.add(g)
    return genes


def _top_rank_genes_from_uns(adata_hex, n_genes: int = 5) -> list[str]:
    if "rank_genes_groups" not in adata_hex.uns:
        return []
    rg = adata_hex.uns["rank_genes_groups"]
    names = rg.get("names", None)
    if names is None:
        return []

    top = []
    if hasattr(names, "dtype") and names.dtype.names is not None:
        # structured array keyed by group
        for g in names.dtype.names:
            top.extend(list(names[g][:n_genes]))
    else:
        # fallback: 2D array-like [rank, group]
        arr = np.asarray(names)
        if arr.ndim == 2:
            top.extend(arr[:n_genes, :].ravel().tolist())

    # unique, preserve order
    seen = set()
    return [x for x in top if x is not None and not (x in seen or seen.add(x))]

def _plot_rank_genes_dotplot_matplotlib(
    adata_hex, n_genes: int = 5,
    cluster_col : str = "leiden",
    cmap: str = "Reds", ax=None,
    figsize=(20, 4), dpi=200
):
    """
    Plotting a dotplot of the top ranked genes from scanpy's rank_genes_groups results, using matplotlib directly.
    """

    # Getting axes if not passed in 
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure

    # Getting leiden clusters 
    if cluster_col not in adata_hex.obs.columns:
        ax.text(0.5, 0.5, f"No {cluster_col} in adata.obs", ha="center", va="center", transform=ax.transAxes)
        return fig, ax

    # Try scran markers first, fall back to scanpy rank_genes_groups
    genes = []
    if 'scran_markers' in adata_hex.uns:
        try:
            genes = _top_genes_from_scran_uns(adata_hex, uns_key='scran_markers', n_genes=n_genes)
            genes = [g for g in genes if g in adata_hex.var_names]
        except Exception:
            genes = []
    if not genes:
        genes = _top_rank_genes_from_uns(adata_hex, n_genes=n_genes)
        genes = [g for g in genes if g in adata_hex.var_names]
    
    # Making sure there are genes
    if len(genes) == 0:
        ax.text(0.5, 0.5, "No ranked genes found in adata.uns['rank_genes_groups']", ha="center", va="center", transform=ax.transAxes)
        return fig, ax

    # Getting cluster assignments per cell
    clusters = adata_hex.obs[cluster_col].astype("category").cat.categories.tolist()
    X = adata_hex[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    # Building the dotplot data: for each cluster and gene, compute % expressing and mean expression
    dot_x, dot_y, dot_s, dot_c = [], [], [], []
    leiden_vals = adata_hex.obs[cluster_col].astype(str).values
    num_clust = len(clusters)

    for yi, cl in enumerate(clusters):
        idx = np.where(leiden_vals == str(cl))[0]
        if len(idx) == 0:
            continue
        x_sub = X[idx, :]
        pct_expr = (x_sub > 0).mean(axis=0) * 100.0
        mean_expr = x_sub.mean(axis=0)

        for xi in range(len(genes)):
            dot_x.append(xi)
            dot_y.append(num_clust - yi - 1)
            dot_s.append(10 + pct_expr[xi] * 2.0)  # size ~ % expressing
            dot_c.append(mean_expr[xi])            # color ~ mean expr

    # Plotting 
    sca = ax.scatter(dot_x, dot_y, s=dot_s, c=dot_c, cmap=cmap, edgecolors="k", linewidths=0.1)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=12)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=12)
    ax.set_title("DEGs", fontsize=14)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mean expression")
    
    return fig, ax

def plot_hex_qc(grid, cmap="RdYlGn_r"):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), dpi=150)

    # panel 1: tz_count
    grid.plot(ax=axes[0], column="tz_count", cmap=cmap, legend=True)
    axes[0].set_title("total transcripts")
    axes[0].axis("off")

    # panel 2: density
    grid['log_tz_count'] = np.log1p(grid['tz_count'])
    grid.plot(ax=axes[1], column="log_tz_count", cmap=cmap, legend=True)
    axes[1].set_title("log(total transcripts+1)")
    axes[1].axis("off")

    # panel 3: n_genes
    grid.plot(ax=axes[2], column="n_genes", cmap=cmap, legend=True)
    axes[2].set_title("n_genes")
    axes[2].axis("off")

    # panel 4: QC pass/fail over all hexes
    if "filtered" in grid.columns:
        qc_val = (~grid["filtered"]).astype(int)  # 1=pass, 0=fail
        qc_cmap = ListedColormap(["#d9d9d9", "#1b9e77"])  # fail, pass
        grid.assign(qc_pass=qc_val).plot(
            ax=axes[3], column="qc_pass", cmap=qc_cmap, legend=False
        )
        handles = [
            mpatches.Patch(color="#1b9e77", label="pass"),
            mpatches.Patch(color="#d9d9d9", label="fail"),
        ]
        axes[3].legend(handles=handles, title="QC", loc="upper left", bbox_to_anchor=(1.02, 1))
        axes[3].set_title("QC pass/fail")
    else:
        axes[3].set_title("QC pass/fail (missing)")
    axes[3].axis("off")

    plt.tight_layout()
    return fig

def plot_hex_clusters(grid, adata_hex, cmap="RdYlGn_r", n_genes=5):
    
    fig1, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=200)
    if "X_umap" in adata_hex.obsm:
        ax = axes[0]
        plot_categorical(
            adata_hex, cluster_col="leiden", coord_base="umap", show=False, text_anno=True, coding=True, ax=ax
        )
        ax.set_title("UMAP - leiden")

        ax = axes[1]
        plot_continuous(adata_hex, color_by="tz_count", coord_base="umap", show=False, ax=ax)
        ax.set_title("UMAP - tz_count")

    color_dict = {}
    for _l, _c in zip(adata_hex.obs['leiden'].cat.categories, adata_hex.uns['leiden_colors']):
        color_dict[_l] = _c
    grid['leiden_color'] = grid['leiden'].map(color_dict)
    
    grid_filt = grid[~grid["filtered"]] if "filtered" in grid.columns else grid
    if "leiden" in adata_hex.obs.columns:
        grid_filt = grid_filt.copy()
        grid_filt["leiden"] = adata_hex.obs["leiden"]

    fig2, axes2 = plt.subplots(1, 3, figsize=(12, 4), dpi=200)
    ax = axes2[0]
    grid_filt.plot(ax=ax, column="tz_count", cmap=cmap, legend=True).axis("off")
    ax.set_title("total transcripts")

    ax = axes2[1]
    g = grid.copy()
    if "filtered" in g.columns:
        g[g["filtered"]].plot(ax=ax, color="#efefef", linewidth=0)
    if "leiden" in g.columns:
        for cl, col in color_dict.items():
            sub = g[g["leiden"] == cl]
            if len(sub) > 0:
                sub.plot(ax=ax, color=col, linewidth=0)
        handles = [mpatches.Patch(color=color_dict[k], label=str(k)) for k in adata_hex.obs["leiden"].astype("category").cat.categories]
        ax.legend(handles=handles, title="Leiden", loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
        ax.set_title("Spatial (Leiden)")
    else:
        ax.text(0.5, 0.5, "Leiden unavailable", ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")

    ax = axes2[2]
    grid_filt.plot(ax=ax, column="n_genes", cmap=cmap, legend=True).axis("off")
    ax.set_title("n_genes")

    fig_dot, ax_dot = plt.subplots(1, 1, figsize=(20, 4), dpi=200)
    _plot_rank_genes_dotplot_matplotlib(adata_hex, n_genes=n_genes, cmap="Reds", ax=ax_dot)
    
    return fig1, fig2, fig_dot