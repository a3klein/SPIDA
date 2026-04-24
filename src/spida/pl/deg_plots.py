import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_dotplot(
    adata: ad.AnnData,
    summary_df: pd.DataFrame,
    celltype_col: str,
    gene_col: str | None = None,
    layer: str | None = None,
    celltypes: list[str] | None = None,
    palette: dict[str, str] | None = None,
    cmap: str = 'Reds',
    figsize: tuple | None = None,
    dpi: int = 150,
    ax=None,
) -> tuple:
    """
    General dotplot from summarize_scran_results or summarize_scanpy_results.

    Dot size = % expressing, dot color = mean expression. Genes are grouped on
    the x-axis by cell type of origin, with optional palette-colored tick labels
    and vertical separators between groups. Font sizes scale with figure size
    and data density.

    Parameters
    ----------
    adata : AnnData
    summary_df : pd.DataFrame
        Output of summarize_scran_results or summarize_scanpy_results.
        Must have a 'cell_type' column and a gene-list column.
    celltype_col : str
        obs column used for the y-axis (all unique values in adata).
    gene_col : str, optional
        Column in summary_df with gene lists. Auto-detected if None:
        scran  → first column matching 'top_*_markers';
        scanpy → 'top_upregulated'.
    layer : str, optional
        Layer for expression values. None → adata.X.
    celltypes : list[str], optional
        Subset of rows in summary_df to include. None → all rows.
    palette : dict[str, str], optional
        Mapping {cell_type: color} for coloring x-tick gene labels by their
        cell-type of origin. None → all gene labels drawn in black.
    cmap : str
        Colormap for mean expression (dot color).
    figsize : tuple, optional
        Auto-computed if None.
    dpi : int
    ax : matplotlib Axes, optional

    Returns
    -------
    fig, ax
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if gene_col is None:
        for col in summary_df.columns:
            if col.startswith('top_') and col.endswith('_markers'):
                gene_col = col
                break
        if gene_col is None and 'top_upregulated' in summary_df.columns:
            gene_col = 'top_upregulated'
        if gene_col is None:
            raise ValueError(
                "Cannot auto-detect gene column. Pass gene_col explicitly "
                "(e.g. 'top_upregulated' or 'top_10_markers')."
            )

    df = summary_df.copy()
    if celltypes is not None:
        df = df[df['cell_type'].isin(celltypes)]
    ct_keys = df['cell_type'].tolist()

    gene_ct_pairs = []
    seen = set()
    for _, row in df.iterrows():
        ct = row['cell_type']
        genes_list = row[gene_col]
        if not isinstance(genes_list, list):
            continue
        for gene in genes_list:
            if gene not in seen and gene in adata.var_names:
                gene_ct_pairs.append((gene, ct))
                seen.add(gene)

    if not gene_ct_pairs:
        print('No genes found in adata.var_names — check gene_col or celltypes.')
        return None, None

    genes       = [g for g, _ in gene_ct_pairs]
    gene_origin = [ct for _, ct in gene_ct_pairs]

    gene_idx = [adata.var_names.get_loc(g) for g in genes]
    X = adata.layers[layer] if layer is not None else adata.X
    X_sub = X[:, gene_idx].toarray() if sp.issparse(X) else np.asarray(X)[:, gene_idx]

    ct_labels = adata.obs[celltype_col].astype(str).values
    all_cts   = sorted(adata.obs[celltype_col].unique())
    n_g, n_c  = len(genes), len(all_cts)

    if figsize is None:
        figsize = (max(10, n_g * 0.45 + 2), max(4, n_c * 0.32 + 2))

    gene_fontsize  = float(np.clip(figsize[0] * 0.65 / max(n_g, 1) * 72 / 1.2, 5, 14))
    ct_fontsize    = float(np.clip(figsize[1] * 0.75 / max(n_c, 1) * 72 / 1.2, 5, 14))
    title_fontsize = float(np.clip(figsize[0] * 0.55, 8, 18))

    dot_x, dot_y, dot_s, dot_c = [], [], [], []
    for yi, cl in enumerate(reversed(all_cts)):
        mask = ct_labels == cl
        if mask.sum() == 0:
            continue
        x_sub     = X_sub[mask]
        pct_expr  = (x_sub > 0).mean(axis=0) * 100.0
        mean_expr = x_sub.mean(axis=0)
        for xi in range(n_g):
            dot_x.append(xi)
            dot_y.append(yi)
            dot_s.append(10 + pct_expr[xi] * 2.5)
            dot_c.append(mean_expr[xi])

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=False)
    else:
        fig = ax.figure

    sca = ax.scatter(
        dot_x, dot_y, s=dot_s, c=dot_c, cmap=cmap,
        vmin=0, vmax=np.percentile(dot_c, 99),
        edgecolors='k', linewidths=0.2, alpha=0.9, zorder=3,
    )

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=45, fontsize=gene_fontsize, ha='right')
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(list(reversed(all_cts)), fontsize=ct_fontsize)
    ax.set_xlim(-0.5, n_g - 0.5)
    ax.set_ylim(-0.5, n_c - 0.5)
    ax.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)

    for tick, ct in zip(ax.get_xticklabels(), gene_origin):
        tick.set_color(palette[ct] if palette and ct in palette else 'black')

    prev_ct = None
    for xi, ct in enumerate(gene_origin):
        if ct != prev_ct and xi > 0:
            ax.axvline(xi - 0.5, color='gray', lw=0.7, linestyle='--', alpha=0.5, zorder=1)
        prev_ct = ct

    ax.set_title(
        f'Top markers — {gene_col} ({n_g} genes, {len(ct_keys)} cell types)',
        fontsize=title_fontsize,
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='2%', pad=0.08)
    cbar_label_fs = float(np.clip(title_fontsize * 0.65, 6, 12))
    cbar_tick_fs  = float(np.clip(title_fontsize * 0.55, 5, 10))
    cb = fig.colorbar(sca, cax=cax)
    cb.set_label('Mean expression', fontsize=cbar_label_fs)
    cb.ax.tick_params(labelsize=cbar_tick_fs)

    origin_cts_shown = list(dict.fromkeys(gene_origin))

    if palette:
        ct_handles = [
            mpatches.Patch(color=palette.get(ct, 'black'), label=ct)
            for ct in origin_cts_shown
        ]
        leg_ct = ax.legend(
            handles=ct_handles, title='Gene origin',
            bbox_to_anchor=(1.05, 1.0), loc='upper left',
            bbox_transform=ax.transAxes,
            frameon=False, fontsize=max(6, ct_fontsize * 0.85),
            title_fontsize=max(7, ct_fontsize),
            ncols=max(1, len(ct_handles) // 25 + 1),
        )
        ax.add_artist(leg_ct)

    size_handles = [
        plt.scatter([], [], s=10 + p * 2.5, c='lightgray', edgecolors='k',
                    linewidths=0.3, label=f'{p}%')
        for p in [25, 50, 75, 100]
    ]
    ax.legend(
        handles=size_handles, title='% expressing',
        bbox_to_anchor=(1.05, 0.35), loc='upper left',
        bbox_transform=ax.transAxes,
        frameon=False, fontsize=max(6, ct_fontsize * 0.85),
        title_fontsize=max(7, ct_fontsize),
    )

    if own_fig:
        fig.tight_layout()

    return fig, ax


def plot_volcano(
    deseq2_results: pd.DataFrame,
    logfc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
    top_n_labels: int = 10,
    label_col: str = 'names',
    figsize: tuple = (8, 6),
    dpi: int = 150,
    ax=None,
) -> tuple:
    """
    Volcano plot for call_degs_pydeseq2 results.

    Parameters
    ----------
    deseq2_results : pd.DataFrame
        Output of call_degs_pydeseq2. Must have columns: log2FoldChange,
        padj, significant, and label_col.
    logfc_threshold : float
        Vertical threshold lines drawn at ±logfc_threshold.
    pval_threshold : float
        Horizontal threshold line drawn at -log10(pval_threshold).
    top_n_labels : int
        Number of top significant genes to label (ranked by |log2FC|).
    label_col : str
        Column to use for gene labels (default 'names').
    figsize : tuple
    dpi : int
    ax : matplotlib Axes, optional

    Returns
    -------
    fig, ax
    """
    try:
        from adjustText import adjust_text
        _has_adjusttext = True
    except ImportError:
        _has_adjusttext = False

    df = deseq2_results.copy()

    min_nonzero = df.loc[df['padj'] > 0, 'padj'].min()
    df['padj_clipped']   = df['padj'].clip(lower=min_nonzero)
    df['neg_log10_padj'] = -np.log10(df['padj_clipped'])

    sig    = df['significant']
    up     = sig & (df['log2FoldChange'] > 0)
    down   = sig & (df['log2FoldChange'] < 0)
    nonsig = ~sig

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure

    ax.scatter(
        df.loc[nonsig, 'log2FoldChange'], df.loc[nonsig, 'neg_log10_padj'],
        c='lightgray', s=12, alpha=0.6, linewidths=0, zorder=2, label='Not significant',
    )
    ax.scatter(
        df.loc[up, 'log2FoldChange'], df.loc[up, 'neg_log10_padj'],
        c='#d62728', s=18, alpha=0.8, linewidths=0.2, edgecolors='k', zorder=3,
        label=f'Up ({up.sum()})',
    )
    ax.scatter(
        df.loc[down, 'log2FoldChange'], df.loc[down, 'neg_log10_padj'],
        c='#1f77b4', s=18, alpha=0.8, linewidths=0.2, edgecolors='k', zorder=3,
        label=f'Down ({down.sum()})',
    )

    ax.axhline(-np.log10(pval_threshold), color='k', lw=0.8, linestyle='--', alpha=0.6)
    ax.axvline( logfc_threshold, color='k', lw=0.8, linestyle='--', alpha=0.6)
    ax.axvline(-logfc_threshold, color='k', lw=0.8, linestyle='--', alpha=0.6)

    top_genes = (
        df[sig]
        .assign(abs_lfc=df.loc[sig, 'log2FoldChange'].abs())
        .sort_values('abs_lfc', ascending=False)
        .head(top_n_labels)
    )

    if _has_adjusttext:
        texts = [
            ax.text(row['log2FoldChange'], row['neg_log10_padj'],
                    row[label_col], fontsize=7, alpha=0.9)
            for _, row in top_genes.iterrows()
        ]
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            expand=(1.2, 1.4),
        )
    else:
        for _, row in top_genes.iterrows():
            ax.annotate(
                row[label_col],
                xy=(row['log2FoldChange'], row['neg_log10_padj']),
                xytext=(4, 2), textcoords='offset points',
                fontsize=7, alpha=0.9,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            )

    ax.set_xlabel('log₂ Fold Change', fontsize=11)
    ax.set_ylabel('-log₁₀(adjusted p-value)', fontsize=11)

    if 'comparison' in df.columns:
        title = df['comparison'].iloc[0].replace('_vs_', ' vs ')
    else:
        title = 'DESeq2 volcano'
    ax.set_title(title, fontsize=12)

    ax.legend(frameon=False, fontsize=9, loc='upper center')
    ax.grid(True, alpha=0.2, zorder=0)

    n_sig = sig.sum()
    ax.annotate(
        f'significant: {n_sig}    up: {up.sum()}    down: {down.sum()}',
        xy=(0.5, -0.12), xycoords='axes fraction',
        ha='center', va='top', fontsize=8,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
    )

    if own_fig:
        plt.tight_layout()

    return fig, ax
