# author: Amit Klein - a3klein@ucsd.edu
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

try:
    import scranpy as _scranpy
except ImportError:
    _scranpy = None

try:
    from pydeseq2.dds import DeseqDataSet as _DeseqDataSet
    from pydeseq2.default_inference import DefaultInference as _DefaultInference
    from pydeseq2.ds import DeseqStats as _DeseqStats
except ImportError:
    _DeseqDataSet = _DefaultInference = _DeseqStats = None

from spida.utilities.ad_utils import _remove_small_clusters, _downsample_ref_clusters


# ---------------------------------------------------------------------------
# Scanpy-based DEGs  (one-vs-rest per cell type)
# ---------------------------------------------------------------------------

def call_degs_scanpy(
    adata: ad.AnnData,
    celltype_col: str,
    layer: str | None = None,
    min_cells: int = 10,
    max_cells: int = 50000,
    min_genes: int = 100,
    logfc_threshold: float = 0.25,
    pval_threshold: float = 0.05,
    method: str = 'wilcoxon',
    correction_method: str = 'benjamini-hochberg',
    n_genes: int | None = None,
    save_results: bool = True,
    output_dir: str | None = None,
    uns_key: str | None = None,
    verbose: bool = False,
) -> dict[str, pd.DataFrame]:
    """
    Perform differential gene expression analysis for each cell type against
    all others using scanpy's rank_genes_groups.

    Parameters
    ----------
    adata : AnnData
    celltype_col : str
        obs column with cell type annotations.
    layer : str, optional
        Layer to use. None → adata.X.
    min_cells : int
        Minimum cells per cell type.
    max_cells : int
        Maximum cells per group (downsamples if exceeded).
    min_genes : int
        Minimum genes to test.
    logfc_threshold : float
        Minimum absolute log fold change for significance.
    pval_threshold : float
        Adjusted p-value threshold for significance.
    method : str
        Test method ('wilcoxon', 't-test', 't-test_overestim_var').
    correction_method : str
        Multiple testing correction ('benjamini-hochberg', 'bonferroni').
    n_genes : int, optional
        Number of top genes to return per cell type. None → all genes.
    save_results : bool
    output_dir : str, optional
    uns_key : str, optional
        If provided, stacked results are stored in adata.uns[uns_key] in a
        format serializable to h5ad and zarr (no list-valued columns).

    Returns
    -------
    dict[str, pd.DataFrame]
        Cell type names → DEG result DataFrames.
    """
    import scanpy as sc

    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Column '{celltype_col}' not found in adata.obs")

    celltype_counts = adata.obs[celltype_col].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= min_cells].index.tolist()

    if verbose:
        logger.info(f"Found {len(valid_celltypes)} cell types with >= {min_cells} cells")
        removed = [ct for ct in celltype_counts.index if ct not in valid_celltypes]
        if removed:
            logger.info(f"Removed cell types: {removed}")

    if len(valid_celltypes) < 2:
        raise ValueError("Need at least 2 cell types with sufficient cells")

    adata_filtered = adata[adata.obs[celltype_col].isin(valid_celltypes)].copy()
    sc.settings.verbosity = 1
    deg_results = {}

    try:
        for celltype in valid_celltypes:
            if verbose:
                logger.info(f"Processing cell type: {celltype}")

            adata_filtered.obs['temp_groups'] = adata_filtered.obs[celltype_col].apply(
                lambda x: celltype if x == celltype else 'others'
            )
            group_counts = adata_filtered.obs['temp_groups'].value_counts()
            if len(group_counts) < 2:
                continue

            _max = min(max_cells, group_counts[celltype], group_counts['others'])
            keep_cells = set()
            for grp in [celltype, 'others']:
                subset = adata_filtered.obs[adata_filtered.obs['temp_groups'] == grp]
                keep_cells.update(
                    subset.sample(_max).index if len(subset) > _max else subset.index
                )
            adata_ds = adata_filtered[adata_filtered.obs.index.isin(keep_cells)]

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    sc.tl.rank_genes_groups(
                        adata_ds,
                        groupby='temp_groups',
                        groups=[celltype],
                        reference='others',
                        method=method,
                        use_raw=False,
                        layer=layer,
                        n_genes=adata_ds.n_vars if n_genes is None else n_genes,
                        tie_correct=True,
                    )
                    result = sc.get.rank_genes_groups_df(
                        adata_ds, group=celltype, pval_cutoff=1.0, log2fc_min=None
                    )
                    if len(result) > 0:
                        method_map = {'benjamini-hochberg': 'fdr_bh', 'bonferroni': 'bonferroni'}
                        mc = method_map.get(correction_method)
                        if mc:
                            _, pvals_adj, _, _ = multipletests(result['pvals'], method=mc)
                            result['pvals_adj'] = pvals_adj
                        else:
                            result['pvals_adj'] = result['pvals']

                        result['significant'] = (
                            (result['pvals_adj'] < pval_threshold) &
                            (result['logfoldchanges'].abs() > logfc_threshold)
                        )
                        result = result.sort_values(
                            ['significant', 'pvals_adj', 'logfoldchanges'],
                            ascending=[False, True, False],
                        )
                        result['cell_type'] = celltype
                        result['comparison'] = f'{celltype}_vs_others'
                        deg_results[celltype] = result

                        if verbose:
                            logger.info(f"  Found {result['significant'].sum()} significant DEGs")
                except Exception as e:
                    logger.error(f"Error processing {celltype}: {e}")
                    continue
    finally:
        if 'temp_groups' in adata_filtered.obs.columns:
            adata_filtered.obs.drop('temp_groups', axis=1, inplace=True)

    if deg_results and uns_key is not None:
        adata.uns[uns_key] = pd.concat(deg_results.values(), ignore_index=True)

    if save_results and output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        for celltype, df in deg_results.items():
            fname = f"DEGs_{celltype.replace('/', '_').replace(' ', '_')}_vs_others.csv"
            df.to_csv(out / fname, index=False)
        if deg_results:
            pd.concat(deg_results.values(), ignore_index=True).to_csv(
                out / f'DEGs_all_celltypes_{celltype_col}.csv', index=False
            )

    return deg_results


def call_degs_by_celltype(*args, **kwargs) -> dict[str, pd.DataFrame]:
    """Deprecated. Use call_degs_scanpy instead."""
    warnings.warn(
        "call_degs_by_celltype is deprecated; use call_degs_scanpy instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return call_degs_scanpy(*args, **kwargs)


def summarize_scanpy_results(
    deg_results: dict[str, pd.DataFrame],
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Summary table of scanpy DEG results (output of call_degs_scanpy).

    Parameters
    ----------
    deg_results : dict
        Output of call_degs_scanpy.
    top_n : int
        Number of top genes to list per cell type.

    Returns
    -------
    pd.DataFrame
    """
    rows = []
    for celltype, df in deg_results.items():
        n_sig  = df['significant'].sum()
        n_up   = ((df['significant']) & (df['logfoldchanges'] > 0)).sum()
        n_down = ((df['significant']) & (df['logfoldchanges'] < 0)).sum()
        top_up = (
            df[(df['significant']) & (df['logfoldchanges'] > 0)]
            .head(top_n)['names'].tolist()
        )
        top_down = (
            df[(df['significant']) & (df['logfoldchanges'] < 0)]
            .head(top_n)['names'].tolist()
        )
        rows.append({
            'cell_type':          celltype,
            'total_genes_tested': len(df),
            'significant_genes':  n_sig,
            'upregulated':        n_up,
            'downregulated':      n_down,
            'top_upregulated':    top_up,
            'top_downregulated':  top_down,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# scranpy-based marker genes  (pairwise all-vs-all)
# ---------------------------------------------------------------------------

def call_degs_scran(
    adata: ad.AnnData,
    celltype_col: str,
    layer: str | None = None,
    block_col: str | None = None,
    min_cells: int = 10,
    max_cells_per_cluster: int | None = None,
    num_threads: int = 4,
    save_results: bool = False,
    output_dir: str | None = None,
    uns_key: str | None = None,
    verbose: bool = False,
) -> dict[str, pd.DataFrame]:
    """
    Marker gene detection using scranpy's score_markers (pairwise all-vs-all).

    Parameters
    ----------
    adata : AnnData
    celltype_col : str
        obs column with cell type labels.
    layer : str, optional
        Layer to use. None → adata.X. Log-normalised counts recommended
        (e.g. 'volume_norm'). Sparse matrices are passed through as-is.
    block_col : str, optional
        obs column for batch/donor blocking (e.g. 'donor').
    min_cells : int
        Cell types with fewer cells are excluded.
    max_cells_per_cluster : int, optional
        If set, downsample each cluster to this many cells before scoring.
    num_threads : int
    save_results : bool
    output_dir : str, optional
    uns_key : str, optional
        If provided, stacked results are stored in adata.uns[uns_key] in a
        format serializable to h5ad and zarr (no list-valued columns).
    verbose : bool

    Returns
    -------
    dict[str, pd.DataFrame]
        Keys are cell type names. Each DataFrame has one row per gene with
        columns for auc, cohens_d, delta_mean, delta_detected
        (each with min / mean / median / max / min_rank across all pairwise
        comparisons).
    """
    if _scranpy is None:
        raise ImportError("scranpy is required for call_degs_scran. Install it first.")
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Column '{celltype_col}' not found in adata.obs")

    adata_sub = _remove_small_clusters(adata, col=celltype_col, min_cells=min_cells)
    if max_cells_per_cluster is not None:
        adata_sub = _downsample_ref_clusters(
            adata_sub, col=celltype_col, max_cells=max_cells_per_cluster
        )

    if verbose:
        n_cts = adata_sub.obs[celltype_col].nunique()
        print(
            f'Running score_markers: {adata_sub.n_obs} cells, '
            f'{adata_sub.n_vars} genes, {n_cts} cell types'
        )
        if max_cells_per_cluster:
            print(f'  (downsampled to ≤{max_cells_per_cluster} cells per cluster)')

    X      = adata_sub.layers[layer] if layer is not None else adata_sub.X
    groups = adata_sub.obs[celltype_col].values
    block  = adata_sub.obs[block_col].values if block_col is not None else None

    raw = _scranpy.score_markers(X.T, groups=groups, block=block, num_threads=num_threads)

    STAT_KEYS    = ['auc', 'cohens_d', 'delta_mean', 'delta_detected']
    SUMMARY_COLS = ['min', 'mean', 'median', 'max', 'min_rank']
    gene_names   = adata_sub.var_names.tolist()
    valid_cts    = list(raw['group_ids'])
    deg_results  = {}

    for ct in valid_cts:
        record = {'names': gene_names}
        for stat in STAT_KEYS:
            if stat not in raw.keys():
                continue
            stat_df = raw[stat][ct].to_pandas()
            for col in SUMMARY_COLS:
                if col in stat_df.columns:
                    record[f'{stat}_{col}'] = stat_df[col].values
        df = pd.DataFrame(record)
        df['cell_type'] = ct
        if 'auc_mean' in df.columns:
            df = df.sort_values('auc_mean', ascending=False)
        deg_results[ct] = df.reset_index(drop=True)

    if uns_key is not None:
        adata.uns[uns_key] = pd.concat(deg_results.values(), ignore_index=True)

    if save_results and output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        for ct, df in deg_results.items():
            fname = f"scran_{ct.replace('/', '_').replace(' ', '_')}.csv"
            df.to_csv(out / fname, index=False)
        pd.concat(deg_results.values(), ignore_index=True).to_csv(
            out / f'scran_all_{celltype_col}.csv', index=False
        )
        if verbose:
            print(f'Results saved to {out}')

    return deg_results


def summarize_scran_results(
    scran_results: dict[str, pd.DataFrame],
    top_n: int = 10,
    auc_threshold: float = 0.7,
    min_rank_threshold: int = 50,
) -> pd.DataFrame:
    """
    Summary table of scranpy marker gene results (output of call_degs_scran).

    Selects top markers per cell type using auc_mean >= auc_threshold and
    auc_min_rank <= min_rank_threshold, ranked by min_rank then auc_mean.

    Parameters
    ----------
    scran_results : dict
        Output of call_degs_scran.
    top_n : int
        Max markers to report per cell type.
    auc_threshold : float
        Minimum mean AUC (0.5–1.0).
    min_rank_threshold : int
        Maximum min_rank across pairwise comparisons.

    Returns
    -------
    pd.DataFrame
    """
    rows = []
    for ct, df in scran_results.items():
        has_auc      = 'auc_mean' in df.columns and 'auc_min_rank' in df.columns
        has_cohens_d = 'cohens_d_mean' in df.columns

        if has_auc:
            top = (
                df.loc[
                    (df['auc_mean'] >= auc_threshold) &
                    (df['auc_min_rank'] <= min_rank_threshold)
                ]
                .sort_values(['auc_min_rank', 'auc_mean'], ascending=[True, False])
                .head(top_n)
            )
            top_genes = top['names'].tolist()
            n_markers = len(top)
        else:
            top_genes = df['names'].head(top_n).tolist()
            n_markers = None

        row = {
            'cell_type':       ct,
            'n_genes_tested':  len(df),
            f'n_markers (auc>={auc_threshold}, rank<={min_rank_threshold})': n_markers,
            f'top_{top_n}_markers': top_genes,
        }
        if has_auc and not top.empty:
            row['mean_auc_top_markers'] = top['auc_mean'].mean().round(3)
        if has_cohens_d and not top.empty:
            row['mean_cohens_d_top_markers'] = top['cohens_d_mean'].mean().round(3)
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# pydeseq2 pseudobulk DESeq2  (two-group comparison)
# ---------------------------------------------------------------------------

def call_degs_pydeseq2(
    adata: ad.AnnData,
    celltype_col: str,
    group1: str,
    group2: str,
    sample_col: str,
    layer: str | None = None,
    design_formula: str | None = None,
    min_samples: int = 3,
    lfc_shrink: bool = True,
    pval_threshold: float = 0.05,
    logfc_threshold: float = 0.5,
    n_cpus: int = 1,
    save_results: bool = False,
    output_dir: str | None = None,
    uns_key: str | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Pseudobulk DESeq2 comparison between two specific cell types.

    Cells are summed per sample_col × group into pseudobulk count profiles,
    then DESeq2 is run comparing group1 vs group2 (group2 is the reference).
    Requires pydeseq2 >= 0.5.

    Parameters
    ----------
    adata : AnnData
    celltype_col : str
        obs column with cell type labels.
    group1 : str
        Numerator cell type (positive LFC = higher in group1).
    group2 : str
        Reference cell type.
    sample_col : str
        obs column for pseudobulk aggregation (e.g. 'donor').
    layer : str, optional
        Layer with raw integer counts. None → adata.X. 'counts' recommended.
    design_formula : str, optional
        Patsy-style DESeq2 formula. Default: '~ group'.
        Covariates must be obs columns, e.g. '~ brain_region + group'.
    min_samples : int
        Min donors/samples required per group; raises ValueError if not met.
    lfc_shrink : bool
        Apply apeGLM LFC shrinkage after fitting. Default True.
    pval_threshold, logfc_threshold : float
        Thresholds for the 'significant' flag (uses padj and log2FoldChange).
    n_cpus : int
    save_results : bool
    output_dir : str, optional
    uns_key : str, optional
        If provided, results DataFrame is stored in adata.uns[uns_key].
    verbose : bool

    Returns
    -------
    pd.DataFrame
        DESeq2 results: names, baseMean, log2FoldChange, lfcSE, stat,
        pvalue, padj, significant, group1, group2, comparison.
    """
    if _DeseqDataSet is None:
        raise ImportError("pydeseq2 is required for call_degs_pydeseq2. Install it first.")

    for g in [group1, group2]:
        if g not in adata.obs[celltype_col].values:
            raise ValueError(f"'{g}' not found in adata.obs['{celltype_col}']")

    if design_formula is None:
        design_formula = '~ group'

    formula_terms = [t.strip() for t in design_formula.replace('~', '').split('+')]
    covariate_cols = [
        t for t in formula_terms if t != 'group' and t in adata.obs.columns
    ]

    mask      = adata.obs[celltype_col].isin([group1, group2])
    adata_sub = adata[mask]

    X = adata_sub.layers[layer] if layer is not None else adata_sub.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(float)

    obs         = adata_sub.obs.reset_index(drop=True)
    samples_arr = obs[sample_col].astype(str).values
    groups_arr  = obs[celltype_col].astype(str).values

    pb_counts, pb_meta, pb_index = [], [], []
    for samp in np.unique(samples_arr):
        for grp in [group1, group2]:
            cell_mask = (samples_arr == samp) & (groups_arr == grp)
            if cell_mask.sum() == 0:
                continue
            pb_counts.append(X[cell_mask].sum(axis=0))
            meta_row = {'sample_id': samp, 'group': grp}
            for cov in covariate_cols:
                vals = obs.loc[cell_mask, cov]
                meta_row[cov] = vals.mode().iloc[0]
            pb_meta.append(meta_row)
            pb_index.append(f'{samp}__{grp}')

    counts_df = pd.DataFrame(
        np.round(np.array(pb_counts)).astype(int),
        columns=adata_sub.var_names,
        index=pb_index,
    )
    meta_df = pd.DataFrame(pb_meta, index=pb_index)
    meta_df['group'] = pd.Categorical(meta_df['group'], categories=[group2, group1])

    per_group = meta_df.groupby('group', observed=True).size()
    if verbose:
        print('Pseudobulk samples per group:')
        print(per_group.to_string())
    if (per_group < min_samples).any():
        raise ValueError(
            f'Insufficient samples for DESeq2 (min_samples={min_samples}): '
            f'{per_group.to_dict()}'
        )

    nonzero   = counts_df.columns[counts_df.sum(axis=0) > 0]
    counts_df = counts_df[nonzero]
    if verbose:
        print(f'Genes after zero-filtering: {len(nonzero)} / {adata_sub.n_vars}')

    inference = _DefaultInference(n_cpus=n_cpus)
    dds = _DeseqDataSet(
        counts=counts_df,
        metadata=meta_df,
        design=design_formula,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    stat_res = _DeseqStats(
        dds,
        contrast=['group', group1, group2],
        inference=inference,
    )
    stat_res.summary()

    if lfc_shrink:
        lfc_cols = stat_res.LFC.columns.tolist()
        if verbose:
            print(f'Available LFC coefficients: {lfc_cols}')
        coeff = next(
            (c for c in lfc_cols if 'group' in c.lower() and 'intercept' not in c.lower()),
            None,
        )
        if coeff is None:
            warnings.warn('Could not find group coefficient for LFC shrinkage; skipping.')
        else:
            if verbose:
                print(f'Applying LFC shrinkage with coeff: {coeff}')
            try:
                stat_res.lfc_shrink(coeff=coeff)
            except Exception as e:
                warnings.warn(f'LFC shrinkage failed ({e}); returning unshrunk estimates.')

    results = stat_res.results_df.copy()
    results.index.name = 'names'
    results = results.reset_index()
    results['significant'] = (
        (results['padj'] < pval_threshold) &
        (results['log2FoldChange'].abs() > logfc_threshold)
    )
    results['group1']     = group1
    results['group2']     = group2
    results['comparison'] = f'{group1}_vs_{group2}'
    results = results.sort_values(
        ['significant', 'padj', 'log2FoldChange'],
        ascending=[False, True, False],
    ).reset_index(drop=True)

    if uns_key is not None:
        adata.uns[uns_key] = results

    if uns_key is not None:
        adata.uns[uns_key] = results

    if save_results and output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        fname = f'deseq2_{group1}_vs_{group2}.csv'.replace('/', '_').replace(' ', '_')
        results.to_csv(out / fname, index=False)
        if verbose:
            print(f'Saved to {out / fname}')

    return results


# ---------------------------------------------------------------------------
# Reconstruct summaries from adata.uns (after h5ad / zarr round-trip)
# ---------------------------------------------------------------------------

def summarize_scran_results_from_uns(
    adata: ad.AnnData,
    uns_key: str = 'scran_markers',
    top_n: int = 10,
    auc_threshold: float = 0.7,
    min_rank_threshold: int = 50,
) -> pd.DataFrame:
    """
    Reconstruct a summary table from scran marker results stored in adata.uns.

    Equivalent to calling summarize_scran_results on the original dict output
    of call_degs_scran. Works after h5ad / zarr serialization round-trips.

    Parameters
    ----------
    adata : AnnData
    uns_key : str
        Key in adata.uns where the stacked DataFrame was stored.
    top_n, auc_threshold, min_rank_threshold
        Passed through to summarize_scran_results.
    """
    if uns_key not in adata.uns:
        raise KeyError(f"'{uns_key}' not found in adata.uns")
    stacked = adata.uns[uns_key].copy()
    results = {
        ct: grp.reset_index(drop=True)
        for ct, grp in stacked.groupby('cell_type', sort=False)
    }
    return summarize_scran_results(
        results, top_n=top_n,
        auc_threshold=auc_threshold, min_rank_threshold=min_rank_threshold,
    )


def summarize_scanpy_results_from_uns(
    adata: ad.AnnData,
    uns_key: str = 'scanpy_markers',
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Reconstruct a summary table from scanpy DEG results stored in adata.uns.

    Equivalent to calling summarize_scanpy_results on the original dict output
    of call_degs_scanpy. Works after h5ad / zarr serialization round-trips.
    The 'significant' column is cast to bool if it was stored as int (h5ad
    stores bool as uint8).

    Parameters
    ----------
    adata : AnnData
    uns_key : str
        Key in adata.uns where the stacked DataFrame was stored.
    top_n : int
    """
    if uns_key not in adata.uns:
        raise KeyError(f"'{uns_key}' not found in adata.uns")
    stacked = adata.uns[uns_key].copy()
    stacked['significant'] = stacked['significant'].astype(bool)
    results = {
        ct: grp.reset_index(drop=True)
        for ct, grp in stacked.groupby('cell_type', sort=False)
    }
    return summarize_scanpy_results(results, top_n=top_n)
