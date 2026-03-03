# author: Amit Klein - a3klein@ucsd.edu - methods generated with copilot 
import os
import sys
from pathlib import Path
import warnings
import logging

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)

# #### USAGE: 
# if 'Subclass_name' in adata.obs.columns:
#     print("=== Running DEG analysis at Subclass level ===")
#     deg_results_subclass = call_degs_by_celltype(
#         adata=adata,
#         celltype_col='Subclass_name',
#         min_cells=25,
#         logfc_threshold=0.25,
#         pval_threshold=0.05,
#         method='wilcoxon',
#         correction_method='benjamini-hochberg',
#         save_results=True,
#         output_dir='{filepath}'
#     )    
#     # Show summary
#     summary_subclass = summarize_deg_results(deg_results_subclass, top_n=20)
#     print("\nSubclass-level DEG Summary:")
#     print(summary_subclass)
#     summary_subclass.to_csv("{filepath}", index=False)

def call_degs_by_celltype(
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
    verbose: bool = False
) -> dict[str, pd.DataFrame]:
    """
    Perform differential gene expression analysis for each cell type against all others.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    celltype_col : str
        Column name in adata.obs containing cell type annotations
        Options: 'Neighborhood_name', 'Class_name', 'Subclass_name', 'Group_name'
    layer : str, optional
        Layer to use for expression data. If None, uses adata.X
    min_cells : int
        Minimum number of cells required per cell type
    min_genes : int
        Minimum number of genes to test
    logfc_threshold : float
        Minimum log fold change threshold
    pval_threshold : float
        P-value threshold for significance
    method : str
        Statistical test method ('wilcoxon', 't-test', 't-test_overestim_var')
    correction_method : str
        Multiple testing correction method ('benjamini-hochberg', 'bonferroni')
    n_genes : int, optional
        Number of top genes to return per cell type
    save_results : bool
        Whether to save results to files
    output_dir : str, optional
        Directory to save results
        
    Returns:
    --------
    dict
        Dictionary with cell type names as keys and DEG results as DataFrames
    """
    
    # Validate inputs
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Column '{celltype_col}' not found in adata.obs")
    
    # Filter out cell types with too few cells
    celltype_counts = adata.obs[celltype_col].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= min_cells].index.tolist()
    
    if verbose: 
        logger.info(f"Found {len(valid_celltypes)} cell types with >= {min_cells} cells:")
        for ct in valid_celltypes:
            logger.info(f"  {ct}: {celltype_counts[ct]} cells")
        logger.info("removed_celltypes: %s", 
                    [ct for ct in celltype_counts.index if ct not in valid_celltypes])

    if len(valid_celltypes) < 2:
        raise ValueError("Need at least 2 cell types with sufficient cells")
    
    # Filter adata to only include valid cell types
    adata_filtered = adata[adata.obs[celltype_col].isin(valid_celltypes)].copy()
    
    # Set up scanpy settings
    sc.settings.verbosity = 1
    
    # Store original groups for restoration
    original_groups = adata_filtered.obs.get('temp_groups', None)
    
    deg_results = {}
    
    try:
        for celltype in valid_celltypes:
            if verbose: logger.info(f"\nProcessing cell type: {celltype}")
            
            # Create binary grouping: current cell type vs all others
            adata_filtered.obs['temp_groups'] = adata_filtered.obs[celltype_col].apply(
                lambda x: celltype if x == celltype else 'others'
            )
            
            # Check if we have both groups
            group_counts = adata_filtered.obs['temp_groups'].value_counts()
            if len(group_counts) < 2:
                logger.info(f"  Skipping {celltype}: insufficient groups")
                continue

            if verbose: 
                logger.info(f"  {celltype}: {group_counts[celltype]} cells")
                logger.info(f"  others: {group_counts['others']} cells")
                    
            max_cells = min(max_cells, group_counts[celltype], group_counts['others'])
            keep_cells = set()
            if group_counts['others'] > max_cells: 
                keep_cells.update(adata_filtered.obs[adata_filtered.obs['temp_groups'] == "others"].sample(max_cells).index)
            else: 
                keep_cells.update(adata_filtered.obs[adata_filtered.obs['temp_groups'] == "others"].index)
            if group_counts[celltype] > max_cells:
                keep_cells.update(adata_filtered.obs[adata_filtered.obs['temp_groups'] == celltype].sample(max_cells).index)
            else: 
                keep_cells.update(adata_filtered.obs[adata_filtered.obs['temp_groups'] == celltype].index)
            adata_ds = adata_filtered[adata_filtered.obs.index.isin(keep_cells)]
            group_counts = adata_ds.obs['temp_groups'].value_counts()

            if verbose: 
                logger.info(f"  Total cells after filtering: {adata_ds.n_obs}")
                logger.info(f"\t  {celltype}: {group_counts[celltype]} cells")
                logger.info(f"\t  others: {group_counts['others']} cells")
                
                logger.info("Running differential expression analysis")
            # Perform differential expression
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
                        tie_correct=True
                    )
                    
                    # Extract results
                    result = sc.get.rank_genes_groups_df(
                        adata_ds, 
                        group=celltype,
                        pval_cutoff=1.0,  # Get all genes, filter later
                        log2fc_min=None   # Get all genes, filter later
                    )
                    
                    # Add multiple testing correction
                    if len(result) > 0:
                        if correction_method == 'benjamini-hochberg':
                            rejected, pvals_corrected, _, _ = multipletests(
                                result['pvals'], 
                                method='fdr_bh'
                            )
                            result['pvals_adj'] = pvals_corrected
                        elif correction_method == 'bonferroni':
                            rejected, pvals_corrected, _, _ = multipletests(
                                result['pvals'], 
                                method='bonferroni'
                            )
                            result['pvals_adj'] = pvals_corrected
                        else:
                            result['pvals_adj'] = result['pvals']
                        
                        # Add significance flags
                        result['significant'] = (
                            (result['pvals_adj'] < pval_threshold) & 
                            (abs(result['logfoldchanges']) > logfc_threshold)
                        )
                        
                        # Sort by adjusted p-value and log fold change
                        result = result.sort_values(['significant', 'pvals_adj', 'logfoldchanges'], 
                                                  ascending=[False, True, False])
                        
                        # Add cell type information
                        result['cell_type'] = celltype
                        result['comparison'] = f"{celltype}_vs_others"
                        
                        deg_results[celltype] = result
                        
                        n_sig = result['significant'].sum()
                        if verbose: logger.info(f"  Found {n_sig} significant DEGs")

                except Exception as e:
                    logger.error(f"  Error processing {celltype}: {str(e)}")
                    continue
    
    finally:
        # Clean up temporary column
        if 'temp_groups' in adata_filtered.obs.columns:
            adata_filtered.obs.drop('temp_groups', axis=1, inplace=True)
    
    # Save results if requested
    if save_results and output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save individual cell type results
        for celltype, result_df in deg_results.items():
            filename = f"DEGs_{celltype.replace('/', '_').replace(' ', '_')}_vs_others.csv"
            result_df.to_csv(output_path / filename, index=False)
        
        # Save combined results
        if deg_results:
            combined_df = pd.concat(deg_results.values(), ignore_index=True)
            combined_filename = f"DEGs_all_celltypes_{celltype_col}.csv"
            combined_df.to_csv(output_path / combined_filename, index=False)
            if verbose: logger.info(f"\nSaved results to {output_path}")

    return deg_results


def summarize_deg_results(deg_results: dict[str, pd.DataFrame], 
                         top_n: int = 10) -> pd.DataFrame:
    """
    Create a summary table of DEG results.
    
    Parameters:
    -----------
    deg_results : dict
        Results from call_degs_by_celltype
    top_n : int
        Number of top genes to show per cell type
        
    Returns:
    --------
    pd.DataFrame
        Summary table
    """
    summary_data = []
    
    for celltype, result_df in deg_results.items():
        n_total = len(result_df)
        n_significant = result_df['significant'].sum()
        n_upregulated = ((result_df['significant']) & 
                        (result_df['logfoldchanges'] > 0)).sum()
        n_downregulated = ((result_df['significant']) & 
                          (result_df['logfoldchanges'] < 0)).sum()
        
        # Get top upregulated genes
        top_up = result_df[
            (result_df['significant']) & 
            (result_df['logfoldchanges'] > 0)
        ].head(top_n)['names'].tolist()
        
        # Get top downregulated genes
        top_down = result_df[
            (result_df['significant']) & 
            (result_df['logfoldchanges'] < 0)
        ].head(top_n)['names'].tolist()
        
        summary_data.append({
            'cell_type': celltype,
            'total_genes_tested': n_total,
            'significant_genes': n_significant,
            'upregulated': n_upregulated,
            'downregulated': n_downregulated,
            'top_upregulated': top_up,  # Show all
            'top_downregulated': top_down  # Show all
        })
    
    return pd.DataFrame(summary_data)



# Plotting Functionality by copilot, not the most informative in my opinion
def plot_deg_summary(deg_results: dict[str, pd.DataFrame], 
                     title: str = "DEG Analysis Summary",
                     figsize: tuple = (12, 8)):
    """
    Create summary plots for DEG results.
    """
    summary_df = summarize_deg_results(deg_results)
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle(title, fontsize=16)
    
    # Number of significant genes per cell type
    axes[0, 0].bar(range(len(summary_df)), summary_df['significant_genes'])
    axes[0, 0].set_title('Significant DEGs per Cell Type')
    axes[0, 0].set_ylabel('Number of DEGs')
    axes[0, 0].set_xticks(range(len(summary_df)))
    axes[0, 0].set_xticklabels(summary_df['cell_type'], rotation=45, ha='right')
    
    # Up vs Down regulated
    x = range(len(summary_df))
    axes[0, 1].bar(x, summary_df['upregulated'], label='Upregulated', alpha=0.7)
    axes[0, 1].bar(x, -summary_df['downregulated'], label='Downregulated', alpha=0.7)
    axes[0, 1].set_title('Up vs Down Regulated DEGs')
    axes[0, 1].set_ylabel('Number of DEGs')
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(summary_df['cell_type'], rotation=45, ha='right')
    axes[0, 1].legend()
    axes[0, 1].axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # Total genes tested
    axes[1, 0].bar(range(len(summary_df)), summary_df['total_genes_tested'])
    axes[1, 0].set_title('Total Genes Tested per Cell Type')
    axes[1, 0].set_ylabel('Number of Genes')
    axes[1, 0].set_xticks(range(len(summary_df)))
    axes[1, 0].set_xticklabels(summary_df['cell_type'], rotation=45, ha='right')
    
    # Percentage significant
    pct_sig = (summary_df['significant_genes'] / summary_df['total_genes_tested']) * 100
    axes[1, 1].bar(range(len(summary_df)), pct_sig)
    axes[1, 1].set_title('Percentage of Significant DEGs')
    axes[1, 1].set_ylabel('Percentage (%)')
    axes[1, 1].set_xticks(range(len(summary_df)))
    axes[1, 1].set_xticklabels(summary_df['cell_type'], rotation=45, ha='right')
    
    plt.tight_layout()
    plt.show()
    
    return summary_df


def plot_volcano(deg_result: pd.DataFrame, 
                 celltype: str,
                 logfc_threshold: float = 0.25,
                 pval_threshold: float = 0.05,
                 top_n: int = 10):
    """
    Create a volcano plot for a specific cell type's DEG results.
    """
    plt.figure(figsize=(10, 8))
    
    # Calculate -log10(p-value)
    neg_log_pval = -np.log10(deg_result['pvals_adj'])
    
    # Create scatter plot
    plt.scatter(deg_result['logfoldchanges'], neg_log_pval, 
               c='lightgray', alpha=0.6, s=20)
    
    # Highlight significant genes
    significant = deg_result['significant']
    plt.scatter(deg_result.loc[significant, 'logfoldchanges'], 
               neg_log_pval[significant], 
               c='red', alpha=0.8, s=30)
    
    # Add threshold lines
    plt.axhline(y=-np.log10(pval_threshold), color='blue', linestyle='--', alpha=0.7)
    plt.axvline(x=logfc_threshold, color='blue', linestyle='--', alpha=0.7)
    plt.axvline(x=-logfc_threshold, color='blue', linestyle='--', alpha=0.7)
    
    # Label top genes
    top_genes = deg_result[significant].head(top_n)
    for idx, row in top_genes.iterrows():
        plt.annotate(row['names'], 
                    (row['logfoldchanges'], -np.log10(row['pvals_adj'])),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(Adjusted P-value)')
    plt.title(f'Volcano Plot: {celltype} vs Others')
    plt.grid(True, alpha=0.3)
    
    # Add text box with statistics
    n_sig = significant.sum()
    n_up = ((deg_result['significant']) & (deg_result['logfoldchanges'] > 0)).sum()
    n_down = ((deg_result['significant']) & (deg_result['logfoldchanges'] < 0)).sum()
    
    textstr = f'Significant DEGs: {n_sig}\nUpregulated: {n_up}\nDownregulated: {n_down}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.show()