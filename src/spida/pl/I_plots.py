from pathlib import Path
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

from ._scatteplot import plot_categorical, plot_continuous
from ._spatial import plot_spatial_continuous, plot_spatial_categorical


def plot_allcools_joint_embeddings(
    adata : ad.AnnData, 
    exp_name : str, 
    seg_name : str,
    donor : str,
    ref_cell_type_column : str = "all_annot",
    pdf_file=None,
): 
    plt.rcParams['axes.facecolor'] = 'white'

    fig, ax = plt.subplots(2, 2, figsize=(8, 8), dpi=200)

    plot_categorical(adata, cluster_col="integrated_leiden", coord_base='X_integrated_umap', show=False, ax=ax[0,0])
    ax[0,0].set_title("ALL - Integrated Leiden Clusters")
    plot_categorical(adata, cluster_col=ref_cell_type_column, coord_base='X_integrated_umap', coding=True, text_anno=True, show=False, ax=ax[0,1])
    ax[0,1].set_title("REF - Subclass annotations")

    # plot_categorical(adata, cluster_col=None, coord_base='X_integrated_umap', color='lightgrey', show=False, ax=ax[1,0])
    plot_categorical(adata[adata.obs['Modality'] == "query"], cluster_col="integrated_leiden", coord_base='X_integrated_umap', show=False, ax=ax[1,0])
    ax[1,0].set_title("MERFISH - Integrated Leiden Clusters")
    plot_categorical(adata[adata.obs['Modality'] == "query"], cluster_col="all_annot", coord_base='X_integrated_umap', text_anno=True, coding=True, show=False, ax=ax[1,1])
    ax[1,1].set_title("MERFISH - Subclass annotations")
    plt.suptitle(f"{exp_name} - {seg_name} - {donor}", fontsize=16)

    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()


def plot_allcools_spatial_annot(
    adata : ad.AnnData, 
    exp_name : str, 
    seg_name : str,
    donor : str,
    ref_cell_type_column : str = "Subclass",
    pdf_file=None,
): 
    plt.rcParams['axes.facecolor'] = 'white'

    fig, ax = plt.subplots(2, 2, figsize=(8, 8), dpi=200)

    plot_categorical(adata, cluster_col=f"allcools_{ref_cell_type_column}", coord_base='X_base_umap', text_anno=True, coding=True, show=False, ax=ax[0,0])
    ax[0,0].set_title("Subclass Annot")
    plot_continuous(adata, color_by=f"allcools_{ref_cell_type_column}_transfer_score", coord_base='X_base_umap', show=False, ax=ax[0,1], cmap='viridis', hue_norm=[0, 1])
    ax[0,1].set_title("Subclass Transfer Score")

    plot_categorical(adata, cluster_col=f"allcools_{ref_cell_type_column}", coord_base='spatial', show=False, ax=ax[1,0])
    ax[1,0].set_title("Subclass Annot")
    plot_continuous(adata, color_by=f"allcools_{ref_cell_type_column}_transfer_score", coord_base='spatial', show=False, ax=ax[1,1], cmap='viridis', hue_norm=[0, 1])
    ax[1,1].set_title("Subclass Transfer Score")
    plt.suptitle(f"{exp_name} - {seg_name} - {donor}", fontsize=16)

    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()


    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    sns.histplot(adata.obs[f'allcools_{ref_cell_type_column}_transfer_score'], bins=100, kde=True, ax=ax)
    ax.set_title("Subclass Transfer Score Distribution")
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    bars = ax.bar(
        adata.obs[f'allcools_{ref_cell_type_column}'].value_counts().index,
        adata.obs[f'allcools_{ref_cell_type_column}'].value_counts().values,
        color='lightblue', edgecolor='black', alpha=1
    )
    ax.bar_label(bars, padding=3, fontsize=8, color='black', weight='bold')
    ax.set_title("Subclass Counts")
    ax.set_xlabel("Subclass")
    ax.set_ylabel("Count")
    ax.set_ylim((0, ax.get_ylim()[1] * 1.05))
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()
    

def plot_allcools_c2c(
    adata_comb,
    adata_merf,
    ref_group,
    query_group,
    conf_mat,
    ref_cell_type_column : str = "supercluster_name",
    pdf_file=None,
): 
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    sns.heatmap(conf_mat, vmin=0.1, vmax=0.9, cmap='cividis', cbar=None, ax=ax, yticklabels=True)
    cumsum_c, cumsum_r = 0, 0
    for xx in ref_group.sort_values().unique():
        rlen, clen = (query_group==xx).sum(), (ref_group==xx).sum()
        ax.plot([cumsum_c, cumsum_c + clen], [cumsum_r, cumsum_r], c='w', linewidth=0.5)
        ax.plot([cumsum_c, cumsum_c + clen], [cumsum_r + rlen, cumsum_r + rlen], c='w', linewidth=0.5)
        ax.plot([cumsum_c, cumsum_c], [cumsum_r, cumsum_r + rlen], c='w', linewidth=0.5)
        ax.plot([cumsum_c + clen, cumsum_c + clen], [cumsum_r, cumsum_r + rlen], c='w', linewidth=0.5)
        cumsum_c += clen
        cumsum_r += rlen

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()

    fig, ax = plt.subplots(1, 3,figsize=(20, 6), dpi=200)
    plot_categorical(adata_comb, coord_base="integrated_umap", cluster_col=f"c2c_allcools_label_{ref_cell_type_column}", show=False, coding=True, text_anno=True, ax=ax[0])
    plot_categorical(adata_merf, coord_base="base_umap", cluster_col=f"c2c_allcools_label_{ref_cell_type_column}", show=False, coding=True, text_anno=True, ax=ax[1])
    plot_categorical(adata_merf, coord_base="spatial", cluster_col=f"c2c_allcools_label_{ref_cell_type_column}", show=False, coding=True, ax=ax[2])
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    vc = adata_merf.obs[f'c2c_allcools_label_{ref_cell_type_column}'].value_counts()
    bars = ax.bar(x = vc.index, height=vc.values, color='lightblue', edgecolor='black', alpha=1)
    # texts = ax.bar_label(bars, padding=3, fontsize=8, color='black', weight='bold')
    texts=[]
    for j,rect in enumerate(bars):
        left = rect.get_x()+0.5
        top = rect.get_y()+rect.get_height()
        texts.append(ax.text(left,top,'%i'%vc.values[j], ha='center', va='bottom', weight='bold', fontsize=8))
    # adjustText.adjust_text(
    #     texts, add_objects=bars, only_move='y+', ax=ax,
    #     arrowprops=dict(arrowstyle="-", color="k", alpha=1))
    ax.set_xticklabels(vc.index, rotation=45, ha='right')
    ax.set_title(f"{ref_cell_type_column} distribution (n = {vc.sum()} / {adata_merf.shape[0]})")
    ax.set_ylabel("Number of cells")
    ax.set_xlabel(f"Annotated {ref_cell_type_column}")
    plt.tight_layout()
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()

def plot_degs(
    adata : ad.AnnData,
    groupby : str,
    layer : str = 'volume_norm',
): 
    """
    For plotting DEGs after integration
    """
    # TODO: Implement the DEG plots as part of the sc plots utility package. 
    raise NotImplementedError("This function is not yet implemented.")

    
    import scanpy as sc
    sc.tl.rank_genes_groups(adata, layer=layer, groupby=groupby, method='wilcoxon', groups=groups)

