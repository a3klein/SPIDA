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
    pdf_file=None,
): 
    plt.rcParams['axes.facecolor'] = 'white'

    fig, ax = plt.subplots(2, 2, figsize=(8, 8), dpi=200)

    plot_categorical(adata, cluster_col="integrated_leiden", coord_base='X_integrated_umap', show=False, ax=ax[0,0])
    ax[0,0].set_title("Integrated Leiden Clusters")
    plot_categorical(adata, cluster_col="all_annot", coord_base='X_integrated_umap', coding=True, text_anno=True, show=False, ax=ax[0,1])
    ax[0,1].set_title("All Subclass annotations")

    # plot_categorical(adata, cluster_col=None, coord_base='X_integrated_umap', color='lightgrey', show=False, ax=ax[1,0])
    plot_categorical(adata[adata.obs['Modality'] == "query"], cluster_col="integrated_leiden", coord_base='X_integrated_umap', show=False, ax=ax[1,0])
    ax[1,0].set_title("Integrated Leiden Clusters")
    plot_categorical(adata[adata.obs['Modality'] == "query"], cluster_col="all_annot", coord_base='X_integrated_umap', text_anno=True, coding=True, show=False, ax=ax[1,1])
    ax[1,1].set_title("All Subclass annotations")
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
    pdf_file=None,
): 
    plt.rcParams['axes.facecolor'] = 'white'

    fig, ax = plt.subplots(2, 2, figsize=(8, 8), dpi=200)

    plot_categorical(adata, cluster_col="allcools_Subclass", coord_base='X_base_umap', text_anno=True, coding=True, show=False, ax=ax[0,0])
    ax[0,0].set_title("Subclass Annot")
    plot_continuous(adata, color_by="allcools_Subclass_transfer_score", coord_base='X_base_umap', show=False, ax=ax[0,1], cmap='viridis', hue_norm=[0, 1])
    ax[0,1].set_title("Subclass Transfer Score")

    plot_categorical(adata, cluster_col="allcools_Subclass", coord_base='spatial', show=False, ax=ax[1,0])
    ax[1,0].set_title("Subclass Annot")
    plot_continuous(adata, color_by="allcools_Subclass_transfer_score", coord_base='spatial', show=False, ax=ax[1,1], cmap='viridis', hue_norm=[0, 1])
    ax[1,1].set_title("Subclass Transfer Score")
    plt.suptitle(f"{exp_name} - {seg_name} - {donor}", fontsize=16)

    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()


    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    sns.histplot(adata.obs['allcools_Subclass_transfer_score'], bins=100, kde=True, ax=ax)
    ax.set_title("Subclass Transfer Score Distribution")
    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    bars = ax.bar(
        adata.obs['allcools_Subclass'].value_counts().index,
        adata.obs['allcools_Subclass'].value_counts().values,
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

