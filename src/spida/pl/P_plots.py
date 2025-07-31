from pathlib import Path
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

from spida._constants import CELL_X, CELL_Y
from spida.pl import plot_scatter # TODO: depreceate this! 
from spida._utilities import _region_to_donor
from ._scatteplot import plot_categorical, plot_continuous


def plot_feature_distribution(
    df: pd.DataFrame,
    feature: str,
    feature_alias: str = None,
    min_val=None,
    max_val=None,
    num_bins: int = 100,
    pdfFile=None,
    ax=None,
):
    """
    attributes:
        feature(str) : the feature to plot
        feature_alias(str) : the name to plot for the feature
        min_val(float) : the minimum value for the feature
        max_val(float) : the maximum value for the feature
        pdfFile(str) : if not None, save the figure to this pdf
        ax(Matplotlib.Axes) : if not None, plot the figure on this ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)

    sns.histplot(df[feature].to_numpy(), bins=num_bins, ax=ax)
    if min_val:
        ax.axvline(min_val, color="red", linestyle="--")
    if max_val:
        ax.axvline(max_val, color="red", linestyle="--")

    ax.set_title(f"{feature} Distribution")
    if feature_alias:
        ax.set_title(f"{feature_alias} Distribution")

    ax.set_xlabel(feature)
    return ax


def plot_filt_cells(df_feature: pd.DataFrame, title="", pdfFile=None, ax=None):
    """
    attributes:
        pdfFile(str) : if not None, save the figure to this pdf
        ax(Matplotlib.Axes) : if not None, plot the figure on this ax
    """
    filt_palette = {True: "#507B58", False: "#540202"}
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)

    # Plotting the filtered cells
    sns.scatterplot(
        ax=ax,
        data=df_feature,
        x=CELL_X,
        y=CELL_Y,
        hue="pass_qc",
        s=1,
        palette=filt_palette,
    )
    ax.set_aspect("equal")
    ax.set_title(title)
    return ax


def plot_filtering(adata: ad.AnnData, exp_name: str, reg_name: str, prefix_name: str):
    """
    Plot filtering results (summary of filtering features with their cutoffs)

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the filtered data.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    """
    plt.rcParams["axes.facecolor"] = "white"
    fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(16, 8), constrained_layout=True)

    plot_feature_distribution(
        adata.obs,
        feature="volume",
        feature_alias="Cell Volume",
        min_val=adata.uns["cutoffs"]["volume_min"],
        max_val=adata.uns["cutoffs"]["volume_max"],
        ax=axes[0][0],
    )
    plot_feature_distribution(
        adata.obs,
        feature="nCount_RNA",
        feature_alias="RNA Count",
        min_val=adata.uns["cutoffs"]["n_count_min"],
        max_val=adata.uns["cutoffs"]["n_count_max"],
        ax=axes[0][1],
    )
    plot_feature_distribution(
        adata.obs,
        feature="nFeature_RNA",
        feature_alias="Num Unique RNAs",
        min_val=adata.uns["cutoffs"]["n_gene_min"],
        max_val=adata.uns["cutoffs"]["n_gene_max"],
        ax=axes[0][2],
    )
    plot_feature_distribution(
        adata.obs,
        feature="nBlank",
        feature_alias="Blank Count",
        min_val=adata.uns["cutoffs"]["n_blank_min"],
        max_val=adata.uns["cutoffs"]["n_blank_max"],
        ax=axes[1][0],
    )
    plot_feature_distribution(
        adata.obs,
        feature="nCount_RNA_per_Volume",
        feature_alias="Transcripts Per Cell Size",
        min_val=adata.uns["cutoffs"]["n_count_per_volume_min"],
        max_val=adata.uns["cutoffs"]["n_count_per_volume_max"],
        ax=axes[1][1],
    )
    plot_filt_cells(adata.obs, ax=axes[1][2], title="Filtered Cells")
    plt.suptitle(
        "QC_Summary - Exp %s; - Region %s; - Segmentation %s"
        % (exp_name, reg_name, prefix_name)
    )

    return fig, axes


def plot_setup(
    adata: ad.AnnData, exp_name: str, reg_name: str, prefix_name: str, pdf_file=None
):
    """
    Plot the setup results (summary of setup features with their cutoffs).

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the setup data.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    pdf_file (str, optional): Path to save the plot as a PDF file. Defaults to None.
    """
    plt.rcParams["axes.facecolor"] = "white"
    donor = _region_to_donor(reg_name)
    # Plot the UMAP embeddings
    plot_scatter(
        adata,
        title=f"{donor} - UMAP",
        colors=["volume", "nCount_RNA", "leiden"],
        ncols=3,
        coord_base="umap",
        pdf_file=pdf_file,
    )
    # Plot the TSNE embeddings
    plot_scatter(
        adata,
        title=f"{donor} - UMAP",
        colors=["volume", "nCount_RNA", "leiden"],
        ncols=3,
        coord_base="tsne",
        pdf_file=pdf_file,
    )
    # Plot the spatial coordinates
    plot_scatter(
        adata,
        title=f"{donor} - UMAP",
        colors=["volume", "nCount_RNA", "leiden"],
        ncols=3,
        coord_base="spatial",
        pdf_file=pdf_file,
    )


def plot_resolvi(
    adata: ad.AnnData, exp_name: str, reg_name: str, prefix_name: str, pdf_file=None
):
    """
    Plot the results of the RESOLVI clustering.
    Parameters:
    adata (anndata.AnnData): The AnnData object containing the RESOLVI clustering results.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    pdf_file (str, optional): Path to save the plot as a PDF file. Defaults to None.
    """
    plt.rcParams["axes.facecolor"] = "white"
    donor = _region_to_donor(reg_name)
    # Plot the UMAP embeddings
    fig, axes = plt.subplots(3, 3, figsize=(20, 20), dpi=300)
    plot_categorical(adata, coord_base="umap", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[0][0])
    axes[0][0].set_title(f"tsne - {donor} - leiden Clustering")
    plot_categorical(adata, coord_base="resolvi_umap", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[0][1])
    axes[0][1].set_title(f"tsne - {donor} - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="corr_umap", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[0][2])
    axes[0][2].set_title(f"tsne - {donor} - corrected expression leiden Clustering")
    # Plot the TSNE embeddings
    plot_categorical(adata, coord_base="tsne", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[1][0])
    axes[1][0].set_title(f"tsne - {donor} - leiden Clustering")
    plot_categorical(adata, coord_base="resolvi_tsne", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[1][1])
    axes[1][1].set_title(f"tsne - {donor} - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="corr_tsne", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[1][2])
    axes[1][2].set_title(f"tsne - {donor} - corrected expression leiden Clustering")
    # Plot the spatial coordinates
    plot_categorical(adata, coord_base="spatial", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[2][0])
    axes[2][0].set_title(f"{donor} - Leiden Clustering")
    plot_categorical(adata, coord_base="spatial", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[2][1])
    axes[2][1].set_title(f"{donor} - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="spatial", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[2][2])
    axes[2][2].set_title(f"{donor} - corrected expression leiden Clustering")
    
    if pdf_file:
        pdf_file.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
        plt.close()


def plot_doublets(adata, save_file, ax=None):
    arr = (
        adata.obsm["X_spatial"] if "X_spatial" in adata.obsm else adata.obsm["spatial"]
    )
    ds = min(200000 / adata.shape[0], 4)  # the point size to use for plotting

    filt_palette = {True: "#ff0000", False: "#cccccc"}

    if not ax:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=300, constrained_layout=True)
    sns.scatterplot(
        x=arr[:, 0],
        y=arr[:, 1],
        s=ds,
        hue=adata.obs["doublet_bool"],
        marker=".",
        ax=ax,
        palette=filt_palette,
        linewidth=0,
        alpha=0.9,
        legend=False,
        rasterized=True,
    )
    ax.set_aspect("equal")
    ax.set_title("Doublet Prediction")
    ax.axis("off")

    if save_file:
        fig.savefig(save_file, bbox_inches="tight")
        plt.close(fig)
    else:
        return ax

def plot_dataset(
    adata,
    save_path:str|Path=None,
    show=True,
    gene=None,
):
    """
    Plot the dataset.
    Parameters:
    adata (anndata.AnnData): The AnnData object containing the dataset.
    save_path (str, optional): Path to save the plot. If None, the plot is shown but not saved.
    show (bool, optional): Whether to show the plot. Defaults to True.
    """
    if isinstance(save_path, str):
        save_path = Path(save_path)
    plt.rcParams["axes.facecolor"] = "white"

    # plot region / donor / replicate metadata
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 10), dpi=300)
    plot_categorical(adata, coord_base="umap", cluster_col="region", text_anno=True, coding=True, show=False, ax=axes[0])
    axes[0].set_title(f"Region")
    plot_categorical(adata, coord_base="umap", cluster_col="donor", text_anno=False, coding=True, show=False, ax=axes[1])
    axes[1].set_title(f"Donor")
    plot_categorical(adata, coord_base="umap", cluster_col="replicate", text_anno=False, coding=True, show=False, ax=axes[2])
    axes[2].set_title(f"Replicate")
    if save_path:
        fig.savefig(save_path / "umap_meta.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)
    
    # plot tsne metadata
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 10), dpi=300)
    plot_categorical(adata, coord_base="tsne", cluster_col="region", text_anno=True, coding=True, show=False, ax=axes[0])
    axes[0].set_title(f"Region")
    plot_categorical(adata, coord_base="tsne", cluster_col="donor", text_anno=False, coding=True, show=False, ax=axes[1])
    axes[1].set_title(f"Donor")
    plot_categorical(adata, coord_base="tsne", cluster_col="replicate", text_anno=False, coding=True, show=False, ax=axes[2])
    axes[2].set_title(f"Replicate")
    if save_path:
        fig.savefig(save_path / "tsne_meta.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)    
    
    # plot umap / tsne / spatial coordinates 
    # Plot the UMAP embeddings
    fig, axes = plt.subplots(3, 3, figsize=(20, 20), dpi=300)
    plot_categorical(adata, coord_base="umap", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[0][0])
    axes[0][0].set_title(f"tsne - leiden Clustering")
    plot_categorical(adata, coord_base="resolvi_umap", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[0][1])
    axes[0][1].set_title(f"tsne - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="corr_umap", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[0][2])
    axes[0][2].set_title(f"tsne - corrected expression leiden Clustering")
    # Plot the TSNE embeddings
    plot_categorical(adata, coord_base="tsne", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[1][0])
    axes[1][0].set_title(f"tsne - leiden Clustering")
    plot_categorical(adata, coord_base="resolvi_tsne", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[1][1])
    axes[1][1].set_title(f"tsne - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="corr_tsne", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[1][2])
    axes[1][2].set_title(f"tsne - corrected expression leiden Clustering")
    # Plot the spatial coordinates
    plot_categorical(adata, coord_base="spatial", cluster_col="leiden", text_anno=False, coding=True, show=False, ax=axes[2][0])
    axes[2][0].set_title(f"Leiden Clustering")
    plot_categorical(adata, coord_base="spatial", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[2][1])
    axes[2][1].set_title(f"RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="spatial", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[2][2])
    axes[2][2].set_title(f"corrected expression leiden Clustering")
    if save_path:
        fig.savefig(save_path / "categorical.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)

    
    if gene is None: 
        gene = "MOBP"
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(20, 20), dpi=300)
    axes = axes.flatten()
    plot_continuous(adata, coord_base="tsne", color_by=gene, layer="generated_expression", ax=axes[0], show=False, title=f"tsne - {gene} - generated expression")
    plot_continuous(adata, coord_base="umap", color_by=gene, layer="generated_expression", ax=axes[1], show=False, title=f"umap - {gene} - generated expression")
    plot_continuous(adata, coord_base="spatial", color_by=gene, layer="generated_expression", ax=axes[2], show=False, title=f"spatial - {gene} - generated expression")
    plot_continuous(adata, coord_base="tsne", color_by=gene, layer="raw", ax=axes[3], show=False, title=f"tsne - {gene} - raw expression")
    plot_continuous(adata, coord_base="umap", color_by=gene, layer="raw", ax=axes[4], show=False, title=f"umap - {gene} - raw expression")
    plot_continuous(adata, coord_base="spatial", color_by=gene, layer="raw", ax=axes[5], show=False, title=f"spatial - {gene} - raw expression")
    if save_path:
        fig.savefig(save_path / "continuous.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)