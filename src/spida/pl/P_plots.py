from pathlib import Path
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

from spida._constants import CELL_X, CELL_Y
from spida.pl import plot_scatter # TODO: depreceate this! 
from spida._utilities import _region_to_donor
from ._scatteplot import plot_categorical, plot_continuous
from ._spatial import plot_spatial_continuous, plot_spatial_categorical


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
        adata.obs[adata.obs["pass_qc_pre"]] if "pass_qc_pre" in adata.obs.columns else adata.obs,
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
    fig, axes = plt.subplots(3, 3, figsize=(20, 20), dpi=300)
    plot_continuous(adata, coord_base="base_umap", color_by="volume", ax=axes[0][0], show=False, title="umap - Cell Volume")
    plot_continuous(adata, coord_base="base_umap", color_by="nCount_RNA", ax=axes[0][1], show=False, title="umap - Number of Transcripts")
    plot_categorical(adata, coord_base="base_umap", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=axes[0][2])
    axes[0][2].set_title(f"umap - Leiden Clusters")    

    plot_continuous(adata, coord_base="base_tsne", color_by="volume", ax=axes[1][0], show=False, title="tsne - Cell Volume")
    plot_continuous(adata, coord_base="base_tsne", color_by="nCount_RNA", ax=axes[1][1], show=False, title="tsne - Number of Transcripts")
    plot_categorical(adata, coord_base="base_tsne", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=axes[1][2])
    axes[1][2].set_title(f"tsne - Leiden Clusters")

    plot_continuous(adata, coord_base="spatial", color_by="volume", ax=axes[2][0], show=False, title="spatial - Cell Volume")
    plot_continuous(adata, coord_base="spatial", color_by="nCount_RNA", ax=axes[2][1], show=False, title="spatial - Number of Transcripts")
    plot_categorical(adata, coord_base="spatial", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=axes[2][2])
    axes[2][2].set_title(f"spatial - Leiden Clusters")

    if pdf_file:
        pdf_file.savefig(fig, dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close()

    # # Plot the UMAP embeddings
    # plot_scatter(
    #     adata,
    #     title=f"{donor} - UMAP",
    #     colors=["volume", "nCount_RNA", "leiden"],
    #     ncols=3,
    #     coord_base="umap",
    #     pdf_file=pdf_file,
    # )
    # # Plot the TSNE embeddings
    # plot_scatter(
    #     adata,
    #     title=f"{donor} - UMAP",
    #     colors=["volume", "nCount_RNA", "leiden"],
    #     ncols=3,
    #     coord_base="tsne",
    #     pdf_file=pdf_file,
    # )
    # # Plot the spatial coordinates
    # plot_scatter(
    #     adata,
    #     title=f"{donor} - UMAP",
    #     colors=["volume", "nCount_RNA", "leiden"],
    #     ncols=3,
    #     coord_base="spatial",
    #     pdf_file=pdf_file,
    # )


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
    axes[0][0].set_title(f"umap - {donor} - leiden Clustering")
    plot_categorical(adata, coord_base="resolvi_umap", cluster_col="resolvi_leiden", text_anno=False, coding=True, show=False, ax=axes[0][1])
    axes[0][1].set_title(f"umap - {donor} - RESOLVI leiden Clustering")
    plot_categorical(adata, coord_base="corr_umap", cluster_col="corr_leiden", text_anno=False, coding=True, show=False, ax=axes[0][2])
    axes[0][2].set_title(f"umap - {donor} - corrected expression leiden Clustering")
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

# TODO: Reconcile with salk-ucsd differences
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

    experiments = adata.obs['brain_region'].unique()
    donors = adata.obs['donor'].unique()
    replicates = adata.obs['replicate'].unique()

    # this part is just for now: 
    experiment_palette = {}
    experiment_palette['CAB'] = '#F5867F'
    experiment_palette['CAH'] = '#AB4642'
    experiment_palette['CAT'] = '#430300'
    experiment_palette['PU'] = '#F98F34'
    experiment_palette['GP'] = '#6BBC46'
    experiment_palette['GPe'] = '#007600'
    experiment_palette['MGM1'] = '#FF2600'
    experiment_palette['NAC'] = '#0C4E9B'
    experiment_palette['STH'] = '#6B98C4'
    experiment_palette['SUBTH'] = '#6B98C4'

    donor_palette = {
        'UWA7648': '#D87C79',
        'UCI4723': '#7A4300',
        'UCI2424': '#D7A800',
        'UCI5224': '#AB4CAA'
    }

    replicate_palette = {
        "ucsd" : '#039BE5',
        "salk" : '#FFD54F'
    }
    if "experiment_palette" not in adata.uns:
        adata.uns["experiment_palette"] = experiment_palette
    if "donor_palette" not in adata.uns:
        adata.uns["donor_palette"] = donor_palette
    if "replicate_palette" not in adata.uns:
        adata.uns["replicate_palette"] = replicate_palette

    if "brain_region_colors" not in adata.uns:
        cols = []
        for exp_order in adata.obs['brain_region'].cat.categories: 
            cols.append(experiment_palette[exp_order])
        adata.uns['brain_region_colors'] = cols
    if "donor_colors" not in adata.uns:
        cols = []
        for exp_order in adata.obs['donor'].cat.categories: 
            cols.append(donor_palette[exp_order])
        adata.uns['donor_colors'] = cols
    if "replicate_colors" not in adata.uns:
        cols = []
        for exp_order in adata.obs['replicate'].cat.categories: 
            cols.append(replicate_palette[exp_order])
        adata.uns['replicate_colors'] = cols



    # plot region / donor / replicate metadata
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 10), dpi=300)
    plot_categorical(adata, coord_base="umap", cluster_col="brain_region", text_anno=True, coding=True, show=False, ax=axes[0])
    axes[0].set_title(f"Brain Region")
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
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20, 20), dpi=300)
    axes = axes.flatten()
    plot_categorical(adata, coord_base="tsne", cluster_col="brain_region", text_anno=True, coding=True, show=False, ax=axes[0])
    axes[0].set_title(f"Brain Region")
    plot_categorical(adata, coord_base="tsne", cluster_col="donor", text_anno=False, coding=True, show=False, ax=axes[1])
    axes[1].set_title(f"Donor")
    plot_categorical(adata, coord_base="tsne", cluster_col="replicate", text_anno=False, coding=True, show=False, ax=axes[2])
    axes[2].set_title(f"Replicate")
    plot_categorical(adata, coord_base="resolvi_tsne", cluster_col="brain_region", text_anno=True, coding=True, show=False, ax=axes[3])
    axes[3].set_title(f"Brain Region (resolvi_tsne)")
    plot_categorical(adata, coord_base="resolvi_tsne", cluster_col="donor", text_anno=False, coding=True, show=False, ax=axes[4])
    axes[4].set_title(f"Donor (resolvi_tsne)")
    plot_categorical(adata, coord_base="resolvi_tsne", cluster_col="replicate", text_anno=False, coding=True, show=False, ax=axes[5])
    axes[5].set_title(f"Replicate (resolvi_tsne)")
    plot_categorical(adata, coord_base="corr_tsne", cluster_col="brain_region", text_anno=True, coding=True, show=False, ax=axes[6])
    axes[6].set_title(f"Brain Region (corrected expression tsne)")
    plot_categorical(adata, coord_base="corr_tsne", cluster_col="donor", text_anno=False, coding=True, show=False, ax=axes[7])
    axes[7].set_title(f"Donor (corrected expression tsne)")
    plot_categorical(adata, coord_base="corr_tsne", cluster_col="replicate", text_anno=False, coding=True, show=False, ax=axes[8])
    axes[8].set_title(f"Replicate (corrected expression tsne)")
    if save_path:
        fig.savefig(save_path / "tsne_meta.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)    
    
    # plot umap / tsne / spatial coordinates 
    # Plot the UMAP embeddings
    fig, axes = plt.subplots(2, 3, figsize=(20, 20), dpi=300)
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
    plot_continuous(adata, coord_base="tsne", color_by=gene, layer="counts", ax=axes[3], show=False, title=f"tsne - {gene} - raw expression")
    plot_continuous(adata, coord_base="umap", color_by=gene, layer="counts", ax=axes[4], show=False, title=f"umap - {gene} - raw expression")
    plot_continuous(adata, coord_base="spatial", color_by=gene, layer="counts", ax=axes[5], show=False, title=f"spatial - {gene} - raw expression")
    if save_path:
        fig.savefig(save_path / "continuous.png", bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
        plt.close(fig)

    # Plot the spatial coordinates
    plot_spatial_categorical(adata, experiments=experiments, donors=donors, replicates=replicates, 
                             color_key="leiden", combined_legend=True,
                             show=False, output=save_path / "spatial_leiden" if save_path else None)
    plot_spatial_categorical(adata, experiments=experiments, donors=donors, replicates=replicates, 
                             color_key="resolvi_leiden", combined_legend=True,
                             show=False, output=save_path / "spatial_resolvi_leiden" if save_path else None)
    plot_spatial_categorical(adata, experiments=experiments, donors=donors, replicates=replicates, 
                             color_key="corr_leiden", combined_legend=True,
                             show=False, output=save_path / "spatial_corr_leiden" if save_path else None)