import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import spatialdata as sd
import spatialdata_plot as sdp

from spida.pl import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

from spida.pl.transcript_plots import _plot_rank_genes_dotplot_matplotlib
from spida.utilities.boxes import generate_multiple_random_boxes, visualize_boxes
from spida.pl.P_plots import plot_feature_distribution, plot_filt_cells



def plot_load_images(sdata, image_key, out_dir):
    image_channels = sd.models.get_channel_names(sdata[image_key])
    image_scale_keys = list(sdata[image_key].keys())
    cs = 'global'
        
    max_int = (
        sdata[image_key][image_scale_keys[-1]]["image"]
        .max(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )
    min_int = (
        sdata[image_key][image_scale_keys[-1]]["image"]
        .min(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )


    for i, channel in enumerate(image_channels):
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
        norm = Normalize(vmin=min_int[channel], vmax=max_int[channel] * 0.5)
        sdata.pl.render_images(
            image_key, channel=channel, cmap="grey", norm=norm
        ).pl.show(ax=ax, title=channel, coordinate_systems=cs, colorbar=False)
        plt.savefig(os.path.join(out_dir, f"load_{channel}.svg"))
        plt.close()

def plot_decon_images(sdata, decon_image_key, out_dir):
    image_channels = sd.models.get_channel_names(sdata[decon_image_key])
    image_scale_keys = list(sdata[decon_image_key].keys())
    cs = 'global'
        
    max_int = (
        sdata[decon_image_key][image_scale_keys[-1]]["image"]
        .max(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )
    min_int = (
        sdata[decon_image_key][image_scale_keys[-1]]["image"]
        .min(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )


    for i, channel in enumerate(image_channels):
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
        norm = Normalize(vmin=min_int[channel], vmax=max_int[channel] * 0.5)
        sdata.pl.render_images(
            decon_image_key, channel=channel, cmap="grey", norm=norm
        ).pl.show(ax=ax, title=channel, coordinate_systems=cs, colorbar=False)
        plt.savefig(os.path.join(out_dir, f"decon_{channel}.svg"))
        plt.close()


def plot_tz_qc(grid, out_dir, cmap="RdYlGn_r"): 
    # panel 1: tz_count
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    grid.plot(ax=ax, column="tz_count", cmap=cmap, legend=True)
    ax.set_title("total transcripts")
    ax.axis("off")
    plt.savefig(os.path.join(out_dir, f"tz_qc_raw_count.svg"))
    plt.close()

    # panel 2: density
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    grid['log_tz_count'] = np.log1p(grid['tz_count'])
    grid.plot(ax=ax, column="log_tz_count", cmap=cmap, legend=True)
    ax.set_title("log(total transcripts+1)")
    ax.axis("off")
    plt.savefig(os.path.join(out_dir, f"tz_qc_log_count.svg"))
    plt.close()

    # panel 3: n_genes
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    grid.plot(ax=ax, column="n_genes", cmap=cmap, legend=True)
    ax.set_title("n_genes")
    ax.axis("off")
    plt.savefig(os.path.join(out_dir, f"tz_qc_n_genes.svg"))
    plt.close()

    # panel 4: QC pass/fail over all hexes
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    if "filtered" in grid.columns:
        qc_val = (~grid["filtered"]).astype(int)  # 1=pass, 0=fail
        qc_cmap = ListedColormap(["#d9d9d9", "#1b9e77"])  # fail, pass
        grid.assign(qc_pass=qc_val).plot(
            ax=ax, column="qc_pass", cmap=qc_cmap, legend=False
        )
        handles = [
            mpatches.Patch(color="#1b9e77", label="pass"),
            mpatches.Patch(color="#d9d9d9", label="fail"),
        ]
        ax.legend(handles=handles, title="QC", loc="upper left", bbox_to_anchor=(1.02, 1))
        ax.set_title("QC pass/fail")
    else:
        ax.set_title("QC pass/fail (missing)")
    ax.axis("off")
    plt.savefig(os.path.join(out_dir, f"tz_qc_pass_fail.svg"))
    plt.close()


def plot_tz_hex_qc(adata_hex, out_dir):
    grid = adata_hex.obs.copy()
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    plot_categorical(
        adata_hex, cluster_col="leiden", coord_base="umap", show=False, text_anno=True, coding=True, ax=ax
    )
    ax.set_title("UMAP - leiden")
    plt.savefig(os.path.join(out_dir, f"tz_qc_hex_umap_leiden.svg"))
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    plot_continuous(adata_hex, color_by="tz_count", coord_base="umap", show=False, ax=ax)
    ax.set_title("UMAP - tz_count")
    plt.savefig(os.path.join(out_dir, f"tz_qc_hex_umap_tz_count.svg"))
    plt.close()

    color_dict = {}
    for _l, _c in zip(adata_hex.obs['leiden'].cat.categories, adata_hex.uns['leiden_colors']):
        color_dict[_l] = _c
    grid['leiden_color'] = grid['leiden'].map(color_dict)

    grid_filt = grid[~grid["filtered"]] if "filtered" in grid.columns else grid
    if "leiden" in adata_hex.obs.columns:
        grid_filt = grid_filt.copy()
        grid_filt["leiden"] = adata_hex.obs["leiden"]

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
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
    plt.savefig(os.path.join(out_dir, f"tz_qc_hex_spatial_leiden.svg"))
    plt.close() 

    n_genes = 5
    fig_dot, ax_dot = plt.subplots(1, 1, figsize=(20, 4), dpi=200)
    _plot_rank_genes_dotplot_matplotlib(adata_hex, n_genes=n_genes, cmap="Reds", ax=ax_dot)
    plt.savefig(os.path.join(out_dir, f"tz_qc_hex_umap_leiden_dotplot.svg"))
    plt.close()


def plot_segload(sdata, image_key, shapes_key, segmentation_key, out_dir):
    cs='global'
    # Subsetting for faster cropping
    sub_sdata = sdata.subset([image_key, shapes_key])

    # Generating image metadata
    image_scale_keys = list(sub_sdata[image_key].keys())
    max_int = (
        sub_sdata[image_key][image_scale_keys[-1]]["image"]
        .max(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )
    min_int = (
        sub_sdata[image_key][image_scale_keys[-1]]["image"]
        .min(["x", "y"])
        .compute()
        .to_dataframe()
        .to_dict()["image"]
    )
    norm = Normalize(vmin=min_int["DAPI"], vmax=max_int["DAPI"] * 0.5)

    # generating boxes
    c, frame_height, frame_width = sub_sdata[image_key][image_scale_keys[0]][
        "image"
    ].shape
    boxes = generate_multiple_random_boxes(
        frame_width, frame_height, 1500, 1500, num_boxes=20, avoid_overlap=True
    )

    keep_boxes = []
    for b_cont, _b in enumerate(boxes):
        xmin, ymin, width, height = _b
        xmax = xmin + width
        ymax = ymin + height

        try:
            cropped_sdata = sub_sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=[xmin, ymin],
                max_coordinate=[xmax, ymax],
                target_coordinate_system=cs,
            )
        except Exception as e:
            print(f"Error cropping data: {e}")
            continue
        if shapes_key not in cropped_sdata:
            print(f"Shapes key {shapes_key} not found in cropped_sdata")
            continue

        keep_boxes.append(_b)
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
        pp = cropped_sdata.pl.render_images(
            image_key, channel="DAPI", cmap="grey", norm=norm
        )
        pp = pp.pl.render_shapes(
            shapes_key, fill_alpha=0, cmap="jet", outline_alpha=1, outline_color="red"
        )
        pp.pl.show(ax=ax, title="Shapes", coordinate_systems=cs, colorbar=False)
        ax.axis("off")
        ax.set_title("Shapes - " + str(b_cont + 1))
        plt.savefig(os.path.join(out_dir, f"segload_{segmentation_key}_shapes_box_{b_cont+1}.svg"))
        plt.close()

    # plotting
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)

    pp = sub_sdata.pl.render_images(
        image_key, channel="DAPI", cmap="grey", norm=norm
    ).pl.show(ax=ax, title="DAPI", coordinate_systems=cs, colorbar=False)
    visualize_boxes(
        frame_width=frame_width,
        frame_height=frame_height,
        boxes=boxes,
        ax=ax,
        show_plot=False,
    )
    plt.savefig(os.path.join(out_dir, f"segload_{segmentation_key}_boxes_on_dapi.svg"))
    plt.close()


def plot_seg_qc(adata, segmentation_key, out_dir):
    plt.rcParams["axes.facecolor"] = "white"

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_feature_distribution(
        adata.obs,
        feature="volume",
        feature_alias="Cell Volume",
        min_val=adata.uns["cutoffs"]["volume_min"],
        max_val=adata.uns["cutoffs"]["volume_max"],
        ax=ax,
    )
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_volume_distribution.svg"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_feature_distribution(
        adata.obs,
        feature="nCount_RNA",
        feature_alias="RNA Count",
        min_val=adata.uns["cutoffs"]["n_count_min"],
        max_val=adata.uns["cutoffs"]["n_count_max"],
        ax=ax,
    )
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_rna_count_distribution.svg"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_feature_distribution(
        adata.obs,
        feature="nFeature_RNA",
        feature_alias="Num Unique RNAs",
        min_val=adata.uns["cutoffs"]["n_gene_min"],
        max_val=adata.uns["cutoffs"]["n_gene_max"],
        ax=ax,
    )
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_n_gene_distribution.svg"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_feature_distribution(
        adata.obs,
        feature="nBlank",
        feature_alias="Blank Count",
        min_val=adata.uns["cutoffs"]["n_blank_min"],
        max_val=adata.uns["cutoffs"]["n_blank_max"],
        ax=ax,
    )
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_blank_count_distribution.svg"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_feature_distribution(
        adata.obs[adata.obs["pass_qc_pre"]] if "pass_qc_pre" in adata.obs.columns else adata.obs,
        feature="nCount_RNA_per_Volume",
        feature_alias="Transcripts Per Cell Size",
        min_val=adata.uns["cutoffs"]["n_count_per_volume_min"],
        max_val=adata.uns["cutoffs"]["n_count_per_volume_max"],
        ax=ax,
    )
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_n_count_per_volume_distribution.svg"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    plot_filt_cells(adata.obs, ax=ax, title="Filtered Cells")
    plt.savefig(os.path.join(out_dir, f"seg_qc_{segmentation_key}_filtered_cells.svg"))
    plt.close()

def plot_cell_cluster(adata, segmentation_key, out_dir): 
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_continuous(adata, coord_base="base_umap", color_by="volume", ax=ax, show=False, title="umap - Cell Volume")
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_umap_volume.svg"))
        plt.close()
        
        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_continuous(adata, coord_base="base_umap", color_by="nCount_RNA", ax=ax, show=False, title="umap - Number of Transcripts")
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_umap_nCount_RNA.svg"))
        plt.close()

        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_continuous(adata, coord_base="base_umap", color_by="nFeature_RNA", ax=ax, show=False, title="umap - Number of Genes")
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_umap_nFeature_RNA.svg"))
        plt.close()
        
        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_continuous(adata, coord_base="base_umap", color_by="nCount_RNA_per_Volume", ax=ax, show=False, title="umap - Number of Transcripts per Volume")
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_umap_nCount_RNA_per_Volume.svg"))
        plt.close()
        
        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_categorical(adata, coord_base="base_umap", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=ax)
        ax.set_title(f"umap - Leiden Clusters")    
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_umap_leiden.svg"))
        plt.close()

        fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
        plot_categorical(adata, coord_base="spatial", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=ax)
        ax.set_title(f"spatial - Leiden Clusters")    
        plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_spatial_leiden.svg"))
        plt.close()

        for _cat in adata.obs['base_leiden'].cat.categories:
            fig, ax = plt.subplots(1 ,1, figsize=(5, 5), dpi=300)
            categorical_scatter(adata, coord_base='spatial', color='lightgrey', ax=ax)
            adata_sub = adata[adata.obs['base_leiden'] == _cat].copy()
            plot_categorical(adata_sub, coord_base="spatial", cluster_col="base_leiden", text_anno=False, coding=True, show=False, ax=ax)
            ax.set_title(f"spatial - Leiden Cluster {_cat}")
            plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_spatial_leiden_cluster_{_cat}.svg"))
            plt.close()

def plot_seg_clust_dotplot(adata, segmentation_key, out_dir):
    # This should not be here no more. 
    # import scanpy as sc    
    # sc.tl.rank_genes_groups(adata, groupby="base_leiden", method="t-test_overestim_var")
    n_genes = 5
    fig_dot, ax_dot = plt.subplots(1, 1, figsize=(20, 4), dpi=200, constrained_layout=True)
    _plot_rank_genes_dotplot_matplotlib(adata, cluster_col="base_leiden", n_genes=n_genes, cmap="Reds", ax=ax_dot)
    plt.savefig(os.path.join(out_dir, f"seg_clust_{segmentation_key}_leiden_dotplot.svg"))
    plt.close()