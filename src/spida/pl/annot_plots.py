from pathlib import Path
import json
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns
from spida.pl import plot_categorical, plot_continuous, categorical_scatter

def plot_annot(
    adata, 
    hierarchy_list,
    annot_name,
    palette_path : str | Path | dict | None = None, 
    pdf_path : str | Path | None = None,
    plot_filt : bool = True,
    filt_column : str | None = "filt",
    score_ext : str | None = "transfer_score",
    umap_col : str = "base_umap",
    single_clust_plot_thr : int = 15,
):
    """
    Plotting results from an annotation run
    """
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams['axes.facecolor'] = 'white'

    if filt_column is not None and f'{annot_name}_{filt_column}' not in adata.obs.columns:
        raise ValueError(f"Filter column '{filt_column}' not found in adata.obs")
    # TODO: move to utilities as read_palette()
    if palette_path is not None: 
        if isinstance(palette_path, (str, Path)):
            if not Path(palette_path).exists():
                raise FileNotFoundError(f"Palette file not found at {palette_path}")
            if Path(palette_path).suffix == '.json':
                with open(palette_path) as f:
                    palette = json.load(f)
            elif Path(palette_path).suffix in ['.csv', '.tsv']:
                palette_df = pd.read_csv(palette_path, sep='\t' if Path(palette_path).suffix == '.tsv' else ',')
                if 'category' not in palette_df.columns or 'color' not in palette_df.columns:
                    raise ValueError("Palette CSV/TSV must contain 'category' and 'color' columns")
                palette = dict(zip(palette_df['category'], palette_df['color']))
            else:
                raise ValueError("Unsupported palette file format. Use .json, .csv, or .tsv")
        elif isinstance(palette_path, dict):
            palette = palette_path
        else:
            raise ValueError("palette_path must be a string, Path, or dict")

    num_cols = len(hierarchy_list)
    if score_ext is not None: 
        num_rows = 2
    else:
        num_rows = 1

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 15))
        for i, _cc in enumerate(hierarchy_list):
            ax=axes[0][i] if num_rows > 1 else axes[i]
            # If too many types do not annotate text
            text_annot = False if len(adata.obs[f"{annot_name}_{_cc}"].unique()) > 15 else True
            pp = palette[_cc] if palette_path is not None else None # get palette
            # plot
            plot_categorical(adata, cluster_col=f"{annot_name}_{_cc}", coord_base=umap_col, show=False, ax=ax, text_anno=text_annot, coding=True, palette_path=pp)
            ax.set_title(_cc)
            if num_rows > 1:
                ax = axes[1][i]
                plot_continuous(adata, color_by=f"{annot_name}_{_cc}_{score_ext}", coord_base=umap_col, show=False, ax=ax, hue_portion=1, cmap="RdYlBu_r")
                ax.set_title(f"{_cc} transfer score")
        plt.suptitle("Annotation results, UMAP", fontsize=16)
        pdf.savefig(fig)
        plt.close(fig)


        fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 15))
        for i, _cc in enumerate(hierarchy_list):
            ax=axes[0][i] if num_rows > 1 else axes[i]
            pp = palette[_cc] if palette_path is not None else None
            plot_categorical(adata, cluster_col=f"{annot_name}_{_cc}", coord_base='spatial', show=False, ax=ax, text_anno=False, coding=True, palette_path=pp)
            ax.set_title(_cc)
            if num_rows > 1:
                ax = axes[1][i]
                plot_continuous(adata, color_by=f"{annot_name}_{_cc}_{score_ext}", coord_base='spatial', show=False, ax=ax, hue_portion=1, cmap="RdYlBu_r")
                ax.set_title(f"{_cc} transfer score")
        plt.suptitle("Annotation results, Spatial", fontsize=16)
        pdf.savefig(fig)
        plt.close(fig)    


        if score_ext is not None: 
            df_annot = adata.obs[[f"{annot_name}_{_cc}" for _cc in hierarchy_list] + [f"{annot_name}_{_cc}_{score_ext}" for _cc in hierarchy_list]].copy()
            for _cc in hierarchy_list:
                fig, ax = plt.subplots(figsize=(8, 4))
                pp = palette[_cc] if palette_path is not None else None
                sns.violinplot(
                    data=df_annot, y=f"{annot_name}_{_cc}_{score_ext}", x=f"{annot_name}_{_cc}",
                    hue=f"{annot_name}_{_cc}", ax=ax, inner=None, alpha=0.7,
                    palette=pp
                )
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=6)
                ax.set_ylim(0,1.3)
                xaxis_order = [tick.get_text() for tick in ax.get_xticklabels()]
                text = df_annot[f"{annot_name}_{_cc}"].value_counts().reindex(xaxis_order).values.astype(str)
                for i, t in enumerate(text):
                    ax.text(i, ax.get_ylim()[1]-0.01, t, ha='center', va='top', fontsize=6)
                ax.set_title(f"{_cc} transfer score distribution")
                pdf.savefig(fig)
                plt.close(fig)
        
        if plot_filt and filt_column is not None:
            df_annot = adata.obs[[f"{annot_name}_{_cc}" for _cc in hierarchy_list] + [f"{annot_name}_{_cc}_{score_ext}" for _cc in hierarchy_list] + [f"{annot_name}_{filt_column}"]].copy()
            df_annot = df_annot.loc[df_annot[f"{annot_name}_{filt_column}"] == True].copy()
            for _cc in hierarchy_list:
                df_annot[f"{annot_name}_{_cc}"] = df_annot[f"{annot_name}_{_cc}"].astype("category").cat.remove_unused_categories()
                fig, ax = plt.subplots(figsize=(8, 4))
                pp = palette[_cc] if palette_path is not None else None
                sns.violinplot(
                    data=df_annot, y=f"{annot_name}_{_cc}_{score_ext}", x=f"{annot_name}_{_cc}",
                    hue=f"{annot_name}_{_cc}", ax=ax, inner=None, alpha=0.7,
                    palette=pp
                )
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=6)
                # ax.set_ylim(0,1.3)
                xaxis_order = [tick.get_text() for tick in ax.get_xticklabels()]
                text = df_annot[f"{annot_name}_{_cc}"].value_counts().reindex(xaxis_order).values.astype(str)
                for i, t in enumerate(text):
                    ax.text(i, ax.get_ylim()[1]-0.01, t, ha='center', va='top', fontsize=6)
                ax.set_title(f"{_cc} Filtered transfer score distribution", fontsize=10)
                pdf.savefig(fig)
                plt.close(fig)

        for _cc in hierarchy_list:
            if plot_filt and filt_column is not None:
                adata = adata[adata.obs[f"{annot_name}_{filt_column}"] == True].copy()
            if adata.obs[f"{annot_name}_{_cc}"].nunique() <= single_clust_plot_thr:
                ncols = 4
                nrows = int(np.ceil(len(adata.obs[f"{annot_name}_{_cc}"].cat.categories) / ncols))
                pp = palette[_cc] if palette_path is not None else None
                fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
                for ct, ax in zip(adata.obs[f"{annot_name}_{_cc}"].cat.categories, axes.flatten()):
                    adata_sub = adata[adata.obs[f"{annot_name}_{_cc}"] == ct].copy()
                    categorical_scatter(adata, coord_base="spatial", color="lightgrey", ax=ax)
                    plot_categorical(adata_sub, cluster_col=f"{annot_name}_{_cc}", coord_base='spatial', show=False, text_anno=False, ax=ax)
                    ax.set_title(ct)

                for ax in axes.flatten()[len(adata.obs[f"{annot_name}_{_cc}"].cat.categories):]:
                    ax.axis('off')
                # plt.tight_layout()
                plt.suptitle(f"{_cc} spatial distribution; Filtered = {plot_filt and filt_column is not None}", fontsize=16, y=1.02)
                pdf.savefig(fig)
                plt.close(fig)