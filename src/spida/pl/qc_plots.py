### Author: Amit Klein; Date: 07/29/2025; email: a3klein@ucsd.edu

# This file contains code for generating the QC plots for the aggregated experiments. 

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
donor_palette = {'UWA7648': '#D87C79','UCI4723': '#7A4300','UCI2424': '#D7A800','UCI5224': '#AB4CAA'}



# def plot_qc_dataset(): 
def _plot_violin_QC(adata, metric, title:str|None=None, cutoffs_dict = None, top_cut_name:str=None, bottom_cut_name:str=None): 
    """Plots a violin plot for a given metric in the adata object."""
    plot_cutoffs = True if cutoffs_dict is not None else False

    plots = adata.obs[['experiment', 'donor', 'replicate']].drop_duplicates().sort_values(by=['experiment', 'donor']).values
    n_plots = len(plots)
    toplot = adata.obs
    order = plots[:, 0] + "_" + plots[:, 1] + "_" + plots[:, 2]
    toplot['indexer'] = toplot['experiment'].astype(str) + "_" + toplot['donor'].astype(str) + "_" + toplot['replicate'].astype(str)
    toplot['indexer'] = toplot['indexer'].astype('category')

    fig, ax = plt.subplots(figsize=(30, 5), nrows = 1, ncols = 1)

    sns.violinplot(data = toplot, x='indexer', y=metric, order=order, hue="donor",
                   ax=ax, linewidth=0, palette=donor_palette
    )
    ax.set_xlabel("")
    ax.set_ylabel(metric)

    xticks = ax.get_xticks()
    ticks = []
    tick_labels = []
    pre_exp = plots[0, 0]
    pre_tick = xticks[0]
    for tick_mark, (_exp, _donor, _replicates) in zip(xticks, plots): 
        if plot_cutoffs: 
            cutoffs = cutoffs_dict[f"{_exp}_{_donor}_tables"]['cutoffs']
            top_boundary = cutoffs[top_cut_name]
            bottom_boundary = cutoffs[bottom_cut_name]
            xmin = (tick_mark) / len(xticks)
            xmax = (tick_mark + 1) / len(xticks)
            ax.axhline(top_boundary, xmin=xmin, xmax=xmax, color='red', linestyle='--', linewidth=2)
            ax.axhline(bottom_boundary, xmin=xmin, xmax=xmax, color='red', linestyle='--', linewidth=2)
        if (pre_exp != _exp): 
            ax.axvline(tick_mark - 0.5, color='black', linestyle='--', linewidth=2)
            ticks.append( (pre_tick + tick_mark - 1) / 2 )
            tick_labels.append(pre_exp)
            pre_exp = _exp
            pre_tick = tick_mark

    ticks.append( (pre_tick + tick_mark) / 2 )
    tick_labels.append(pre_exp)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)

    ax.legend(title="Donor", bbox_to_anchor=(1, 1), loc='upper left', fontsize=12)

    if title is not None: 
        ax.set_title(title, fontsize=20)

    plt.show()

def plot_violin_QC(
        adata, qc_metric:str|list[str],
        title:str|list[str],
        cutoffs_dict=None,
        top_cut_names:str|list[str]=None,
        bottom_cut_names:str|list[str]=None
    ):
    """Plots a violin plot for a given QC metric."""
    mpl.style.use('default')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    # mpl.rcParams['font.family'] = 'sans-serif'
    # mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['figure.dpi'] = 80
    mpl.rcParams['savefig.dpi']=300
    
    if isinstance(qc_metric, str):
        qc_metric = [qc_metric]

    if isinstance(title, str):
        title = [title] * len(qc_metric)

    for _metric, _title, _top_cut_name, _bottom_cut_name in zip(qc_metric, title, top_cut_names, bottom_cut_names):
        _plot_violin_QC(adata, _metric, _title, cutoffs_dict, _top_cut_name, _bottom_cut_name)