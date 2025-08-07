# For plotting spatial data (both continuous and categorical)
# Relies on the categorical and scatter plots from _scatterplot
# Author: Amit Klein

import os
import anndata as ad
import numpy as np
from matplotlib.colors import Colormap
import matplotlib.pyplot as plt
from ._utils import (
    fig_make_colorbar, fig_make_legend
)
from ..utilities.sd_utils import _get_obs_or_gene
from .palettes import add_color_scheme
from ._scatteplot import plot_categorical, plot_continuous

def _determine_layout(
    experiments : list | str = [],
    donors : list | str = [], 
    replicates : list | str = [],
): 
    """ Determine the layout of the plots based on the number of replicates, donors, and experiments."""
    
    if isinstance(experiments, str):
        experiments = [experiments]
    if isinstance(donors, str):
        donors = [donors]
    if isinstance(replicates, str): 
        replicates = [replicates]

    num_repl = len(replicates)
    num_donors = len(donors)
    num_experiments = len(experiments)

    # print(f"Number of replicates: {num_repl}, donors: {num_donors}, experiments: {num_experiments}")

    # This nested if statement determines the order of dimensions based on the number of replicates, donors, and experiments.
    if num_repl > 1: 
        num_cols = num_repl
        if num_donors > 1: 
            name_order = {0 : "brain_region", 1 : "donor", 2 : "replicate"}
            data_order = {0 : experiments, 1 : donors, 2 : replicates}
        else: 
            if num_experiments > 1:
                name_order = {0 : "donor", 1 : "brain_region", 2 : "replicate"}
                data_order = {0 : donors, 1 : experiments, 2 : replicates}
                
            else:
                name_order = {0 : "brain_region", 1 : "donor", 2 : "replicate"}
                data_order = {0 : experiments, 1 : donors, 2 : replicates}
    else: 
        if num_donors > 1: 
            name_order = {0 : "replicate", 1 : "brain_region", 2 : "donor"}
            data_order = {0 : replicates, 1 : experiments, 2 : donors}
        else:
            if num_experiments > 1:
                name_order = {0 : "replicate", 1 : "donor", 2 : "brain_region"}
                data_order = {0 : replicates, 1 : donors, 2 : experiments}
            else:
                name_order = {0 : "replicate", 1 : "donor", 2 : "brain_region"}
                data_order = {0 : replicates, 1 : donors, 2 : experiments}

    return name_order, data_order

def plot_spatial_continuous(
    adata : ad.AnnData, 
    experiments : str | list, 
    donors : str | list, 
    replicates : str | list = ['salk'], 
    color_key : str = 'leiden',
    layer : str = None,
    cmap : str | Colormap = 'Wistia', 
    combined_legend : bool = False,
    hspace : float = 0.4,
    wspace : float = 0.2,
    left : float = 0.1,
    right : float = 0.9,
    bottom : float = 0.1,
    top : float = 0.9,
    pmax : int = 99,
    pmin : int = None,
    output : str = None, 
    show : bool = False,
    **kwargs
): 
    """ 
    Plot spatial continuous data with given adata and layout. 

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the spatial data.
    experiments : str | list
        Name(s) of the experiment(s) to plot.
    donors : str | list
        Name(s) of the donor(s) to plot. 
    replicates : str | list, optional
        Name(s) of the replicate(s) to plot (default is ['salk']).
    color_key : str, optional
        Key in adata.obs to color the points by (default is 'leiden').
    layer : str, optional
        Layer in adata to use for the plot (default is None, which is .X).
    cmap : str | mpl.colors.Colormap, optional
        Colormap to use for the plot (default is 'Wistia').
    combined_legend : bool, optional
        Whether to combine legends across subplots (default is False).
    hspace : float, optional
        Height space between subplots (default is 0.4).
    wspace : float, optional
        Width space between subplots (default is 0.2).
    left : float, optional
        Left margin of the figure (default is 0.1).
    right : float, optional
        Right margin of the figure (default is 0.9).
    bottom : float, optional
        Bottom margin of the figure (default is 0.1).
    top : float, optional
        Top margin of the figure (default is 0.9).
    pmax : int, optional
        Maximum percentile for color normalization (default is 99).
    pmin : int, optional
        Minimum percentile for color normalization (default is None, which is 100 - pmax).
    **kwargs : dict
        Additional keyword arguments to pass to the plotting function (plot_continuous).
    """
    name_order, data_order = _determine_layout(
        experiments=experiments,
        donors=donors,
        replicates=replicates,
    )

    num_plots, num_rows, num_cols = (len(data_order[0]), len(data_order[1]), len(data_order[2]))
    # print(f"Number of plots: {num_plots}, rows: {num_rows}, cols: {num_cols}")
    # print(f"Order of dimensions: {name_order}")

    if kwargs is None:
        kwargs = {}

    adata, _drop_col = _get_obs_or_gene(adata, color_key, layer)
    if combined_legend: 
        if pmin is None:
            pmin = 100 - pmax
        hue_norm = np.percentile(adata.obs[color_key], [pmin, pmax])
        kwargs['hue_norm'] = hue_norm
        kwargs['colorbar'] = False


    for i in range(num_plots): 
        _level1 = data_order[0][i]

        fig, axes = plt.subplots(ncols=num_cols, nrows=num_rows, figsize=(5*num_cols, 5*num_rows), dpi=300) #, tight_layout=True)
        for j in range(num_rows): 
            _level2 = data_order[1][j]
            for k in range(num_cols): 
                _level3 = data_order[2][k]
                ax = axes[j, k] if num_rows > 1 else axes[k] if num_cols > 1 else axes
                adata_sub = adata[(adata.obs[name_order[0]] == _level1) & (adata.obs[name_order[1]] == _level2) & (adata.obs[name_order[2]] == _level3)]
                if adata_sub.shape[0] == 0:
                    ax.set_visible(False)
                    continue
                plot_continuous(
                    adata_sub,
                    ax=ax,
                    coord_base="spatial",
                    color_by=color_key,
                    cmap=cmap,
                    layer=layer,
                    show=False,
                    **kwargs
                    )
                
                ax.set_title(f"{_level1} - {_level2} - {_level3}")
        if combined_legend: 
            fig_make_colorbar(
                fig=fig,
                vmin=kwargs['hue_norm'][0],
                vmax=kwargs['hue_norm'][1],
                label=color_key,
                cmap=cmap,
                labelsize=10,
                ticklabel_size=10,
                cbar_width=6,
            )
        plt.suptitle(f"{_level1} - {color_key}", fontsize=16)
        plt.subplots_adjust(top = top, bottom = bottom, left = left, right = right, hspace = hspace, wspace = wspace)
        
        # handle output for each plot made
        if output is not None:
            plot_output = f"{output}_{_level1}_{color_key}.png" 
            plt.savefig(os.path.expanduser(plot_output),bbox_inches='tight',dpi=300)
        if show:
            plt.show()
    
    # wrapping up
    if _drop_col: 
        adata.obs.drop(columns=[color_key], inplace=True)


def plot_spatial_categorical(
    adata : ad.AnnData, 
    experiments : str | list, 
    donors : str | list, 
    replicates : str | list = ['salk'], 
    color_key : str = "leiden",
    coord_base : str = "spatial",
    combined_legend : bool = False,
    output : str = None, 
    show : bool = False,
    **kwargs
): 
    """
    Plot spatial categorical data for given experiments, donors, and replicates.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the spatial data.
    experiments : str | list
        Name(s) of the experiment(s) to plot.
    donors : str | list
        Name(s) of the donor(s) to plot.
    replicates : str | list, optional
        Name(s) of the replicate(s) to plot (default is ['salk']).
    color_key : str, optional
        Key in adata.obs to color the points by (default is 'leiden').
    coord_base : str, optional
        Coordinate base for the plot (default is 'spatial').
    combined_legend : bool, optional
        Whether to combine legends across subplots (default is False).
    **kwargs : dict
        Additional keyword arguments to pass to the plotting function (plot_categorical).
    """
    if color_key not in adata.obs.columns:
        raise ValueError(f"Color key '{color_key}' not found in adata.obs. Please provide a valid key.")
    if f"{color_key}_colors" not in adata.uns:
        _ = add_color_scheme(adata, color_key, palette_key=f"{color_key}_colors")

    name_order, data_order = _determine_layout(
        experiments=experiments,
        donors=donors,
        replicates=replicates,
    )

    num_plots, num_rows, num_cols = (len(data_order[0]), len(data_order[1]), len(data_order[2]))
    print(f"Number of plots: {num_plots}, rows: {num_rows}, cols: {num_cols}")
    print(f"Order of dimensions: {name_order}")

    if kwargs is None:
        kwargs = {}

    if combined_legend: 
        import pandas as pd
        show_legend = False
        ret = True
    else: 
        show_legend = True
        ret = False

    for i in range(num_plots): 
        _level1 = data_order[0][i]

        fig, axes = plt.subplots(ncols=num_cols, nrows=num_rows, figsize=(5*num_cols, 5*num_rows), dpi=300) #, tight_layout=True)
        _data_list = []
        _colors = {}
        for j in range(num_rows): 
            _level2 = data_order[1][j]
            for k in range(num_cols): 
                _level3 = data_order[2][k]
                ax = axes[j, k] if num_rows > 1 else axes[k] if num_cols > 1 else axes
                adata_sub = adata[(adata.obs[name_order[0]] == _level1) & (adata.obs[name_order[1]] == _level2) & (adata.obs[name_order[2]] == _level3)]
                if adata_sub.shape[0] == 0:
                    ax.set_visible(False)
                    continue
                _data, _colors = plot_categorical(
                    adata_sub,
                    ax=ax,
                    coord_base=coord_base,
                    cluster_col=color_key,
                    show=False,
                    ret=ret,
                    show_legend=show_legend,
                    **kwargs
                    )
                if combined_legend:
                    _data_list.append(_data)
                    _colors.update(_colors)
                ax.set_title(f"{_level1} - {_level2} - {_level3}")

        if combined_legend: 
            # TODO: a function to determine labelsize! 
            _data = pd.concat(_data_list, ignore_index=True)
            fig_make_legend(fig, _data=_data, hue=color_key, palette_dict=_colors, labelsize=10)
        plt.suptitle(f"{_level1} - {color_key}", fontsize=16)
        plt.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.1, right = 0.9, hspace = 0.4, wspace = 0.2)
        
        if output is not None:
            plot_output = f"{output}_{_level1}_{color_key}.png" 
            plt.savefig(os.path.expanduser(plot_output),bbox_inches='tight',dpi=300)
        if show:
            plt.show()