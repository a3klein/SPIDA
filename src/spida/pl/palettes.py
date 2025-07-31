# From scanpy.plotting.palettes.py

"""Color palettes in addition to matplotlib's palettes."""
from __future__ import annotations

from typing import TYPE_CHECKING

from matplotlib import cm, colors

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

# Colorblindness adjusted vega_10
# See https://github.com/scverse/scanpy/issues/387
vega_10 = list(map(colors.to_hex, cm.tab10.colors))
vega_10_scanpy = vega_10.copy()
vega_10_scanpy[2] = "#279e68"  # green
vega_10_scanpy[4] = "#aa40fc"  # purple
vega_10_scanpy[8] = "#b5bd61"  # kakhi

default_10 = vega_10_scanpy

# default matplotlib 2.0 palette
# see 'category20' on https://github.com/vega/vega/wiki/Scales#scale-range-literals
vega_20 = list(map(colors.to_hex, cm.tab20.colors))

# reorderd, some removed, some added
vega_20_scanpy = [
    # dark without grey:
    *vega_20[0:14:2],
    *vega_20[16::2],
    # light without grey:
    *vega_20[1:15:2],
    *vega_20[17::2],
    # manual additions:
    "#ad494a",
    "#8c6d31",
]
vega_20_scanpy[2] = vega_10_scanpy[2]
vega_20_scanpy[4] = vega_10_scanpy[4]
vega_20_scanpy[7] = vega_10_scanpy[8]  # kakhi shifted by missing grey
# TODO: also replace pale colors if necessary

default_20 = vega_20_scanpy

# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference https://research.wu.ac.at/en/publications/escaping-rgbland-selecting-colors-for-statistical-graphics-26
zeileis_28 = [
    "#023fa5",
    "#7d87b9",
    "#bec1d4",
    "#d6bcc0",
    "#bb7784",
    "#8e063b",
    "#4a6fe3",
    "#8595e1",
    "#b5bbe3",
    "#e6afb9",
    "#e07b91",
    "#d33f6a",
    "#11c638",
    "#8dd593",
    "#c6dec7",
    "#ead3c6",
    "#f0b98d",
    "#ef9708",
    "#0fcfc0",
    "#9cded6",
    "#d5eae7",
    "#f3e1eb",
    "#f6c4e1",
    "#f79cd4",
    # these last ones were added:
    "#7f7f7f",
    "#c7c7c7",
    "#1CE6FF",
    "#336600",
]

default_28 = zeileis_28

# from https://godsnotwheregodsnot.blogspot.com/2012/09/color-distribution-methodology.html
godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#6A3A4C",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
]

default_102 = godsnot_102


### COLOR SCHEME FROM COPILOT
def add_color_scheme(
    adata,
    column: str,
    palette_key: str | None = None,
    palette: str | list[str] | None = None,
    na_color: str = '#CCCCCC'
):
    """
    Add a color palette to AnnData.uns based on an existing categorical column.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object
    column : str
        Name of the column in adata.obs to create colors for
    palette_key : str, optional
        Key name for the palette in adata.uns. If None, uses f"{column}_colors"
    palette : str or list, optional
        Color palette to use. Can be:
        - matplotlib/seaborn palette name (e.g., 'tab10', 'Set1', 'viridis')
        - list of hex colors
        - None (uses default tab10 for categorical, viridis for continuous)
    na_color : str, default '#CCCCCC'
        Color to use for NaN/missing values
        
    Returns
    -------
    None
        Modifies adata.uns in place by adding the color palette
        
    Examples
    --------
    >>> # For categorical data
    >>> add_color_scheme(adata, 'cell_type')
    >>> add_color_scheme(adata, 'cluster', palette='Set1')
    
    >>> # For continuous data (stores colormap name)
    >>> add_color_scheme(adata, 'n_genes', palette='viridis')
    
    >>> # Custom colors
    >>> custom_colors = ['#FF0000', '#00FF00', '#0000FF']
    >>> add_color_scheme(adata, 'condition', palette=custom_colors)
    """
    import pandas as pd
    import anndata as ad
    from matplotlib.colors import rgb2hex
    import seaborn as sns
    
    if column not in adata.obs.columns:
        raise ValueError(f"Column '{column}' not found in adata.obs")
    
    if palette_key is None:
        palette_key = f"{column}_colors"
    
    # Initialize adata.uns if it doesn't exist
    if adata.uns is None:
        adata.uns = {}
    
    # Get the data
    data = adata.obs[column].copy()
    
    # Handle missing values
    has_na = data.isna().any()
    
    # Determine if data is categorical or continuous
    if isinstance(data.dtype, pd.CategoricalDtype) or pd.api.types.is_object_dtype(data):
        # Categorical data
        if isinstance(data.dtype, pd.CategoricalDtype):
            # Use category order if available
            unique_vals = data.cat.categories.tolist()
        else:
            unique_vals = data.dropna().unique()
        
        n_categories = len(unique_vals)
        
        if palette is None:
            # Use tab10 for small numbers, otherwise generate colors
            if n_categories <= 10:
                palette = default_10
            elif n_categories <= 20:
                palette = default_20
            elif n_categories <= 28:
                palette = default_28
            else:
                palette = default_102
            colors = sns.color_palette(palette, n_categories)
        elif isinstance(palette, str):
            colors = sns.color_palette(palette, n_categories)
        elif isinstance(palette, list):
            if len(palette) < n_categories:
                # Extend palette if too short
                colors = (palette * ((n_categories // len(palette)) + 1))[:n_categories]
            else:
                colors = palette[:n_categories]
        else:
            raise ValueError("palette must be a string or list of colors")
        
        # Convert colors to hex
        hex_colors = [rgb2hex(c) if not isinstance(c, str) else c for c in colors]
        
        # Store color palette in adata.uns following scanpy convention
        # This stores colors in the same order as the categories
        adata.uns[palette_key] = hex_colors
                
        # Return the color mapping for convenience
        return dict(zip(unique_vals, hex_colors))
        
    else:
        # Continuous data - store colormap name or custom colormap
        if palette is None:
            palette = 'viridis'
        
        if isinstance(palette, str):
            # Store colormap name for continuous data
            adata.uns[palette_key] = palette
            # print(f"Added colormap '{palette}' as '{palette_key}' to adata.uns")
        else:
            # Store custom colormap as list of colors
            adata.uns[palette_key] = palette
            # print(f"Added custom colormap as '{palette_key}' to adata.uns")
        return palette