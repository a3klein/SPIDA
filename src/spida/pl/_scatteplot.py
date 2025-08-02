# Functions for plotting single cell data in python. 
# Functions come from ALLCools, and Pym3c by Hanqing Liu and Wubin Ding 
# Adapted by Amit Klein

import os
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, Colormap
import seaborn as sns

from ._utils import (
    _extract_coords, _density_based_sample, _auto_size, 
    zoom_min_max, zoom_ax, _take_data_series, level_one_palette,
    _calculate_luminance, tight_hue_range, _text_anno_scatter,
    density_contour, despine, _make_tiny_axis_label, plot_text_legend,
    plot_color_dict_legend, plot_marker_legend, plot_cmap_legend, get_cmap
)
from ..utilities.sd_utils import _get_obs_or_gene

from .palettes import add_color_scheme

def continuous_scatter(
    data : pd.DataFrame | ad.AnnData,
    ax=None,
    coord_base : str ="umap",
    x : str =None,
    y : str =None,
    scatter_kws=None,
    hue : str =None,
    hue_norm : list = None,
    hue_portion : float = 0.95,
    color=None,
    cmap ="viridis",
    colorbar=True,
    size=None,
    size_norm=None,
    size_portion=0.95,
    sizes=None,
    sizebar=True,
    text_anno=None,
    dodge_text=False,
    dodge_kws=None,
    text_kws=None,luminance=0.48,
    text_transform=None,
    axis_format="tiny",
    max_points=50000,
    s : str | float ="auto",
    labelsize=6,
	ticklabel_size=4,
    linewidth=0.5,
    zoomxy=1.05,
    outline=None,
    outline_kws=None,
    outline_pad=2,
    return_fig : bool = False,
    rasterized : str  | bool ="auto",
    cbar_kws=None,
	cbar_width=3,
):
	"""
	Plot scatter on given adata.

	Parameters
	----------
	data : pd.DataFrame | ad.AnnData
		The data to plot
	ax : _type_, optional
		a matplotlib axis if part of a figure, by default None
	coord_base : str, optional
		the coord base to plot ("umap", "tsne", "spatial", ...), by default "umap"
	x : str, optional
		the x coordinate to plot, by default None
	y : str, optional
		the y coordinate to plot, by default None
	scatter_kws : dict, optional
		kwargs passed to sns.scatterplot, by default None
	hue : str, optional
		the column name to use for color encoding, by default None
	hue_norm : list, optional
		the vmin and vmax for the plot colors, by default None
	hue_portion : float, optional
		The cutoffs to use for building the norm object (quantiles) must be between 0-1, by default 0.95
	color : _type_, optional
		_description_, by default None
	cmap : str, optional
		name of colormap or actual colormap, by default "viridis"
	colorbar : bool, optional
		_description_, by default True
	size : _type_, optional
		_description_, by default None
	size_norm : _type_, optional
		_description_, by default None
	size_portion : float, optional
		_description_, by default 0.95
	sizes : _type_, optional
		_description_, by default None
	sizebar : bool, optional
		_description_, by default True
	text_anno : _type_, optional
		_description_, by default None
	dodge_text : bool, optional
		_description_, by default False
	dodge_kws : _type_, optional
		_description_, by default None
	text_kws : _type_, optional
		_description_, by default None
	luminance : float, optional
		_description_, by default 0.48
	text_transform : _type_, optional
		_description_, by default None
	axis_format : str, optional
		_description_, by default "tiny"
	max_points : int, optional
		The maximum number to plot (subsets to this size), by default 50000
	s : str | float, optional
		the size of the dots, if not set it is calculated proportional to the number of points, by default "auto"
	labelsize : int, optional
		_description_, by default 6
	ticklabel_size : int, optional
		_description_, by default 4
	linewidth : float, optional
		_description_, by default 0.5
	zoomxy : float, optional
		_description_, by default 1.05
	outline : _type_, optional
		_description_, by default None
	outline_kws : _type_, optional
		_description_, by default None
	outline_pad : int, optional
		_description_, by default 2
	return_fig : bool, optional
		Whether to return the figure object, by default False
	rasterized : str, optional
		whether to rasterize the resulting image, by default "auto"
	cbar_kws : _type_, optional
		_description_, by default None
	cbar_width : int, optional
		width of colorbar, by default 3 mm

	Returns
	-------
	_type_
		_description_

	Raises
	------
	ValueError
		_description_
	TypeError
		_description_
	"""
	import seaborn as sns
	import copy
	from matplotlib.cm import ScalarMappable
	# init figure if not provided
	if ax is None:
		fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
	else:
		fig = None

	# add coords
	_data, x, y = _extract_coords(data, coord_base, x, y)
	# _data has 2 cols: "x" and "y"

	# down sample plot data if needed.
	if max_points is not None:
		if _data.shape[0] > max_points:
			_data = _density_based_sample(_data, seed=1, size=max_points, coords=["x", "y"])
	n_dots = _data.shape[0]

	# determine rasterized
	if rasterized == "auto":
		if n_dots > 200:
			rasterized = True
		else:
			rasterized = False

	# auto size if user didn't provide one
	if s == "auto":
		s = _auto_size(ax, n_dots)

	# default scatter options
	_scatter_kws = {"linewidth": 0, "s": s, "legend": None, "rasterized": rasterized}
	if color is not None:
		if hue is not None:
			raise ValueError("Only one of color and hue can be provided")
		_scatter_kws["color"] = color
	if scatter_kws is not None:
		_scatter_kws.update(scatter_kws)

	# deal with color
	if hue is not None:
		if isinstance(hue, str):
			_data["hue"] = _take_data_series(data, hue).astype(float)
			colorbar_label = hue
		else:
			_data["hue"] = hue.astype(float)
			colorbar_label = hue.name

		if hue_norm is None:
			# get the smallest range that include "hue_portion" of data
			# hue_norm = tight_hue_range(_data["hue"], hue_portion)
			hue_norm=(_data["hue"].quantile(1-hue_portion),_data["hue"].quantile(hue_portion))
		# cnorm is the normalizer for color
		cnorm = Normalize(vmin=hue_norm[0], vmax=hue_norm[1])
		if isinstance(cmap, str):
			# from here, cmap become colormap object
			cmap = copy.copy(get_cmap(cmap))
			cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5)) # replace with na_color
		else:
			if not isinstance(cmap, ScalarMappable):
				raise TypeError(f"cmap can only be str or ScalarMappable, got {type(cmap)}")
	else:
		hue_norm = None
		cnorm = None
		colorbar_label = ""

	# deal with size
	if size is not None:
		if isinstance(size, str):
			_data["size"] = _take_data_series(data, size).astype(float)
		else:
			_data["size"] = size.astype(float)
		size = "size"

		if size_norm is None:
			# get the smallest range that include "size_portion" of data
			size_norm = tight_hue_range(_data["size"], size_portion)

			# snorm is the normalizer for size
			size_norm = Normalize(vmin=size_norm[0], vmax=size_norm[1])

		# replace s with sizes
		s = _scatter_kws.pop("s")
		if sizes is None:
			sizes = (min(s, 1), s)
	else:
		size_norm = None
		sizes = None

	sns.scatterplot(
		x="x",
		y="y",
		data=_data,
		hue="hue",
		palette=cmap,
		hue_norm=cnorm,
		size=size,
		sizes=sizes,
		size_norm=size_norm,
		ax=ax,
		**_scatter_kws,
	)

	if text_anno is not None:
		if isinstance(text_anno, str):
			_data["text_anno"] = _take_data_series(data, text_anno)
		else:
			_data["text_anno"] = text_anno
		if str(_data["text_anno"].dtype) == "category":
			_data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

		_text_anno_scatter(
			data=_data[["x", "y", "text_anno"]],
			ax=ax,
			x="x",
			y="y",
			dodge_text=dodge_text,
			dodge_kws=dodge_kws,
			text_transform=text_transform,
			anno_col="text_anno",
			text_kws=text_kws,
			luminance=luminance,
		)

	# deal with outline
	if outline:
		if isinstance(outline, str):
			_data["outline"] = _take_data_series(data, outline)
		else:
			_data["outline"] = outline
		_outline_kws = {
			"linewidth": linewidth,
			"palette": None,
			"c": "lightgray",
			"single_contour_pad": outline_pad,
		}
		if outline_kws is not None:
			_outline_kws.update(outline_kws)
		density_contour(ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws)

	# clean axis
	if axis_format == "tiny":
		_make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
	elif (axis_format == "empty") or (axis_format is None):
		despine(ax=ax, left=True, bottom=True)
		ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
	else:
		pass

	return_axes = [ax]

	# make color bar # TODO: replace this with the function from ._utils
	if colorbar and (hue is not None):
		# small ax for colorbar
		# default_cbar_kws=dict(loc="upper left", borderpad=0,width="3%", height="20%") #bbox_to_anchor=(1,1)
		if cbar_kws is None:
			cbar_kws={}
		# for k in default_cbar_kws:
		#     if k not in cbar_kws:
		#         cbar_kws[k]=default_cbar_kws[k]

		mm2inch = 1 / 25.4
		space=0
		legend_width = (
			cbar_width * mm2inch * ax.figure.dpi / ax.figure.get_window_extent().width
		)  # mm to px to fraction
		pad = (space + ax.yaxis.labelpad * 1.2 * ax.figure.dpi / 72) / ax.figure.get_window_extent().width
		# labelpad unit is points
		ax_pos = ax.get_position()
		left = ax_pos.x1 + pad
		ax_legend = ax.figure.add_axes(
			[left, ax_pos.y0 + ax_pos.height * 0.8, legend_width, ax_pos.height * 0.2]
		)  # left, bottom, width, height
		# print("test:",hue_norm)
		# cbar_kws.setdefault('vmin',hue_norm[0])
		# cbar_kws.setdefault('vmax',hue_norm[1])
		cbar_kws['vmin']=hue_norm[0]
		cbar_kws['vmax']=hue_norm[1]
		cbar = plot_cmap_legend(
			ax=ax,
			cax=ax_legend,
			cmap=cmap,
			label=hue,
			kws=cbar_kws.copy(),labelsize=labelsize, 
			linewidth=linewidth,ticklabel_size=ticklabel_size,
			)
		return_axes.append([ax_legend,cbar])

	# make size bar
	if sizebar and (size is not None):
		# TODO plot dot size bar
		pass

	if zoomxy is not None:
		zoom_ax(ax, zoomxy)

	if return_fig:
		return (fig, tuple(return_axes)), _data
	else:
		return

def plot_continuous(
    adata : ad.AnnData | str,
    ax=None,
    coord_base='tsne',
    color_by:str="nCount_RNA",
	cmap:str|Colormap="Wistia",
    hue_portion : float = 0.99,
    layer:str=None,
    figsize=(4, 3.5),
	ncol:int=2,
    output:str=None,
    show:bool=True,
    title:str=None,
    **kwargs
): 
	"""Plotting continuous data from an AnnData object.
	
	Parameters
    ----------
    adata : ad.AnnData | str
		AnnData object or path to an AnnData file.
    ax :
		matplotlib axis to plot on, if None, a new figure and axis will be created.
    coord_base : str
		Coordinate system to use for the plot, default is 'tsne'.
    color_by : str 
		Column name in adata.obs or gene name in adata.var_names to use for coloring the points.
    cmap : str |  matplotlib.colors.Colormap
		Colormap to use for the plot, default is 'Wistia'.
	hue_portion : float
		The portion of the data to use for the color normalization, default is 0.99.
	layer : str
		Name of the layer to use for the plot, default is None (uses adata.X).
    figsize : tuple
		Figure size, default is (4, 3.5).
	ncol : int
		Number of columns in the figure, default is None.
	output : str
		Path to save the output figure, default is None.
    show : bool
		Whether to show the plot, default is True.
    title : str
		Title of the plot, default is None.
    kwargs : dict
        set text_anno=None to plot clustering without text annotations,
        coding=True to plot clustering without code annotations,
        set show_legend=False to remove the legend

    Returns
    -------
	"""
	import os

	if isinstance(adata,str): # getting adata
		adata=ad.read_h5ad(adata,backed='r')

	# check layer
	if layer is not None: # make sure layer exists if it is being used
		if layer not in adata.layers:
			raise ValueError(f"Layer {layer} not found in adata.layers, please check the layer name.")    
		
	adata, _drop_col = _get_obs_or_gene(adata, color_by, layer) # get the column from obs or var
		
	# figure out cmap: 
	if f"{color_by}_cmap" in adata.uns: # check if cmap is in uns
		cmap = adata.uns[f"{color_by}_cmap"]

	if ax is None:
		fig, ax = plt.subplots(figsize=figsize, dpi=300)
	# handle_cmap
	cax = continuous_scatter(
		data=adata,
		ax=ax,
		coord_base=coord_base,
		hue=color_by,
		hue_portion=hue_portion,
		cmap=cmap,
		**kwargs
	)

	if title is not None: # set title
		ax.set_title(title)

	if _drop_col:  # drop the column if it was added for plotting
		adata.obs.drop(columns=[color_by], inplace=True)

	# Handle outputs
	if output is not None:
		plt.savefig(os.path.expanduser(output),bbox_inches='tight',dpi=300)
	if show:
		plt.show()



def categorical_scatter(
    data,
	ax=None,
    coord_base="umap",
    x=None,
    y=None, # coords
    hue=None,
    palette="auto",
    color=None, # color
    text_anno=None,
    text_kws=None,
    luminance=None,
    text_transform=None,
    dodge_text=False,
    dodge_kws=None, # text annotation
    show_legend=False,
    legend_kws=None, # legend
    s="auto",
    size=None,
    sizes=None, # sizes is a dict
    size_norm=None,
    size_portion=0.95, 
    axis_format="tiny",
    max_points=50000,
    labelsize=4,
    linewidth=0.5,
    zoomxy=1.05,
    outline=None,
    outline_pad=3,
    alpha=0.7,
    outline_kws=None,
    scatter_kws=None,
    rasterized="auto",
    coding=False,
	id_marker=True,
    legend_color_text=True,
	rectangle_marker=False,
    marker_fontsize=4,
    marker_pad=0.1,
):
	"""
	This function was copied from ALLCools and made some modifications.
	Plot categorical scatter plot with versatile options.

	Parameters
	----------
	rasterized
		Whether to rasterize the figure.
	return_fig
		Whether to return the figure.
	size_portion
		The portion of the figure to be used for the size norm.
	data
		Dataframe that contains coordinates and categorical variables
	ax
		this function do not generate ax, must provide an ax
	coord_base
		coords name, if provided, will automatically search for x and y
	x
		x coord name
	y
		y coord name
	hue : str
		categorical col name or series for color hue.
	palette : str or dict
		palette for color hue.
	color
		specify single color for all the dots
	text_anno
		categorical col name or series for text annotation.
	text_kws
		kwargs pass to plt.text, see: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
		including bbox, to see parameter for bbox, go to: https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyBboxPatch.html#matplotlib.patches.FancyBboxPatch
		commonly used parameters are: 
		```
		text_kws=dict(fontsize=5,fontweight='black',
					color='black', #color could be a dict, keys are text to be annotated
					bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
								facecolor=(0.8, 0.8, 0.8, 0.2), #facecolor could also be a dict
								alpha=1,linewidth=0.5)
					)
		```
	text_transform
		transform for text annotation.
	dodge_text
		whether to dodge text annotation.
	dodge_kws
		kwargs for dodge text annotation.
	show_legend
		whether to show legend.
	legend_kws
		kwargs for legend.
	s
		single size value of all the dots.
	size
		mappable size of the dots.
	sizes
		mapping size to the sizes value.
	size_norm
		normalize size range for mapping.
	axis_format
		axis format.
	max_points
		maximum number of points to plot.
	labelsize
		label size pass to `ax.text`
	linewidth
		line width pass to `ax.scatter`
	zoomxy
		zoom factor for x and y-axis.
	outline
		categorical col name or series for outline.
	outline_pad
		outline padding.
	outline_kws
		kwargs for outline.
	scatter_kws
		kwargs for scatter.

	Returns
	-------
	if return_fig is True, return the figure and axes.
	else, return None.
	"""
	if ax is None:
		ax = plt.gca()

	# add coords _data has 2 cols: "x" and "y", index are obs_names
	_data, x, y = _extract_coords(data, coord_base, x, y)
	
	# down sample plot data if needed. (MAX_POINTS)
	if max_points is not None:
		if _data.shape[0] > max_points:
			_data = _density_based_sample(_data, seed=1, size=max_points, coords=["x", "y"])
	n_dots = _data.shape[0]

	# determine rasterized
	if rasterized == "auto":
		if n_dots > 200:
			rasterized = True
		else:
			rasterized = False

	# auto size if user didn't provide one
	if s == "auto":
		s = _auto_size(ax, n_dots)

	# default scatter options
	_scatter_kws = {"linewidth": 0, "s": s, "legend": None, "palette": palette, "rasterized": rasterized}
	if color is not None:
		if hue is not None:
			raise ValueError("Only one of color and hue can be provided")
		_scatter_kws["color"] = color
	if scatter_kws is not None:
		_scatter_kws.update(scatter_kws)

	# deal with color
	palette_dict = None
	if hue is not None:
		if isinstance(hue, str):
			_data["hue"] = _take_data_series(data, hue)
		else:
			_data["hue"] = hue.copy()
		_data["hue"] = _data["hue"].astype("category").cat.remove_unused_categories()

		# if the object has get_palette method, use it (AnnotZarr)
		palette = _scatter_kws["palette"]
		# deal with other color palette
		if palette_dict is None:
			if isinstance(palette, str) or isinstance(palette, list):
				palette_dict = level_one_palette(_data["hue"], order=None, palette=palette)
			elif isinstance(palette, dict):
				palette_dict = palette
			else:
				raise TypeError(f"Palette can only be str, list or dict, " f"got {type(palette)}")
		_scatter_kws["palette"] = palette_dict

	# deal with size
	if size is not None:
		if isinstance(size, str):
			_data["size"] = _take_data_series(data, size).astype(float)
		else:
			_data["size"] = size.astype(float)
		size = "size"

		if size_norm is None:
			# get the smallest range that include "size_portion" of data
			size_norm = tight_hue_range(_data["size"], size_portion)

			# snorm is the normalizer for size
			size_norm = Normalize(vmin=size_norm[0], vmax=size_norm[1])

		# discard s from _scatter_kws and use size in sns.scatterplot
		s = _scatter_kws.pop("s")
		if sizes is None:
			sizes = (min(s, 1), s)

	sns.scatterplot(
		x="x",
		y="y",
		data=_data,
		ax=ax,
		hue="hue",
		size=size,
		sizes=sizes,
		size_norm=size_norm,
		**_scatter_kws,
	)

	# deal with text annotation
	code2label=None
	if text_anno is not None:
		# data
		if isinstance(text_anno, str):
			_data["text_anno"] = _take_data_series(data, text_anno)
		else:
			_data["text_anno"] = text_anno.copy()
		if str(_data["text_anno"].dtype) == "category":
			_data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

		# text kws
		text_kws = {} if text_kws is None else text_kws
		default_text_kws = dict(
			color='white',  # color for the text, could be a dict, keys are text to be annotated
			fontweight="bold", #fontsize=labelsize,
			bbox=dict(facecolor=palette_dict, # if None, use default color
				boxstyle='round', #ellipse, round
				edgecolor='white', fill=True, linewidth=linewidth, alpha=alpha))
		# coding & id_marker
		text_anno='text_anno'
		if not coding is None and coding!=False:
			if coding == True:
				_data['code'] = _data['hue'].cat.codes #int
			else:
				assert isinstance(coding,str)
				_data["code"] = _take_data_series(data, coding)
				_data=_data.loc[_data['code'].notna()]
				_data["code"]=_data["code"].astype(int)
			_data["code"] = _data["code"].astype("category").cat.remove_unused_categories()
			text_anno='code'
			_data['color'] = _data['hue'].map(palette_dict)
			code2label=_data.loc[:,['code','hue']].drop_duplicates().set_index('code').hue.to_dict()
			_data['code']=_data['code'].astype(str)
			code_colors=_data.loc[:,['code','color']].drop_duplicates().set_index('code').color.to_dict()
			default_text_kws['bbox']['facecolor']=code_colors # background colors for text annotation
			default_text_kws['bbox']['boxstyle'] = 'circle'
		for k in default_text_kws:
			if k !='bbox':
				text_kws.setdefault(k, default_text_kws[k])
			else:
				if 'bbox' not in text_kws:
					text_kws['bbox']={}
				for k1 in default_text_kws['bbox']:
					text_kws['bbox'].setdefault(k1, default_text_kws['bbox'][k1])

		_text_anno_scatter(
			data=_data[["x", "y", text_anno]],
			ax=ax,
			x="x",
			y="y",
			dodge_text=dodge_text,
			dodge_kws=dodge_kws,
			text_transform=text_transform,
			anno_col=text_anno,
			text_kws=text_kws,
			luminance=luminance,
		)

	# deal with outline
	if not outline is None:
		if isinstance(outline, str):
			_data["outline"] = _take_data_series(data, outline)
		else:
			_data["outline"] = outline.copy()
		_outline_kws = {
			"linewidth": linewidth,
			"palette": None,
			"c": "lightgray",
			"single_contour_pad": outline_pad,
		}
		if outline_kws is not None:
			_outline_kws.update(outline_kws)
		density_contour(ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws)

	# clean axis
	if axis_format == "tiny":
		_make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
	elif (axis_format == "empty") or (axis_format is None):
		despine(ax=ax, left=True, bottom=True)
		ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
	else:
		pass

	# deal with legend
	# TODO: replace this with ax_add_legend() from ._utils
	if show_legend and (hue is not None):
		n_hue = len(palette_dict)
		ncol=1 if n_hue <= 40 else 2 if n_hue <= 100 else 3
		if legend_kws is None:
			legend_kws = {}
		default_lgd_kws = dict(
			ncol=ncol,fontsize=labelsize,
			bbox_to_anchor=(1, 1),loc="upper left",
			# borderpad=0.4, # pad between marker (text) and border
			# labelspacing=0.2, #The vertical space between the legend entries, in font-size units
			# handleheight=0.5, #The height of the legend handles, in font-size units.
			# handletextpad=0.2, # The pad between the legend handle (marker) and text, in font-size units.
			# borderaxespad=0.3, # The pad between the Axes and legend border, in font-size units
			# columnspacing=0.2, #The spacing between columns, in font-size units
			markersize=labelsize #legend_kws["fontsize"],
		)
		for k in default_lgd_kws:
			legend_kws.setdefault(k, default_lgd_kws[k])

		exist_hues = _data["hue"].unique()
		color_dict={hue_name: color for hue_name, color in palette_dict.items() if hue_name in exist_hues}
		
		if not code2label is None and id_marker:
			boxstyle='Circle' if not rectangle_marker else 'Round'
			plot_text_legend(color_dict, code2label, ax, title=hue, 
					color_text=legend_color_text, boxstyle=boxstyle,marker_pad=marker_pad,
					legend_kws=legend_kws,marker_fontsize=marker_fontsize,
					alpha=alpha,luminance=luminance)
		else:
			if rectangle_marker:
				## plot Patch legend (rectangle marker)
				plot_color_dict_legend(
					D=color_dict, ax=ax, title=hue, color_text=legend_color_text, 
					kws=legend_kws,luminance=luminance
				)
			else:
				# plot marker legend (for example, circle marker)
				plot_marker_legend(
					color_dict=color_dict, ax=ax, title=hue, color_text=legend_color_text, 
					marker='o',kws=legend_kws,luminance=luminance
				)

	if zoomxy is not None:
		zoom_ax(ax, zoomxy)

	return _data

def plot_categorical(
    adata : ad.AnnData | str, 
    ax=None,
    coord_base='tsne',
    cluster_col='MajorType',
    palette_path=None,
    axis_format='tiny',
    alpha=0.7,
    coding=True,
    id_marker=True,
    output=None,
    show=True,
    figsize=(4, 3.5),
    sheet_name=None,
    ncol=None,
    fontsize=5,
    legend_fontsize=5,
    legend_kws=None,
    legend_title_fontsize=5,
    marker_fontsize=4,
    marker_pad=0.1,
    linewidth=0.5,
    ret=False,
    text_anno:bool = False,
    text_kws=None,
    **kwargs
):
    """
    Plot cluster.

    Parameters
    ----------
    adata_path :
    ax :
    coord_base :
    cluster_col :
    palette_path :
    coding :
    output :
    show :
    figsize :
    sheet_name :
    ncol :
    fontsize :
    legend_fontsize : int
        legend fontsize, default 5
    legend_kws: dict
        kwargs passed to ax.legend
    legend_title_fontsize: int
        legend title fontsize, default 5
    marker_fontsize: int
        Marker fontsize, default 3
        if id_marker is True, and coding is True. legend marker will be a circle (or rectangle) with code
    linewidth : float
        Line width of the legend marker (circle or rectangle), default 0.5
    text_anon: bool
        Whether to add text annotation for each cluster, default False.
	ret : bool 
		Whether to return the plot and colors, default False.
    kwargs : dict
        set text_anno=None to plot clustering without text annotations,
        coding=True to plot clustering without code annotations,
        set show_legend=False to remove the legend

    Returns
    -------

    """
    if sheet_name is None: # for getting color scheme from excel file
        sheet_name=cluster_col
    if isinstance(adata,str): # getting adata
        adata=ad.read_h5ad(adata,backed='r')
    if not isinstance(adata.obs[cluster_col].dtype, pd.CategoricalDtype): # make sure cluster_col is categorical 
        adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')
    # get palette
    if palette_path is not None:
        if isinstance(palette_path,str):
            colors=pd.read_excel(os.path.expanduser(palette_path),sheet_name=sheet_name,index_col=0).Hex.to_dict()
            keys=list(colors.keys())
            existed_vals=adata.obs[cluster_col].unique().tolist()
            for k in existed_vals:
                if k not in keys:
                    colors[k]='gray'
            for k in keys:
                if k not in existed_vals:
                    del colors[k]
        else:
            colors=palette_path
        adata.uns[cluster_col + '_colors'] = [colors.get(k, 'grey') for k in adata.obs[cluster_col].cat.categories.tolist()]
    else:
        if f'{cluster_col}_colors' not in adata.uns:
            colors = add_color_scheme(adata, cluster_col, palette_key=f"{cluster_col}_colors")
        else:
            colors={cluster:color for cluster,color in zip(adata.obs[cluster_col].cat.categories.tolist(),adata.uns[f'{cluster_col}_colors'])}

    hue=cluster_col
    text_anno = cluster_col if text_anno else None
    text_kws = {} if text_kws is None else text_kws
    text_kws.setdefault("fontsize", fontsize)
    kwargs.setdefault("hue",hue)
    kwargs.setdefault("text_anno", text_anno)
    kwargs.setdefault("text_kws", text_kws)
    kwargs.setdefault("luminance", 0.65)
    kwargs.setdefault("dodge_text", False)
    kwargs.setdefault("axis_format", axis_format)
    kwargs.setdefault("show_legend", True)
    kwargs.setdefault("marker_fontsize", marker_fontsize)
    kwargs.setdefault("marker_pad", marker_pad)
    kwargs.setdefault("linewidth", linewidth)
    kwargs.setdefault("alpha", alpha)
    kwargs["coding"]=coding
    kwargs["id_marker"]=id_marker
    legend_kws={} if legend_kws is None else legend_kws
    default_lgd_kws=dict(
        fontsize=legend_fontsize,
        title=cluster_col,title_fontsize=legend_title_fontsize)
    if ncol is not None:
        default_lgd_kws['ncol']=ncol
    for k in default_lgd_kws:
        legend_kws.setdefault(k, default_lgd_kws[k])
    kwargs.setdefault("dodge_kws", {
            "arrowprops": {
                "arrowstyle": "->",
                "fc": 'grey',
                "ec": "none",
                "connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
            },
            'autoalign': 'xy'})
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
    p = categorical_scatter(
        data=adata[adata.obs[cluster_col].notna(),],
        ax=ax,
        coord_base=coord_base,
        palette=colors,
        legend_kws=legend_kws,
        **kwargs)

    if output is not None:
        plt.savefig(os.path.expanduser(output),bbox_inches='tight',dpi=300)
    if show:
        plt.show()
    
    if ret:
        return p, colors
    else: 
        return None, None