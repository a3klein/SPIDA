# Utility functions for plotting single cell data in python. 
# Functions come from ALLCools, and Pym3c by Hanqing Liu and Wubin Ding 
# Adapted by Amit Klein

import pandas as pd
import anndata as ad
import numpy as np
from matplotlib.colors import (
    Normalize, colorConverter, TwoSlopeNorm, Colormap
)
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text
import matplotlib.pyplot as plt

from .palettes import add_color_scheme

mm2inch = 1 / 25.4

# Custom handler for legend: circle text as marker + label
class TextWithCircleHandler(HandlerBase):
	def __init__(self, marker_text='', label_text='', 
			  text_kws={}, **kwargs):
		HandlerBase.__init__(self, **kwargs)
		self.marker_text = marker_text
		self.text_kws=text_kws

	def create_artists(self, legend, orig_handle,
						xdescent, ydescent, width, height, fontsize, trans):
		# Marker (number with circle)
		self.text_kws.setdefault("fontsize",fontsize)
		# print(self.text_kws)
		shift=2 * self.text_kws['fontsize'] * 0.65 / 72 / mm2inch
		circ_text = Text(
			xdescent + legend.borderaxespad + shift,  height / 2,
			self.marker_text, 
			**self.text_kws
		)
		return [circ_text]

def _extract_coords(
	data, 
	coord_base : str,
	x : str,
	y : str
):
	""" Given a data object (AnnData or xarray or Pandas Dataframe), extract the coordinates for plotting."""

	# import xarray as xr
	import pandas as pd
	if (x is None) or (y is None):
		x = f"{coord_base}_0"
		y = f"{coord_base}_1"

	if isinstance(data, ad.AnnData):
		adata = data
		if coord_base in adata.obsm: 
			_data = pd.DataFrame(
				{
                    "x": adata.obsm[f"{coord_base}"][:, 0],
                    "y": adata.obsm[f"{coord_base}"][:, 1],
                },
                index=adata.obs_names,
            )
		elif f"X_{coord_base}" in adata.obsm: 
			_data = pd.DataFrame(
                {
                    "x": adata.obsm[f"X_{coord_base}"][:, 0],
                    "y": adata.obsm[f"X_{coord_base}"][:, 1],
                },
                index=adata.obs_names,
            )
	# Don't need xarray stuff for spatialdata objects right now
	# elif isinstance(data, xr.Dataset):
	# 	ds = data
	# 	if coord_base not in ds.dims:
	# 		raise KeyError(f"xr.Dataset do not contain {coord_base} dim")
	# 	data_var = {i for i in ds.data_vars.keys() if i.startswith(coord_base)}.pop()
	# 	_data = pd.DataFrame(
	# 		{
	# 			"x": ds[data_var].sel({coord_base: f"{coord_base}_0"}).to_pandas(),
	# 			"y": ds[data_var].sel({coord_base: f"{coord_base}_1"}).to_pandas(),
	# 		}
	# 	)
	else:
		if (x not in data.columns) or (y not in data.columns):
			raise KeyError(f"{x} or {y} not found in columns.")
		_data = pd.DataFrame({"x": data[x], "y": data[y]})
	return _data, x, y

def _density_based_sample(
	data: pd.DataFrame,
	coords: list,
	portion=None,
	size=None,
	seed=None
):
	"""Down sample data based on density, to prevent overplot in dense region and decrease plotting time."""
	from sklearn.neighbors import LocalOutlierFactor
	clf = LocalOutlierFactor(
		n_neighbors=20,
		algorithm="auto",
		leaf_size=30,
		metric="minkowski",
		p=2,
		metric_params=None,
		contamination=0.1,
	)

	# coords should already exist in data, get them by column names list
	data_coords = data[coords]
	clf.fit(data_coords)
	# original score is negative, the larger the denser
	density_score = clf.negative_outlier_factor_
	delta = density_score.max() - density_score.min()
	# density score to probability: the denser the less probability to be picked up
	probability_score = 1 - (density_score - density_score.min()) / delta
	probability_score = np.sqrt(probability_score)
	probability_score = probability_score / probability_score.sum()

	if size is not None:
		pass
	elif portion is not None:
		size = int(data_coords.index.size * portion)
	else:
		raise ValueError("Either portion or size should be provided.")
	if seed is not None:
		np.random.seed(seed)
	selected_cell_index = np.random.choice(
		data_coords.index, size=size, replace=False, p=probability_score
	)  # choice data based on density weights

	# return the down sampled data
	return data.reindex(selected_cell_index)

def _auto_size(ax, n_dots):
    """Auto determine dot size based on ax size and n dots"""
    bbox = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted())
    scale = bbox.width * bbox.height / 14.6  # 14.6 is a 5*5 fig I used to estimate
    n = n_dots / scale  # larger figure means data look sparser
    if n < 500:
        s = 14 - n / 100
    elif n < 1500:
        s = 7
    elif n < 3000:
        s = 5
    elif n < 8000:
        s = 3
    elif n < 15000:
        s = 2
    elif n < 30000:
        s = 1.5
    elif n < 50000:
        s = 1
    elif n < 80000:
        s = 0.8
    elif n < 150000:
        s = 0.6
    elif n < 300000:
        s = 0.5
    elif n < 500000:
        s = 0.4
    elif n < 800000:
        s = 0.3
    elif n < 1000000:
        s = 0.2
    elif n < 2000000:
        s = 0.1
    elif n < 3000000:
        s = 0.07
    elif n < 4000000:
        s = 0.05
    elif n < 5000000:
        s = 0.03
    else:
        s = 0.02
    return s

def zoom_min_max(vmin, vmax, scale):
    """Zoom min and max value."""
    width = vmax - vmin
    width_zoomed = width * scale
    delta_value = (width_zoomed - width) / 2
    return vmin - delta_value, vmax + delta_value

def zoom_ax(ax, zoom_scale, on="both"):
    """Zoom ax on both x and y-axis."""
    on = on.lower()
    xlim = ax.get_xlim()
    xlim_zoomed = zoom_min_max(vmin=xlim[0], vmax=xlim[1],scale=zoom_scale)

    ylim = ax.get_ylim()
    ylim_zoomed = zoom_min_max(vmin=ylim[0], vmax=ylim[1],scale=zoom_scale)

    if (on == "both") or ("x" in on):
        ax.set_xlim(xlim_zoomed)
    if (on == "both") or ("y" in on):
        ax.set_ylim(ylim_zoomed)
		
def _take_data_series(data, k):
	"""Gets the column from the data object (for xarray, AnnData or Pandas)"""
	# Don't need xarray stuff for spatialdata objects right now
	# import xarray as xr
	# if isinstance(data, (xr.Dataset, xr.DataArray)):
	# 	_value = data[k].to_pandas()
	if isinstance(data, ad.AnnData):
		_value = data.obs[k].copy()
	else:
		_value = data[k].copy()
	return _value

# This is a good function
def level_one_palette(
	name_list,
	order=None,
	palette="auto"
):
	"""Generates a color palette for a list of names, with an option to specify the order."""
	import seaborn as sns

	name_set = set(name_list.dropna())
	# TODO: Custom palettes based on Scanpy default palettes? 
	if palette == "auto":
		if len(name_set) < 10:
			palette = "tab10"
		elif len(name_set) < 20:
			palette = "tab20"
		else:
			palette = "rainbow"

	if order is None:
		try:
			order = sorted(name_set)
		except TypeError:
			# name set contains multiple dtype (e.g., str and np.NaN)
			order = list(name_set)
	else:
		if (set(order) != name_set) or (len(order) != len(name_set)):
			raise ValueError("Order is not equal to set(name_list).")

	n = len(order)
	colors = sns.color_palette(palette, n)
	color_palette = {}
	for name, color in zip(order, colors):
		color_palette[name] = color
	return color_palette

def _calculate_luminance(color):
	"""
	Calculate the relative luminance of a color according to W3C standards

	Parameters
	----------
	color : matplotlib color or sequence of matplotlib colors
		Hex code, rgb-tuple, or html color name.
	Returns
	-------
	luminance : float(s) between 0 and 1

	"""
	rgb = colorConverter.to_rgba_array(color)[:, :3]
	rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
	lum = rgb.dot([0.2126, 0.7152, 0.0722])
	try:
		return lum.item()
	except ValueError:
		return lum
	
def tight_hue_range(hue_data, portion):
	"""Automatic select a SMALLEST data range that covers [portion] of the data."""
	hue_data = hue_data[np.isfinite(hue_data)]
	hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
	min_window_right = (
		hue_quantiles.rolling(window=int(portion * 100)).apply(lambda i: i.max() - i.min(), raw=True).idxmin()
	)
	min_window_left = max(0, min_window_right - portion)
	vmin, vmax = tuple(hue_data.quantile(q=[min_window_left, min_window_right]))
	if np.isfinite(vmin):
		vmin = max(hue_data.min(), vmin)
	else:
		vmin = hue_data.min()
	if np.isfinite(vmax):
		vmax = min(hue_data.max(), vmax)
	else:
		vmax = hue_data.max()
	return vmin, vmax

# From ALLCools - Need to install adjust_text to use the dodge_text feature
def _text_anno_scatter(
    data: pd.DataFrame,
    ax,
    x: str,
    y: str,
    dodge_text=False,
    anno_col="text_anno",
    text_kws=None,
    text_transform=None,
    dodge_kws=None,
    luminance=0.48
):
	"""Add text annotation to a scatter plot."""
	import copy
	# prepare kws
	text_kws={} if text_kws is None else text_kws
	text_kws.setdefault("fontsize",5)
	text_kws.setdefault("fontweight","black")
	text_kws.setdefault("ha","center") #horizontalalignment
	text_kws.setdefault("va","center") #verticalalignment
	text_kws.setdefault("color","black") #c
	bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
								facecolor=(0.8, 0.8, 0.8, 0.2),alpha=1,linewidth=0.5)
	text_kws.setdefault("bbox",bbox)
	for key in bbox:
		if key not in text_kws['bbox']:
			text_kws['bbox'][key]=bbox[key]
	# plot each text
	text_list = []
	for text, sub_df in data.groupby(anno_col):
		if text_transform is None:
			text = str(text)
		else:
			text = text_transform(text)
		if text.lower() in ["", "nan"]:
			continue
		_x, _y = sub_df[[x, y]].median()
		
		use_text_kws=copy.deepcopy(text_kws) #text_kws.copy()
		if isinstance(text_kws['bbox']['facecolor'],dict):
			use_text_kws['bbox']['facecolor']=text_kws['bbox']['facecolor'].get(text,'gray')
		if isinstance(text_kws['color'],dict):
			use_color=text_kws['color'].get(text,'black')
			use_text_kws['color']=use_color
		if not luminance is None and not use_text_kws['bbox']['facecolor'] is None:
			lum = _calculate_luminance(use_text_kws['bbox']['facecolor'])
			if lum > luminance:
				use_text_kws['color']='black'
				use_text_kws['bbox']['edgecolor']='black'
		
		text = ax.text(
			_x,
			_y,
			text,
			**use_text_kws
		)
		text_list.append(text)

	if dodge_text:
		try:
			from adjustText import adjust_text

			_dodge_parms = {
				"force_points": (0.02, 0.05),
				"arrowprops": {
					"arrowstyle": "->",
					"fc": "black",
					"ec": "none",
					"connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
				},
				"autoalign": "xy",
			}
			if dodge_kws is not None:
				_dodge_parms.update(dodge_kws)
			adjust_text(text_list, x=data["x"], y=data["y"], **_dodge_parms)
		except ModuleNotFoundError:
			print("Install adjustText package to dodge text, see its github page for help")
	return

def density_contour(
    ax,
    data,
    x,
    y,
    groupby=None,
    c="lightgray",
    single_contour_pad=1,
    linewidth=1,
    palette=None,
):
	"""Add Descriptor"""
	from sklearn.neighbors import LocalOutlierFactor
	_data = data.copy()

	if groupby is not None:
		if isinstance(groupby, str):
			_data["groupby"] = data[groupby]
		else:
			_data["groupby"] = groupby
	else:
		_data["groupby"] = "one group"

	_contour_kws = {"linewidths": linewidth, "levels": (-single_contour_pad,), "linestyles": "dashed"}
	_lof_kws = {"n_neighbors": 25, "novelty": True, "contamination": "auto"}

	xmin, ymin = _data[[x, y]].min()
	xmax, ymax = _data[[x, y]].max()
	xmin, xmax = zoom_min_max(xmin, xmax, 1.2)
	ymin, ymax = zoom_min_max(ymin, ymax, 1.2)

	for group, sub_data in _data[[x, y, "groupby"]].groupby("groupby"):
		xx, yy = np.meshgrid(np.linspace(xmin, xmax, 500), np.linspace(ymin, ymax, 500))
		clf = LocalOutlierFactor(**_lof_kws)
		clf.fit(sub_data.iloc[:, :2].values)
		z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
		z = z.reshape(xx.shape)
		if palette is None:
			_color = c
		else:
			_color = palette[group] if group in palette else c
		# plot contour line(s)
		ax.contour(xx, yy, z, colors=_color, **_contour_kws)
	return

def despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False):
	"""
	Remove the top and right spines from plot(s).

	Parameters
	----------
	fig : matplotlib figure, optional
		Figure to despine all axes of, defaults to the current figure.
	ax : matplotlib axes, optional
		Specific axes object to despine. Ignored if fig is provided.
	top, right, left, bottom : boolean, optional
		If True, remove that spine.

	Returns
	-------
	None

	"""
	if fig is None and ax is None:
		axes = plt.gcf().axes
	elif fig is not None:
		axes = fig.axes
	elif ax is not None:
		axes = [ax]

	for ax_i in axes:
		for side in ["top", "right", "left", "bottom"]:
			is_visible = not locals()[side]
			ax_i.spines[side].set_visible(is_visible)
		if left and not right:  # remove left, keep right
			maj_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.majorTicks)
			min_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.minorTicks)
			ax_i.yaxis.set_ticks_position("right")
			for t in ax_i.yaxis.majorTicks:
				t.tick2line.set_visible(maj_on)
			for t in ax_i.yaxis.minorTicks:
				t.tick2line.set_visible(min_on)

		if bottom and not top:
			maj_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.majorTicks)
			min_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.minorTicks)
			ax_i.xaxis.set_ticks_position("top")
			for t in ax_i.xaxis.majorTicks:
				t.tick2line.set_visible(maj_on)
			for t in ax_i.xaxis.minorTicks:
				t.tick2line.set_visible(min_on)


def _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5):
	"""Make a tiny axis label in the bottom right corner of the axis"""
	# This function assume coord is [0, 1].
	ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
	despine(ax=ax, left=True, bottom=True)

	_arrow_kws = {"width": 0.003, "linewidth": 0, "color": "black"}
	if arrow_kws is not None:
		_arrow_kws.update(arrow_kws)

	ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws, transform=ax.transAxes)
	ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws, transform=ax.transAxes)
	ax.text(
		0.06,
		0.03,
		x.upper().replace("_", " "),
		fontdict={"fontsize": fontsize, "horizontalalignment": "left", "verticalalignment": "center"},
		transform=ax.transAxes,
	)
	ax.text(
		0.03,
		0.06,
		y.upper().replace("_", " "),
		fontdict={
			"fontsize": fontsize,
			"rotation": 90,
			"rotation_mode": "anchor",
			"horizontalalignment": "left",
			"verticalalignment": "center",
		},
		transform=ax.transAxes,
	)
	return


def plot_text_legend(color_dict, code2label, ax=None, title=None, color_text=True, 
					 boxstyle='Circle',marker_pad=0.1,legend_kws=None,marker_fontsize=4,
					 text_kws=None,alpha=0.7,luminance=0.5):
	import copy
	# print(color_dict)
	lgd_kws = legend_kws.copy() if not legend_kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
	# lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.5
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.2
	)  # gap height between two row of legend,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)
	ms = lgd_kws.pop("markersize", 10)

	# text_kws
	if text_kws is None:
		text_kws={}
	default_marker_text_kws=dict(
		bbox=dict(boxstyle=f"{boxstyle},pad={marker_pad}", #Square, Circle, Round
			edgecolor='black',linewidth=0.4,
			fill=True,facecolor='white',alpha=alpha),
		horizontalalignment='center', verticalalignment='center',
		fontsize=marker_fontsize,color='black')
	for k in default_marker_text_kws:
		if k == 'bbox':
			if k not in text_kws:
				text_kws['bbox']=default_marker_text_kws[k]
			else:
				for k1 in default_marker_text_kws['bbox'].keys():
					if k1 not in text_kws['bbox']:
						text_kws['bbox'].setdefault(k1,default_marker_text_kws['bbox'][k1])
		if k not in text_kws:
			text_kws.setdefault(k,default_marker_text_kws[k])
	
	# Create handles and handlers
	handles = []
	handler_map = {}
	for code in sorted([int(i) for i in code2label.keys()]):
		label=code2label[code]
		code_text=str(code)
		handle = Line2D([], [], linestyle=None,
				  label=label)
		handles.append(handle)
		color=color_dict.get(label, 'black')
		text_kws1= copy.deepcopy(text_kws)
		text_kws1['bbox']['facecolor']=color
		lum = _calculate_luminance(color)
		if lum <= 0.1: # for black-like color, use white marker text
			text_kws1['color']='white'
		# print(bbox)
		handler_map[handle] = TextWithCircleHandler(
			marker_text=code_text, label_text=label,
			text_kws=text_kws1,
			)

	# Draw custom legend
	L=ax.legend(handles=handles, handler_map=handler_map, 
			 **lgd_kws)
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	ax.figure.canvas.draw()
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(color_dict[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else color_dict[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(lgd)
	# ax.grid(False)

def plot_color_dict_legend(
	D, ax=None, title=None, color_text=True, 
	kws=None,luminance=0.5
):
	"""
	plot legned for color dict

	Parameters
	----------
	D: a dict, key is categorical variable, values are colors.
	ax: axes to plot the legend.
	title: title of legend.
	color_text: whether to change the color of text based on the color in D.
	label_side: right of left.
	kws: kws passed to plt.legend.

	Returns
	-------
	ax.legend

	"""
	import matplotlib.patches as mpatches
	if ax is None:
		ax = plt.gca()
	lgd_kws = kws.copy() if not kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
	lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.1
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.15
	)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)
	l = [
		mpatches.Patch(color=c, label=l) for l, c in D.items()
	]  # kws:?mpatches.Patch; rasterized=True
	ms = lgd_kws.pop("markersize", 10)
	L = ax.legend(handles=l, **lgd_kws)
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(D[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else D[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(L)
	ax.grid(False)
	return L

def plot_marker_legend(
	color_dict=None, ax=None, title=None, color_text=True, 
	marker='o',kws=None,luminance=0.5
):
	"""
	plot legned for different marker

	Parameters
	----------
	D: a dict, key is categorical variable, values are marker.
	ax: axes to plot the legend.
	title: title of legend.
	color_text: whether to change the color of text based on the color in D.
	label_side: right of left.
	kws: kws passed to plt.legend.

	Returns
	-------
	ax.legend

	"""
	if ax is None:
		ax = plt.gca()

	lgd_kws = kws.copy() if not kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("borderpad", 0.2 * mm2inch * 72)  # 0.1mm
	# lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.1
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4 #0.2 * mm2inch * 72
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.15
	)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)

	ms = lgd_kws.pop("markersize", 10)
	L = [
		Line2D(
			[0],
			[0],
			color=color,
			marker=marker,
			linestyle="None",
			markersize=ms,
			label=l,
		)
		for l,color in color_dict.items()
	]
	L = ax.legend(handles=L, **lgd_kws)
	ax.figure.canvas.draw()
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(color_dict[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else color_dict[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(lgd)
	ax.grid(False)
	return L

def plot_cmap_legend(
	cax=None, ax=None, cmap="turbo", label=None, kws=None,
	labelsize=6, linewidth=0.5,ticklabel_size=4,
):
	"""
	Plot legend for cmap.

	Parameters
	----------
	cax : Axes into which the colorbar will be drawn.
	ax :  axes to anchor.
	cmap : turbo, hsv, Set1, Dark2, Paired, Accent,tab20,exp1,exp2,meth1,meth2
	label : title for legend.
	kws : dict
		kws passed to plt.colorbar (matplotlib.figure.Figure.colorbar).

	Returns
	-------
	cbar: axes of legend

	"""
	label = "" if label is None else label
	cbar_kws = {} if kws is None else kws.copy()
	cbar_kws.setdefault("label", label)
	# cbar_kws.setdefault("aspect",3)
	cbar_kws.setdefault("orientation", "vertical")
	# cbar_kws.setdefault("use_gridspec", True)
	# cbar_kws.setdefault("location", "bottom")
	cbar_kws.setdefault("fraction", 1)
	cbar_kws.setdefault("shrink", 1)
	cbar_kws.setdefault("pad", 0)
	cbar_kws.setdefault("extend", 'both')
	cbar_kws.setdefault("extendfrac", 0.1)
	# print(cbar_kws,kws)
	# print(type(cax))
	vmax = cbar_kws.pop("vmax", 1)
	vmin = cbar_kws.pop("vmin", 0)
	# print(vmin,vmax,'vmax,vmin')
	cax.set_ylim([vmin, vmax])
	# print(cax.get_ylim())
	vcenter= (vmax + vmin) / 2
	center=cbar_kws.pop("center",None)
	if center is None:
		center=vcenter
		m = plt.cm.ScalarMappable(
			norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap
		)
	else:
		m = plt.cm.ScalarMappable(
			norm=TwoSlopeNorm(center,vmin=vmin, vmax=vmax), cmap=cmap
		)
	cbar_kws.setdefault("ticks", [vmin, center, vmax])
	cax.yaxis.set_label_position('right')
	cax.yaxis.set_ticks_position('right')
	cbar = ax.figure.colorbar(m, cax=cax, **cbar_kws)  # use_gridspec=True
	cbar.ax.tick_params(labelsize=ticklabel_size, size=ticklabel_size,width=linewidth) # size is for ticks, labelsize is for the number on ticks (ticklabels)
	cbar.ax.yaxis.label.set_fontsize(labelsize) # colorbar title fontsize
	cbar.ax.grid(False)
	return cbar

def get_cmap(cmap):
	try:
		ret = plt.colormaps.get(cmap)  # matplotlib >= 3.5.1?
	except:
		ret = plt.get_cmap(cmap)  # matplotlib <=3.4.3?
	if ret is None: 
		raise ValueError(f"cmap {cmap} not found, please check the name.")
	return ret

### Added by Amit: 
def fig_plot_marker_legend(
	color_dict=None, ax=None, title=None, color_text=True, 
	marker='o',kws=None,luminance=0.5
):
	"""
	plot legend for different marker

	Parameters
	----------
	D: a dict, key is categorical variable, values are marker.
	ax: axes to plot the legend.
	title: title of legend.
	color_text: whether to change the color of text based on the color in D.
	label_side: right of left.
	kws: kws passed to plt.legend.

	Returns
	-------
	ax.legend

	"""
	if ax is None:
		ax = plt.gca()

	lgd_kws = kws.copy() if not kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "center left"
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("borderpad", 0.2 * mm2inch * 72)  # 0.1mm
	# lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.1
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4 #0.2 * mm2inch * 72
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.5
	)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	# lgd_kws.setdefault("bbox_to_anchor", (0, 1))
	lgd_kws["bbox_to_anchor"] = (0.5, 0.5)
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)

	ms = lgd_kws.pop("markersize", 10)
	L = [
		Line2D(
			[0],
			[0],
			color=color,
			marker=marker,
			linestyle="None",
			markersize=ms,
			label=l,
		)
		for l,color in color_dict.items()
	]
	L = ax.legend(handles=L, **lgd_kws)
	ax.figure.canvas.draw()
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(color_dict[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else color_dict[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	ax.grid(False)
	ax.axis("off")
	return L

def add_legend_axis(
    fig,
    legend_width=3, 
    legend_height=0.6,
    legend_ystart=0.2,
    space=0
): 
    """ 
	Add a legend axis to the right of the current axes in a figure.

	Parameters
	----------
	fig : matplotlib figure
		The figure to which the legend axis will be added.
	legend_width : float, optional
		The width of the legend axis in mm, by default 3 mm.
	legend_height : float, optional		
		The height of the legend axis in fraction of the figure height, by default 0.6.
	legend_ystart : float, optional
		The starting y position of the legend axis in fraction of the figure height, by default 0.2.
	space : float, optional
		The space in mm between the rightmost axis and the legend axis, by default 0.0 mm.
	"""

    ax = fig.get_axes()[0]
    mm2inch = 1 / 25.4
    
    legend_width = (legend_width * mm2inch * fig.dpi / fig.get_window_extent().width)  # mm to px to fraction
    pad = (space + ax.yaxis.labelpad * 1.2 * fig.dpi / 72) / fig.get_window_extent().width # labelpad unit is points
    # Iterate through all axes to find the rightmost position
    rightmost_x = 0
    for _ax in fig.get_axes():
        bbox = _ax.get_position()
        rightmost_x = max(rightmost_x, bbox.x1)
    left = rightmost_x + pad
    # Creating the new legend
    ax_legend = fig.add_axes(
        [left, legend_ystart, legend_width, legend_height]
    )  # left, bottom, width, height
    return ax_legend

def fig_make_legend(
    fig,
    _data, 
    hue,
    palette_dict,
    alpha=0.7,
    code2label=None, 
    id_marker=True, 
    labelsize=4,
    legend_color_text=True,
    legend_kws=None,
    legend_height=0.6,
    legend_width=3, 
    legend_ystart=0.2,
    luminance=None,	
    marker_pad=0.1,
    marker_fontsize=4,
    rectangle_marker=False,     
    space=0
):
    """
    Adding legend to the axis.

    Parameters
    ----------
    fig
       figure to add legend.
    _data
        data to plot legend.
    hue
        hue name for legend.
    palette_dict
        palette dictionary for hue.
    alpha
        alpha value for legend (default 0.7).
    code2label
        dictionary mapping code to label for legend text (default None).
    id_marker
        whether to use id marker for legend (default True).
    labelsize
        label size for legend text (default 4).
    legend_color_text
        whether to use color text for legend (default True).
    legend_kws
		kwargs for legend (default None).
	legend_height
		height of legend in fraction of figure height (default 0.6).
	legend_width
		width of legend in mm (default 3).
	legend_ystart
		starting y position of legend in fraction of figure height (default 0.2).
    luminance
        luminance value for legend (default None).
    marker_pad
        padding between markers in legend (default 0.1).
    marker_fontsize
        font size for markers in legend (default 4).
    rectangle_marker
        whether to use rectangle marker for legend (default False).
	space
		space in mm between the rightmost axis and the legend axis (default 0).
    """
    
    
    n_hue = len(palette_dict)
    ncol=1 if n_hue <= 40 else 2 if n_hue <= 100 else 3
    if legend_kws is None:
        legend_kws = {}
    default_lgd_kws = dict(
        ncol=ncol,fontsize=labelsize,
        bbox_to_anchor=(1, 1),loc="upper left",
        markersize=labelsize #legend_kws["fontsize"],
    )
    for k in default_lgd_kws:
        legend_kws.setdefault(k, default_lgd_kws[k])

    exist_hues = _data["hue"].unique()
    color_dict={hue_name: color for hue_name, color in palette_dict.items() if hue_name in exist_hues}
    
    ax_legend = add_legend_axis(
		fig=fig, legend_width=legend_width, legend_height=legend_height, legend_ystart=legend_ystart, space=space
	)
    if code2label is not None and id_marker:
        boxstyle='Circle' if not rectangle_marker else 'Round'
        plot_text_legend(color_dict, code2label, ax_legend, title=hue, 
                color_text=legend_color_text, boxstyle=boxstyle,marker_pad=marker_pad,
                legend_kws=legend_kws,marker_fontsize=marker_fontsize,
                alpha=alpha,luminance=luminance)
    else:
        if rectangle_marker: 
			## plot Patch legend (rectangle marker)
            plot_color_dict_legend(
                D=color_dict, ax=ax_legend, title=hue, color_text=legend_color_text, 
                kws=legend_kws,luminance=luminance
            )
        else:
            # plot marker legend (for example, circle marker)
            fig_plot_marker_legend(
                color_dict=color_dict, ax=ax_legend, title=hue, color_text=legend_color_text, 
                marker='o',kws=legend_kws,luminance=luminance
            )


# Colorbar for continuous scatter plot
def fig_make_colorbar(
    fig,
    vmin : float,
    vmax : float,
    label : str,
    cmap : str | Colormap = 'Wistia',
    cbar_width=3,
    legend_ystart=0.2, 
    legend_height=0.6,
    labelsize=6,
    linewidth=0.5,
    ticklabel_size=4,
    space=0,
    **cbar_kws
): 
    """ 
	Plot the colorbar for a given axis / figure 
	
	Parameters
    ----------
    fig
        figure to add legend.
    vmin
        minimum value for colorbar.
    vmax
        maximum value for colorbar.
	label
		label for colorbar.
    cmap
        colormap to use for colorbar (default 'Wistia').
	cbar_width
		width of colorbar in mm (default 3).    
	legend_height
		height of legend in fraction of figure height (default 0.6).
	legend_width
		width of legend in mm (default 3).
	legend_ystart
		starting y position of legend in fraction of figure height (default 0.2).
	linewidth
		linewidth for colorbar ticks (default 0.5).
	ticklabel_size
		size of tick labels for colorbar (default 4).
	space
		space in mm between the rightmost axis and the legend axis (default 0).
	cbar_kws
		kwargs for colorbar (default None).
    """
    ax = fig.get_axes()[0]
    
    # make color bar
    # default_cbar_kws=dict(loc="upper left", borderpad=0,width="3%", height="20%") #bbox_to_anchor=(1,1)
    if cbar_kws is None:
        cbar_kws={}
    
    ax_legend = add_legend_axis(
		fig=fig, legend_width=cbar_width, legend_height=legend_height, legend_ystart=legend_ystart, space=space
	)
    cbar_kws['vmin']=vmin
    cbar_kws['vmax']=vmax
    plot_cmap_legend(
        ax=ax,
        cax=ax_legend,
        cmap=cmap,
        label=label,
        kws=cbar_kws.copy(),labelsize=labelsize, 
        linewidth=linewidth,ticklabel_size=ticklabel_size,
        )