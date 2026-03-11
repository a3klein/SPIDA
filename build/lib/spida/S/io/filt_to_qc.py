from __future__ import annotations

import logging
from collections.abc import Iterable

import anndata as ad
import geopandas as gpd
import pandas as pd


logger = logging.getLogger(__package__)


# For inferring points cell-id column if not explicitly provided, in order of priority.
# Enables cross-compatibility across segmentation methods
_POINTS_CELL_ID_CANDIDATES = (
	"cell_id",
	"assignment",
	"cell",
	"EntityID",
	"CELL_ID",
)


# From the transcripts_qc, so geometries can be stored in anndata objects through zarr
def _obs_to_geodf(obs_df: pd.DataFrame) -> gpd.GeoDataFrame:
	"""Build a GeoDataFrame from AnnData `.obs` that stores geometry as WKT or shapely."""
	if "geometry_wkt" in obs_df.columns:
		geom = gpd.GeoSeries.from_wkt(obs_df["geometry_wkt"])
		base = obs_df.drop(columns=["geometry_wkt"]).copy()
		return gpd.GeoDataFrame(base, geometry=geom)
	if "geometry" in obs_df.columns:
		return gpd.GeoDataFrame(obs_df.copy(), geometry="geometry")
	raise ValueError(
		"Could not find geometry in QC table. Expected `geometry_wkt` or `geometry` in `.obs`."
	)


# function to get the qc regions as a geodataframe. 
def _resolve_qc_regions(
	sdata,
	qc_regions: gpd.GeoDataFrame | pd.DataFrame | ad.AnnData | None = None,
	qc_table_key: str | None = None,
	qc_shapes_key: str | None = None,
	qc_filter_col: str = "filtered",
	qc_pass_value: bool | int | str = False,
) -> gpd.GeoDataFrame:
	"""
	Resolve QC pass regions from a provided object or from elements already in `sdata`.

	Priority:
	1) `qc_regions` - an actual object
	2) `qc_shapes_key` from `sdata` - a shapes object in `sdata.shapes` 
	3) `qc_table_key` from `sdata.tables` - a table object in `sdata.tables` with geometry info in `.obs`
	"""
	source = qc_regions
	if source is None and qc_shapes_key is not None:
		source = sdata.shapes[qc_shapes_key]
	if source is None and qc_table_key is not None:
		source = sdata.tables[qc_table_key]

	if source is None:
		raise ValueError(
			"QC regions were not provided. Pass `qc_regions` or provide `qc_shapes_key`/`qc_table_key`."
		)
	# Based on what type of source it is, need to parse the geometry and filter by QC pass column if it exists
	if isinstance(source, ad.AnnData):
		obs = source.obs.copy()
		if qc_filter_col in obs.columns:
			obs = obs.loc[obs[qc_filter_col] == qc_pass_value]
		qc_gdf = _obs_to_geodf(obs)
	elif isinstance(source, gpd.GeoDataFrame):
		qc_gdf = source.copy()
		if qc_filter_col in qc_gdf.columns:
			qc_gdf = qc_gdf.loc[qc_gdf[qc_filter_col] == qc_pass_value]
	elif isinstance(source, pd.DataFrame):
		tmp = source.copy()
		if qc_filter_col in tmp.columns:
			tmp = tmp.loc[tmp[qc_filter_col] == qc_pass_value]
		qc_gdf = _obs_to_geodf(tmp)
	else:
		raise TypeError(
			"`qc_regions` must be a GeoDataFrame, DataFrame, AnnData, or a key into `sdata`."
		)
	# what's returned is just the geometries which pass the qc. So any intersect will work with these geoms. 
	if qc_gdf.empty:
		print("QC regions resolved to an empty GeoDataFrame.")
	return qc_gdf


# self explanatory, based on intersections. 
def _filter_shapes_by_qc(
	shapes: gpd.GeoDataFrame,
	qc_regions: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, set[str]]:
	"""Keep only shapes intersecting the QC pass regions and return kept cell ids."""
	if qc_regions.empty:
		empty = shapes.iloc[0:0].copy()
		return empty, set()

	qc_union = (
		qc_regions.geometry.union_all()
		if hasattr(qc_regions.geometry, "union_all")
		else qc_regions.geometry.unary_union
	)
	keep_mask = shapes.geometry.intersects(qc_union)
	shapes_filt = shapes.loc[keep_mask].copy()
	keep_ids = set(shapes_filt.index.astype(str))
	return shapes_filt, keep_ids



def _resolve_table_cell_series(table: ad.AnnData, table_cell_id_col: str | None = None) -> pd.Series:
	"""Resolve per-row cell ids from a table using explicit column, instance key, or index."""
	if table_cell_id_col is not None:
		if table_cell_id_col not in table.obs.columns:
			raise KeyError(
				f"`table_cell_id_col='{table_cell_id_col}'` not found in table.obs columns."
			)
		return table.obs[table_cell_id_col].astype(str)

	attrs = table.uns.get("spatialdata_attrs", {})
	instance_key = attrs.get("instance_key", None)
	if instance_key is not None and instance_key in table.obs.columns:
		return table.obs[instance_key].astype(str)

	return pd.Series(table.obs_names.astype(str), index=table.obs_names)


def _filter_table_by_cells(
	table: ad.AnnData,
	keep_cell_ids: Iterable[str],
	table_cell_id_col: str | None = None,
) -> ad.AnnData:
	"""Subset table rows to cells that remain after QC filtering."""
	keep_cell_ids = set(str(x) for x in keep_cell_ids) # unique cell ids (bc repeat in shapes (7 z's))
	# Gets the matching cell ids from the anndata object.
	cell_series = _resolve_table_cell_series(table, table_cell_id_col=table_cell_id_col)
	keep_mask = cell_series.isin(keep_cell_ids).to_numpy()
	return table[keep_mask].copy()

# This is hacky, but kudos to codex for making this not a problem with different segmentation methods.
def _resolve_points_cell_id_col(points, points_cell_id_col: str | None = None) -> str:
	"""Infer points cell-id column if not explicitly provided."""
	cols = list(points.columns)
	if points_cell_id_col is not None:
		if points_cell_id_col not in cols:
			raise KeyError(
				f"`points_cell_id_col='{points_cell_id_col}'` not found in points columns: {cols}"
			)
		return points_cell_id_col

	for col in _POINTS_CELL_ID_CANDIDATES:
		if col in cols:
			return col

	raise KeyError(
		"Could not infer points cell-id column. "
		f"Tried {_POINTS_CELL_ID_CANDIDATES} and found columns {cols}."
	)

def _set_unassigned_cells(points, points_cell_id_col: str, keep_cell_ids: Iterable[str]):
    """
    Set points assigned to filtered-out cells as unassigned (`-1`).

    Supports pandas and dask dataframes.
    """
    keep_ids = set(str(x) for x in keep_cell_ids)

    def _update_partition(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        s = out[points_cell_id_col]
        is_numeric = pd.api.types.is_numeric_dtype(s)

        if is_numeric:
            s_num = pd.to_numeric(s, errors="coerce")
            keep_mask = s_num.astype("Int64").astype(str).isin(keep_ids)
            out[points_cell_id_col] = s_num.where(keep_mask, other=-1)
            return out

        s_str = s.astype("string")
        keep_mask = s_str.isin(keep_ids)
        out[points_cell_id_col] = s_str.where(keep_mask, other="-1")
        return out

    # In case of a dask dataframe, doing this in an efficient manner
    if hasattr(points, "map_partitions"):
        meta = points._meta.copy()
        return points.map_partitions(_update_partition, meta=meta)
    return _update_partition(points)

# The runner method
def apply_qc_filter_to_segmentation(
	sdata,
	*,
	table_key: str,
	shapes_key: str,
	points_key: str,
	qc_regions: gpd.GeoDataFrame | pd.DataFrame | ad.AnnData | None = None,
	qc_table_key: str | None = None,
	qc_shapes_key: str | None = None,
	qc_filter_col: str = "filtered",
	qc_pass_value: bool | int | str = False,
	table_cell_id_col: str | None = None,
	points_cell_id_col: str | None = None,
):
    """
    Filter segmentation elements by transcript-density QC regions.

    Effects:
    - Removes shapes outside QC pass regions.
    - Removes corresponding cells from `table`.
    - Sets points assigned to removed cells as unassigned (`-1`) in points cell-id column.
    """
	# Getting only the regions (hexes) that passed qc
    qc_gdf = _resolve_qc_regions(
		sdata,
		qc_regions=qc_regions,
		qc_table_key=qc_table_key,
		qc_shapes_key=qc_shapes_key,
		qc_filter_col=qc_filter_col,
		qc_pass_value=qc_pass_value,
	)

    # Extracting features
    shapes = sdata[shapes_key].copy()
    table = sdata[table_key].copy()
    points = sdata[points_key]

    n_shapes_before = shapes.shape[0]
    n_table_before = table.n_obs

    # Filtering steps (from above functions)
    shapes_filt, keep_cell_ids = _filter_shapes_by_qc(shapes, qc_gdf)
    table_filt = _filter_table_by_cells(
        table,
        keep_cell_ids=keep_cell_ids,
        table_cell_id_col=table_cell_id_col,
    )
    points_cell_id_col = _resolve_points_cell_id_col(points, points_cell_id_col)
    points_filt = _set_unassigned_cells(
        points,
        points_cell_id_col=points_cell_id_col,
        keep_cell_ids=keep_cell_ids,
    )

    # Loading back into the spatialdata object
    sdata[shapes_key] = shapes_filt
    sdata[table_key] = table_filt
    sdata[points_key] = points_filt

    # report
    logger.info(
        "Applied QC segmentation filter (%s): shapes %s -> %s, table cells %s -> %s, kept_ids=%s",
        shapes_key,
        n_shapes_before,
        shapes_filt.shape[0],
        n_table_before,
        table_filt.n_obs,
        len(keep_cell_ids),
    )
	
    return sdata