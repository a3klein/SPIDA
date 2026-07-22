"""Spec-driven readers for loading a segmenter's post-processed output into a
SpatialData object. Collapses the old ``read_vpt.py`` + ``read_proseg.py`` into
one set of helpers parametrized by the cell-id column.

Post-Stage-2 the process step writes segmentation schema files
(``boundaries_micron.parquet``, ``cell_by_gene.csv``, ``detected_transcripts.csv``,
``cell_metadata.csv``) all keyed on ``EntityID`` (proseg ``EntityID == cell``), so
these readers key everything on ``EntityID`` uniformly. proseg additionally keeps
its native counts/transcripts in ``proseg_outputs.zarr`` (v3) or CSVs (v2), read
via the assignment-native path.
"""

from __future__ import annotations

import gzip
from pathlib import Path
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from spatialdata.models import PointsModel, ShapesModel, TableModel
    from spatialdata.transformations import BaseTransformation

import pandas as pd
import anndata as ad
import geopandas as gpd
from shapely.geometry import MultiPolygon

from spida._constants import DEFAULT_PRESET, PROSEG_PRESET, REGION_KEY, GEOMETRY_COL

# cell-id key across all segmentation schema files (proseg cell id == EntityID)
_CELL_KEY = DEFAULT_PRESET["METADATA_CELL_KEY"]  # "EntityID"


def _cast_to_mp(x):
    try:
        return MultiPolygon(x.geoms)
    except Exception:
        return MultiPolygon([x.convex_hull])


def _read_boundaries(boundaries_path: str | Path) -> gpd.GeoDataFrame:
    """Read a boundaries file (parquet or [gzipped] geojson) into a GeoDataFrame."""
    path = str(boundaries_path)
    if path.endswith(".parquet"):
        gdf = gpd.read_parquet(path)
    elif path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            gdf = gpd.read_file(f)
    else:
        gdf = gpd.read_file(path)
    return gdf


def get_polygons(
    boundaries_path: str | Path,
    transformations: dict[str, BaseTransformation],
    *,
    cell_id_col: str = _CELL_KEY,
) -> gpd.GeoDataFrame:
    """Parse a boundaries file into a ShapesModel, indexed by ``cell_id_col`` (str).

    Multi-z segmentation schema boundaries are collapsed downstream by
    ``_cast_multipolygons_to_polygons`` (dedup-to-largest per cell id).
    """
    gdf = _read_boundaries(boundaries_path)
    if gdf.geometry.name != GEOMETRY_COL:
        gdf = gdf.rename_geometry(GEOMETRY_COL)
    gdf.geometry = gdf.geometry.make_valid()
    gdf.geometry = gdf.geometry.map(_cast_to_mp)
    gdf.index = gdf[cell_id_col].astype(str)
    # spatialdata promotes 2D polygons to 3D if any column is literally named "z"
    if "z" in gdf.columns:
        gdf = gdf.rename(columns={"z": "z_layer"})
    return ShapesModel.parse(gdf, transformations=transformations)


def get_table(
    count_path: str | Path,
    obs_path: str | Path,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    region: str,
    *,
    cell_id_col: str = _CELL_KEY,
) -> ad.AnnData:
    """Build the annotated table from a cell-by-gene CSV + an obs (metadata) CSV,
    keyed on ``cell_id_col``. Blank barcodes are moved to ``obsm['blank']``.
    (pandas infers ``.gz`` compression from the filename.)"""
    data = pd.read_csv(count_path, index_col=0)
    obs = pd.read_csv(obs_path, index_col=0)
    obs.index = obs.index.astype(str)
    data.index = data.index.astype(str)

    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = ad.AnnData(data.loc[:, is_gene], dtype=data.values.dtype, obs=obs)
    adata.obsm["blank"] = data.loc[:, ~is_gene]
    if {"center_x", "center_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
    elif {"centroid_x", "centroid_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["centroid_x", "centroid_y"]].values

    return _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col)


def get_table_from_adata(
    adata: ad.AnnData,
    obs_df: pd.DataFrame | None,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    region: str,
    *,
    cell_id_col: str = _CELL_KEY,
) -> ad.AnnData:
    """Build the table from a counts AnnData (e.g. proseg's zarr ``table``),
    optionally overriding obs with ``obs_df`` (the segmentation schema ``cell_metadata.csv``).
    Blank genes go to ``obsm['blank']``; keyed on ``cell_id_col``."""
    adata = adata.copy()
    if obs_df is not None:
        obs_df = obs_df.copy()
        obs_df.index = obs_df.index.astype(str)
        adata.obs = obs_df

    gene = adata.var["gene"] if "gene" in adata.var.columns else pd.Series(adata.var_names, index=adata.var_names)
    is_blank = gene.str.startswith("Blank-").to_numpy()
    blank = adata[:, is_blank].copy()
    adata = adata[:, ~is_blank].copy()
    adata.obs["transcript_count"] = adata.X.toarray().sum(axis=1).flatten()
    adata.obsm["blank"] = blank.X.toarray()
    if {"centroid_x", "centroid_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["centroid_x", "centroid_y"]].values
    adata.uns.pop("spatialdata_attrs", None)  # re-parse from scratch

    return _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col)


def _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col):
    """Attach the standard obs bookkeeping columns and parse into a TableModel."""
    adata.obs["region"] = pd.Series(reg_name, index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(exp_name, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(dataset_id, index=adata.obs_names, dtype="category")
    adata.obs[REGION_KEY] = pd.Series(region, index=adata.obs_names, dtype="category")
    adata.obs[cell_id_col] = adata.obs.index
    return TableModel.parse(
        adata, region_key=REGION_KEY, region=region, instance_key=cell_id_col
    )


def parse_points(transcripts_df, transformations):
    """Parse a proseg transcript table (zarr ``transcripts`` or CSV) into a
    PointsModel using the proseg column preset."""
    gene_col = PROSEG_PRESET["tz_col"]["gene_col"]
    transcripts_df[gene_col] = transcripts_df[gene_col].astype("string")
    return PointsModel.parse(
        transcripts_df,
        coordinates={
            "x": PROSEG_PRESET["tz_col"]["x_col"],
            "y": PROSEG_PRESET["tz_col"]["y_col"],
        },
        transformations=transformations,
        feature_key=gene_col,
    )
