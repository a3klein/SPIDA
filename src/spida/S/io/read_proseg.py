import gzip
from pathlib import Path
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from spatialdata.models import PointsModel, ShapesModel, TableModel  # type: ignore
    from spatialdata.transformations import BaseTransformation

import dask.dataframe as dd
import geopandas as gpd
from shapely.geometry import MultiPolygon
import anndata as ad
import pandas as pd

from spida._constants import PROSEG_PRESET, REGION_KEY, GEOMETRY_COL


def _get_polygons(
    boundaries_path: Path, transformations: dict[str, BaseTransformation]
) -> gpd.GeoDataFrame:
    """
    Get the polygons from the spatialdata object.
    """
    # Reading the polygons file (from boundary path) into a GeoDataFrame.
    if boundaries_path.endswith(".gz"):
        with gzip.open(boundaries_path, "rt") as f:
            cell_polygons = gpd.read_file(f)
    elif boundaries_path.endswith(".parquet"):
        cell_polygons = gpd.read_parquet(boundaries_path)
    else:
        cell_polygons = gpd.read_file(boundaries_path)

    if cell_polygons.geometry.name != GEOMETRY_COL:
        cell_polygons = cell_polygons.rename_geometry(GEOMETRY_COL)
    cell_polygons = cell_polygons[cell_polygons.geometry.is_valid]
    cell_polygons.geometry = cell_polygons.geometry.map(lambda x: MultiPolygon(x.geoms))
    cell_polygons.index = cell_polygons["cell"].astype(str)

    shapes = ShapesModel.parse(cell_polygons, transformations=transformations)
    return shapes


def _get_points(
    transcripts_path: Path, transformations: dict[str, BaseTransformation]
) -> dd.DataFrame:
    """
    Get the points from proseg output for the spatialdata object.
    """
    transcripts_df = dd.read_csv(transcripts_path, compression="gzip")
    transcripts = PointsModel.parse(
        transcripts_df,
        coordinates={
            "x": PROSEG_PRESET["tz_col"]["x_col"],
            "y": PROSEG_PRESET["tz_col"]["y_col"],
        },
        transformations=transformations,
        feature_key=PROSEG_PRESET["tz_col"]["gene_col"],
    )
    return transcripts


def _get_table(
    counts_path: str = None,
    obs_path: str = None,
    reg_name: str = None,
    exp_name: str = None,
    dataset_id: str = None,
    shapes_key: str = None,
):
    """
    Get the table from the spatialdata object.
    """

    ### Get Table Proseg
    data = pd.read_csv(counts_path, compression="gzip")
    obs = pd.read_csv(obs_path, compression="gzip")

    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = ad.AnnData(data.loc[:, is_gene], obs=obs, dtype=data.values.dtype)

    adata.obsm["blank"] = data.loc[:, ~is_gene].values
    adata.obsm["spatial"] = adata.obs[
        [PROSEG_PRESET["meta_cols"]["x_col"], PROSEG_PRESET["meta_cols"]["y_col"]]
    ].values
    adata.obs["region"] = pd.Series(reg_name, index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(exp_name, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(
        dataset_id, index=adata.obs_names, dtype="category"
    )
    adata.obs[REGION_KEY] = pd.Series(
        shapes_key, index=adata.obs_names, dtype="category"
    )
    adata.obs["cell"] = adata.obs.index

    table = TableModel.parse(
        adata, region_key=REGION_KEY, region=shapes_key, instance_key="cell"
    )
    return table
