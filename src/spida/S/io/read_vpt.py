from pathlib import Path
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from spatialdata.models import ShapesModel, TableModel
    from spatialdata.transformations import BaseTransformation

import pandas as pd
import anndata as ad
import geopandas as gpd
from shapely.geometry import MultiPolygon

from spida._constants import DEFAULT_PRESET, GEOMETRY_COL


def cast_to_mp(x): 
    try: 
        return MultiPolygon(x.geoms)
    except Exception: 
        return MultiPolygon([x.convex_hull])

def _get_polygons(boundaries_path: Path, transformations: dict[str, BaseTransformation]) -> gpd.GeoDataFrame:
    geo_df = gpd.read_parquet(boundaries_path)
    if geo_df.geometry.name != GEOMETRY_COL:
        geo_df = geo_df.rename_geometry(GEOMETRY_COL)
    geo_df = geo_df[geo_df[DEFAULT_PRESET['geom_depth_col']] == 0]  # Avoid duplicate boundaries on all z-levels
    geo_df.geometry = geo_df.geometry.make_valid() # Ensure geometries are valid
    # geo_df = geo_df[geo_df.geometry.is_valid]  # Remove invalid geometries
    geo_df.geometry = geo_df.geometry.map(lambda x: cast_to_mp(x))
    geo_df.index = geo_df[DEFAULT_PRESET['METADATA_CELL_KEY']].astype(str)
    return ShapesModel.parse(geo_df, transformations=transformations)


# Adjusting some things from the spatialdata_io merscope default
def _get_table(
    count_path: Path,
    obs_path: Path,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    region: str,
) -> ad.AnnData:
    
    data = pd.read_csv(count_path, index_col=0, dtype={"cell": str})
    obs = pd.read_csv(obs_path, index_col=0, dtype={"EntityID": str})

    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = ad.AnnData(data.loc[:, is_gene], dtype=data.values.dtype, obs=obs)

    adata.obsm["blank"] = data.loc[:, ~is_gene]  # blank fields are excluded from adata.X
    adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
    adata.obs["region"] = pd.Series(reg_name, index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(exp_name, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(dataset_id, index=adata.obs_names, dtype="category")
    adata.obs["cells_region"] = pd.Series(region, index=adata.obs_names, dtype="category")
    adata.obs["EntityID"] = adata.obs.index

    table = TableModel.parse(
        adata,
        region_key="cells_region",
        region=region,
        instance_key="EntityID",
    )
    return table