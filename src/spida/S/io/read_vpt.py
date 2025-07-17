from pathlib import Path
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import BaseTransformation

import geopandas as gpd
from shapely.geometry import MultiPolygon

from spida._constants import DEFAULT_PRESET, GEOMETRY_COL


def _get_polygons(boundaries_path: Path, transformations: dict[str, BaseTransformation]) -> gpd.GeoDataFrame:
    geo_df = gpd.read_parquet(boundaries_path)
    if geo_df.geometry.name != GEOMETRY_COL:
        geo_df = geo_df.rename_geometry(GEOMETRY_COL)
    geo_df = geo_df[geo_df[DEFAULT_PRESET['geom_depth_col']] == 0]  # Avoid duplicate boundaries on all z-levels
    geo_df.geometry = geo_df.geometry.make_valid() # Ensure geometries are valid
    # geo_df = geo_df[geo_df.geometry.is_valid]  # Remove invalid geometries
    geo_df.geometry = geo_df.geometry.map(lambda x: MultiPolygon(x.geoms))
    geo_df.index = geo_df[DEFAULT_PRESET['METADATA_CELL_KEY']].astype(str)

    return ShapesModel.parse(geo_df, transformations=transformations)