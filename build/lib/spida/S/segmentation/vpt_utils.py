import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from pyarrow.parquet import ParquetFile
from fsspec.implementations.local import LocalFileSystem
from typing import Callable
from shapely import affinity, geometry, MultiPolygon
import datetime

counter_digits_num = 8
counter = 0
process_id = ""

# Functions pulled from VPT
def read_micron_to_mosaic_transform(path):
    fs = LocalFileSystem()
    lines = fs.open(path, "r", lambda f: f.readlines())

    def process_line(line: str): 
        return list(map(float, line.split()))

    transform = list(map(process_line, lines))

    if len(transform) != 3 or not all(map(lambda row: len(row) == 3, transform)):
        raise ValueError("Micron to mosaic transform should be a 3x3 matrix")
    
    return transform

def convert_to_multipoly(shape):
    if type(shape) is geometry.Polygon:
        return geometry.multipolygon.MultiPolygon([shape])
    elif type(shape) is geometry.MultiPolygon:
        return shape
    # elif type(shape) is geometry.GeometryCollection:
    #     poly_shapes = [g for g in shape.geoms if type(g) in [geometry.Polygon, geometry.MultiPolygon]]
    #     poly = get_valid_geometry(unary_union(poly_shapes))
    #     return poly if type(poly) is geometry.MultiPolygon else geometry.multipolygon.MultiPolygon([poly])
    else:
        # If type is not Polygon or Multipolygon, the shape
        # is strange / small and should be rejected
        return geometry.MultiPolygon()

def update_geometry(gdf: gpd.GeoDataFrame, callback: Callable, *args, **kwargs): 
    def transform(x) -> geometry.MultiPolygon:
        return convert_to_multipoly(callback(x, *args, **kwargs))

    gdf.geometry = gdf.geometry.apply(transform)

def transform_geoms(gdf: gpd.GeoDataFrame, matrix):
    tform_flat = [*matrix[:2, :2].flatten(), *matrix[:2, 2].flatten()]
    gdf = gdf.copy()
    update_geometry(gdf, affinity.affine_transform, matrix=tform_flat)
    return gdf

def format_experiment_timestamp(timestamp: float):
    time = datetime.datetime.fromtimestamp(timestamp)
    year_day = str(int(time.strftime("%j")) + 100)
    day_time = time.strftime("%H%M%S")
    day_second = str(int(day_time[:2]) * 60**2 + int(day_time[2:4]) * 60 + int(day_time[4:])).zfill(5)
    return "".join([year_day, day_second])

def set_process_id():
    global process_id, counter_digits_num
    if not process_id:
        process_id = format_experiment_timestamp(datetime.datetime.now(datetime.UTC).timestamp())
        counter_digits_num = 12 - len(process_id)

def get_id() -> np.int64:
    global counter
    next_id = np.int64(f"{process_id}{str(counter).zfill(counter_digits_num)}")
    counter += 1
    return next_id


# My function to convert the 3D geom
def _3d_convert_geometry(
    root_dir: str,
    output_dir: str,
    region: str,
    input_boundaries: str = "polygons.parquet",
    output_boundaries: str = "cellpose_micron_space.parquet",
    convert_micron: bool = True,
):
    """
    Convert geometry of cell boundaries from cellpose 3D segmentation into an output compatible with VPT. 
    This is mostly just a wrapper for VPT _convert_geometry however I am removing the overlapping geometry handling
    from VPT since it returns weird results for 3D geometry.
    """

    input_path = f"{output_dir}/{region}/{input_boundaries}"
    output_path = f"{output_dir}/{region}/{output_boundaries}"
    gdf_inp = gpd.read_parquet(input_path)
    
    if convert_micron:
        mic_to_mos_path = f"{root_dir}/{region}/images/micron_to_mosaic_pixel_transform.csv"
        micron_to_mosaic = read_micron_to_mosaic_transform(mic_to_mos_path)
        mosaic_to_micron_mat = np.linalg.inv(micron_to_mosaic)
        gdf_inp = transform_geoms(gdf_inp, mosaic_to_micron_mat)

    set_process_id()

    cell_id = 'EntityID'
    gdf_inp['ParentType'] = None
    gdf_inp['ParentID'] = None
    gdf_inp['Type'] = 'cell'
    gdf_inp['Name'] = None

    entity_ids_map = {i: get_id() for i in gdf_inp[cell_id].unique()}
    gdf_inp[cell_id] = gdf_inp[cell_id].map(entity_ids_map)

    gdf_inp.to_parquet(output_path)