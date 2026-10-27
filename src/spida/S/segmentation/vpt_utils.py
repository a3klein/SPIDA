import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from pyarrow.parquet import ParquetFile
from spida.io.formats.merscope import IMAGES_DIR, TRANSFORM_FILE, read_transform
from typing import Callable
from shapely import affinity, geometry, MultiPolygon
import datetime

counter_digits_num = 8
counter = 0
process_id = ""

# Functions pulled form VPT.
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
    raw_dir: str,
    seg_dir: str,
    input_boundaries: str = "polygons.parquet",
    output_boundaries: str = "cellpose_micron_space.parquet",
    convert_micron: bool = True,
    spacing_z: float = 1.5,
):
    """
    Convert geometry of cell boundaries from cellpose 3D segmentation into an output compatible with VPT. 
    This is mostly just a wrapper for VPT _convert_geometry however I am removing the overlapping geometry handling
    from VPT since it returns weird results for 3D geometry.
    """

    input_path = f"{seg_dir}/{input_boundaries}"
    output_path = f"{seg_dir}/{output_boundaries}"
    gdf_inp = gpd.read_parquet(input_path)
    
    if convert_micron:
        mic_to_mos_path = f"{raw_dir}/{IMAGES_DIR}/{TRANSFORM_FILE}"
        micron_to_mosaic = read_transform(mic_to_mos_path)
        mosaic_to_micron_mat = np.linalg.inv(micron_to_mosaic)
        gdf_inp = transform_geoms(gdf_inp, mosaic_to_micron_mat)

    set_process_id()

    # Derive the segmentation-schema fields from the raw (ID, z) cellpose output
    # (run_cellpose no longer stamps them). Matches the native ingest conventions:
    # ZIndex = z, ZLevel = (z+1)*spacing_z (1-based), one timestamp EntityID per
    # cellpose ID (label consistent across a cell's z-planes).
    cell_id = 'EntityID'
    gdf_inp['ZIndex'] = gdf_inp['z']
    gdf_inp['ZLevel'] = (gdf_inp['z'] + 1) * spacing_z
    gdf_inp['ParentType'] = None
    gdf_inp['ParentID'] = None
    gdf_inp['Type'] = 'cell'
    gdf_inp['Name'] = None

    entity_ids_map = {i: get_id() for i in gdf_inp['ID'].unique()}
    gdf_inp[cell_id] = gdf_inp['ID'].map(entity_ids_map)

    gdf_inp.to_parquet(output_path)