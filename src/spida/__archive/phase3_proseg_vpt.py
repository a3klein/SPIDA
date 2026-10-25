"""Archived in phase 3 of the IO redesign (2026-10-21): VPT-era proseg post-processing, dead since
the VPT removal. Replaced by segmentation/ingest.py `_ingest_layers` (layer polygons -> schema
boundaries; EntityID = proseg's cell id instead of a renumbering) and the native `sum_signals` + the
cell_metadata merge in `process_segmentation_region`. Not imported anywhere in SPIDA;
spida_dev/v1/partition_transcripts/validate_ingest.py still imports `_proseg_3d_to_vpt_parquet`
from its old path (left as is).
"""
# flake8: noqa
# Originally in spida/S/segmentation/backends/proseg.py.

import gzip
import logging
import pathlib

import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__package__)


def _proseg_3d_to_vpt_parquet(shapes_3d: gpd.GeoDataFrame, micron_per_z: float, output_path: str):
    shapes_3d = shapes_3d.rename_geometry("Geometry")
    shapes_3d['z'] = shapes_3d['layer']
    shapes_3d['ID'] = shapes_3d['cell']
    shapes_3d['ZIndex'] = shapes_3d['layer']
    # 1-based ZLevel (every plane owns a full micron_per_z slab), consistent with
    # cellpose/convert-geometry; column name standardized to "ZLevel".
    shapes_3d['ZLevel'] = (shapes_3d['ZIndex'] + 1) * micron_per_z
    shapes_3d['EntityID'] = pd.factorize(shapes_3d['ID'])[0] + 1
    shapes_3d.to_parquet(output_path, index=False)

def add_signals_meta_to_proseg(
    raw_dir: str | pathlib.Path,
    proseg_dir: str | pathlib.Path,
    shapes_filename: str = "cell-polygons-layers.geojson.gz",
    cell_meta_filename: str = "cell-metadata.csv.gz",
    micron_per_z: float = 1.5,
    vpt_bin_path: str | pathlib.Path | None = None,
    **kwargs
):
    """
    Add the sum signals metadata to the proseg cell metadata output.
    ``raw_dir`` is the raw MERSCOPE region dir; ``proseg_dir`` this region's proseg output dir.
    This is necessary for downstream analysis and visualization in VPT.
    """
    # from spida.S.segmentation.vpt_utils import read_micron_to_mosaic_transform, transform_geoms, update_geometry
    # import geopandas as gpd
    # import pandas as pd
    from spida.S.segmentation.vpt import _add_vpt_binary, _cli_sum_signals    
    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    if isinstance(proseg_dir, str):
        proseg_dir = pathlib.Path(proseg_dir)

    # Get paths and make sure they exist
    shapes_3d_path = proseg_dir / shapes_filename
    cell_meta_path = proseg_dir / cell_meta_filename
    assert shapes_3d_path.exists(), f"Shapes file not found at {shapes_3d_path}"
    assert cell_meta_path.exists(), f"Cell metadata file not found at {cell_meta_path}"
    # Read in the shapes and cell metadata
    with gzip.open(shapes_3d_path, "rb") as f:
        shapes_3d = gpd.read_file(f)
    logger.info(f"Read shapes from {shapes_3d_path} with shape {shapes_3d.shape}")
    with gzip.open(cell_meta_path, "rt") as f:
        cell_meta = pd.read_csv(f)
    logger.info(f"Read cell metadata from {cell_meta_path} with shape {cell_meta.shape}")

    # Convert proseg 3D shapes to vpt parquet format
    _proseg_3d_to_vpt_parquet(
        shapes_3d=shapes_3d,
        micron_per_z=micron_per_z,
        output_path= proseg_dir / "proseg_polygons.parquet"
    )
    logger.info(f"Converted proseg 3D shapes to vpt parquet format at {proseg_dir / 'proseg_polygons.parquet'}")
    
    # Sum signals 
    signals = _cli_sum_signals(
       raw_dir = raw_dir,
       seg_dir = proseg_dir,
       input_boundaries = "proseg_polygons.parquet",
       output_signals = "sum_signals.csv",
    )

    # Merge the signals with the cell metadata and save
    signals.index.name = 'EntityID'
    if (cell_meta.index.astype(int) == 0).any(): 
        signals.index = signals.index.astype(int) - 1
    cell_meta = cell_meta.merge(signals, left_on="cell", right_index=True).copy()
    cell_meta.to_csv(proseg_dir / "cell_metadata_with_signals.csv", index=False)
    logger.info(f"Merged signals with cell metadata and saved to {proseg_dir / 'cell_metadata_with_signals.csv'}")
