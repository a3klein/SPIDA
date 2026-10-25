"""
Native output of external segmentation methods, one section per method. A method's spec
(``S/segmentation/backends/registry.py``) names its section in ``native_format``; generic code gets
a reader with ``spec.native_reader(what)`` -> ``read_<native_format>_<what>``.
Readers return plain objects; normalizing them to the segmentation schema is S/segmentation's job.
"""

from __future__ import annotations

import gzip
from pathlib import Path


# cellpose (also mesmer)

def read_cellpose_boundaries(
    path: str | Path,
):
    """
    Pixel-space polygons from ``polygons.parquet`` (columns ID, Geometry, and z for 3D).
    Parameters:
    path (str | Path): The method's boundaries file.
    """
    import geopandas as gpd

    return gpd.read_parquet(path)


# proseg

def read_proseg_boundaries(
    path: str | Path,
):
    """
    Micron-space per-layer polygons from ``cell-polygons-layers.geojson[.gz]``.
    Parameters:
    path (str | Path): The method's boundaries file.
    """
    import geopandas as gpd

    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return gpd.read_file(f)
    return gpd.read_file(path)


def read_proseg_metadata(
    path: str | Path,
    spec,
):
    """
    proseg's per-cell metadata, indexed by EntityID (proseg's cell id).
    Parameters:
    path (str | Path): The method's metadata file (``cell-metadata.csv.gz``).
    spec (SegmentationClass): The proseg spec (its cell-id column).
    """
    import pandas as pd

    meta = pd.read_csv(path).set_index(spec.columns.cell_id)
    meta.index.name = "EntityID"      # proseg cell id == segmentation schema EntityID
    return meta


def read_proseg_assigned(
    seg_dir: str | Path,
    spec,
) -> dict:
    """
    proseg's own transcripts, 2D union cell polygons and counts.
    v3 (``spec.spatialdata_file`` present): ``{"version": 3, "transcripts", "cells", "counts"}`` —
    the elements of proseg_outputs.zarr. v2: ``{"version": 2, "transcripts", "cells", "counts_file"}`` —
    the CSVs and the union geojson.
    Parameters:
    seg_dir (str | Path): The segmentation's folder.
    spec (SegmentationClass): The proseg spec (v2 / v3 file names).
    """
    seg_dir = Path(seg_dir)
    zarr_path = seg_dir / (spec.spatialdata_file or "__none__")
    if spec.spatialdata_file and zarr_path.exists():
        import spatialdata as sd

        proseg = sd.read_zarr(str(zarr_path))
        return {"version": 3, "transcripts": proseg["transcripts"], "cells": proseg["cell_boundaries"],
                "counts": proseg["table"]}

    import dask.dataframe as dd
    import geopandas as gpd

    tx = dd.read_csv(str(seg_dir / spec.transcripts_file), compression="gzip")
    with gzip.open(str(seg_dir / "cell-polygons.geojson.gz"), "rt") as f:  # 2D union (native proseg)
        cells = gpd.read_file(f)
    return {"version": 2, "transcripts": tx, "cells": cells, "counts_file": seg_dir / spec.counts_file}
