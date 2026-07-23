"""Tests for the segmentation ingest normalizer (segmentation/ingest.py).

Covers the three dispatch branches on synthetic inputs: cellpose 2D (delegates to
convert_geometry), cellpose 3D (pixel -> micron + ZIndex/ZLevel/EntityID), and
proseg (micron geojson layers). Asserts the segmentation schema, 1-based ZLevel, and the
per-method EntityID contract.
"""

import gzip
import shutil
from pathlib import Path

import numpy as np
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.ingest import ingest_polygons
from spida.S.segmentation.backends import get_spec

_SCHEMA_COLS = ["ID", "Geometry", "EntityID", "ParentID", "ParentType", "Type",
             "ZIndex", "ZLevel"]
_IDENTITY = "1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n"


def _identity_m2m(tmp_path):
    p = tmp_path / "micron_to_mosaic_pixel_transform.csv"
    p.write_text(_IDENTITY)
    return p


def test_ingest_cellpose_2d(tmp_path):
    # single z-plane pixel polygons -> convert_geometry (replicate to 7 planes)
    polys = gpd.GeoDataFrame(
        {"ID": [0, 1], "z": [0, 0],
         "Geometry": [geom.box(10, 10, 30, 30), geom.box(50, 50, 80, 80)]},
        geometry="Geometry",
    )
    raw = tmp_path / "polygons.parquet"
    polys.to_parquet(raw)
    out = ingest_polygons(get_spec("cellpose"), raw, tmp_path / "b.parquet",
                          micron_to_mosaic=_identity_m2m(tmp_path))
    assert list(out.columns) == _SCHEMA_COLS
    assert sorted(out["ZIndex"].unique()) == list(range(7))
    assert np.allclose(sorted(out["ZLevel"].unique()),
                       [1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5])  # 1-based
    assert out["EntityID"].nunique() == 2                      # 2 cells


def test_ingest_cellpose_3d(tmp_path):
    # multi-z pixel polygons; a cell present on z=0 and z=1
    polys = gpd.GeoDataFrame(
        {"ID": [5, 5, 9], "z": [0, 1, 0],
         "Geometry": [geom.box(10, 10, 30, 30), geom.box(11, 11, 31, 31),
                      geom.box(50, 50, 70, 70)]},
        geometry="Geometry",
    )
    raw = tmp_path / "polygons.parquet"
    polys.to_parquet(raw)
    out = ingest_polygons(get_spec("cellpose"), raw, tmp_path / "b.parquet",
                          micron_to_mosaic=_identity_m2m(tmp_path))
    assert list(out.columns) == _SCHEMA_COLS
    assert sorted(out["ZIndex"].unique()) == [0, 1]
    # ZLevel 1-based: z0 -> 1.5, z1 -> 3.0
    z0 = out[out["ZIndex"] == 0]["ZLevel"].unique()
    z1 = out[out["ZIndex"] == 1]["ZLevel"].unique()
    assert np.allclose(z0, [1.5]) and np.allclose(z1, [3.0])
    # cell 5 spans two planes but must have ONE EntityID; total 2 cells
    assert out["EntityID"].nunique() == 2


def _write_proseg_geojson_gz(tmp_path):
    layers = gpd.GeoDataFrame(
        {"cell": [1, 1, 2], "layer": [0, 1, 0],
         "geometry": [geom.box(100, 100, 140, 140), geom.box(101, 101, 141, 141),
                      geom.box(200, 200, 230, 230)]},
        geometry="geometry", crs=None,
    )
    gj = tmp_path / "cell-polygons-layers.geojson"
    layers.to_file(gj, driver="GeoJSON")
    gz = tmp_path / "cell-polygons-layers.geojson.gz"
    with open(gj, "rb") as fi, gzip.open(gz, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return gz


def test_ingest_proseg(tmp_path):
    gz = _write_proseg_geojson_gz(tmp_path)
    out = ingest_polygons(get_spec("proseg", "3"), gz, tmp_path / "b.parquet")
    assert list(out.columns) == _SCHEMA_COLS
    assert sorted(out["ZIndex"].unique()) == [0, 1]
    # 1-based ZLevel, no coordinate transform (already micron): geometry stays ~micron
    assert np.allclose(sorted(out["ZLevel"].unique()), [1.5, 3.0])
    assert out["EntityID"].nunique() == 2      # 2 proseg cells
    # cell 1 (2 planes) -> 2 rows same EntityID; cell 2 -> 1 row
    eids = out.groupby("EntityID").size().sort_values().tolist()
    assert eids == [1, 2]
    # micron space preserved (box near 100-230, not scaled to pixels)
    assert out.geometry.total_bounds[2] < 300
