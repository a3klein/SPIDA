"""Tests for the pure-Python VPT `convert-geometry` reimplementation.

Fixtures (built with the real vpt binary on a 222-cell tile-seam window that
contains overlaps):
  convert_polygons.parquet          - raw cellpose input (mosaic-pixel, 2D)
  convert_micron_to_mosaic.csv       - 3x3 transform
  expected_convert_micron_space.parquet - VPT output (micron-space 3D)

EntityIDs are timestamp-based (non-deterministic), so equivalence is asserted on
geometry: cells are centroid-matched and compared by symmetric-difference area.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.segmentation_utils import (
    convert_geometry,
    harmonize_polygons,
)

DATA = Path(__file__).parent / "data" / "segmentation_utils"
RAW = DATA / "convert_polygons.parquet"
M2M = DATA / "convert_micron_to_mosaic.csv"
EXPECTED = DATA / "expected_convert_micron_space.parquet"

pytestmark = pytest.mark.skipif(
    not RAW.exists(), reason="convert-geometry fixture not present"
)


def _z0(gdf):
    g = gdf[gdf["ZIndex"] == 0].reset_index(drop=True).copy()
    g["geometry"] = g.geometry.make_valid()
    return g


def test_convert_matches_vpt_geometry(tmp_path):
    out = convert_geometry(RAW, tmp_path / "out.parquet", M2M)
    expected = gpd.read_parquet(EXPECTED)

    # z-planes and ZLevels
    assert sorted(out["ZIndex"].unique()) == list(range(7))
    assert np.allclose(sorted(out["ZLevel"].unique()),
                       [1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5])
    assert list(out.columns) == [
        "ID", "Geometry", "EntityID", "ParentID", "ParentType", "Type",
        "ZIndex", "ZLevel",
    ]

    m, v = _z0(out), _z0(expected)
    assert len(m) == len(v)  # same cell count after harmonize

    from scipy.spatial import cKDTree
    mc = np.c_[m.geometry.centroid.x, m.geometry.centroid.y]
    vc = np.c_[v.geometry.centroid.x, v.geometry.centroid.y]
    dist, idx = cKDTree(mc).query(vc, k=1)
    assert (dist < 1e-3).all()  # every VPT cell has a matching mine cell

    # Tolerance policy: the vast majority of cells match VPT to float precision;
    # a small fraction of contested tile-seam slivers may differ (order-sensitive,
    # see convert_geometry docstring). Assert median ~0 and only a tiny fraction >1%.
    mg, vg = m.geometry.values, v.geometry.values
    ratios = np.array([
        vg[vi].symmetric_difference(mg[idx[vi]]).area / max(vg[vi].area, 1e-9)
        for vi in range(len(v))
    ])
    assert np.median(ratios) < 1e-9              # typical cell is bit-exact
    assert (ratios > 0.01).mean() < 0.01         # <1% of cells differ by >1%


def test_convert_resolves_overlaps(tmp_path):
    """The raw input has overlapping cells; the output must have none."""
    out = convert_geometry(RAW, tmp_path / "out.parquet", M2M)
    z0 = _z0(out)
    z0["i"] = np.arange(len(z0))
    sj = gpd.sjoin(z0, z0, predicate="intersects", how="inner")
    sj = sj[sj["i_left"] < sj["i_right"]]
    g = z0.geometry.values
    inter = gpd.GeoSeries(g[sj["i_left"].to_numpy()]).intersection(
        gpd.GeoSeries(g[sj["i_right"].to_numpy()]), align=False).area.to_numpy()
    assert int((inter > 1e-6).sum()) == 0


def test_output_written_to_disk(tmp_path):
    p = tmp_path / "boundaries.parquet"
    convert_geometry(RAW, p, M2M)
    assert p.exists()
    on_disk = gpd.read_parquet(p)
    assert on_disk["EntityID"].dtype == np.int64
    assert on_disk["EntityID"].nunique() == len(on_disk) // 7


# --------------------------------------------------------------------------- #
# synthetic harmonize checks                                                  #
# --------------------------------------------------------------------------- #
def _sq(x0, y0, s=1.0):
    return geom.MultiPolygon([geom.box(x0, y0, x0 + s, y0 + s)])


def test_harmonize_trims_larger_cell():
    """Two overlapping cells (<50%): larger is trimmed, smaller kept intact."""
    g = gpd.GeoDataFrame(
        {"EntityID": [1, 2],
         "Geometry": [_sq(0, 0, 3.0), _sq(2.5, 0, 1.0)]},  # big overlaps small
        geometry="Geometry",
    )
    before = g.geometry.area.tolist()
    out = harmonize_polygons(g, min_distance=0, min_area=0)
    assert len(out) == 2
    a = dict(zip(out["EntityID"], out.geometry.area))
    assert a[2] == pytest.approx(before[1])   # small cell unchanged
    assert a[1] < before[0]                    # large cell trimmed


def test_harmonize_merges_large_overlap():
    """Overlap > 50% of the smaller cell -> smaller merged into larger, dropped."""
    g = gpd.GeoDataFrame(
        {"EntityID": [1, 2],
         "Geometry": [_sq(0, 0, 3.0), _sq(0.1, 0.1, 1.0)]},  # small ~inside big
        geometry="Geometry",
    )
    out = harmonize_polygons(g, min_distance=0, min_area=0)
    assert list(out["EntityID"]) == [1]  # smaller (2) deprecated


def test_harmonize_size_filter():
    """Cells with area <= min_area are dropped."""
    g = gpd.GeoDataFrame(
        {"EntityID": [1, 2], "Geometry": [_sq(0, 0, 10.0), _sq(50, 50, 0.5)]},
        geometry="Geometry",
    )
    out = harmonize_polygons(g, min_distance=2, min_area=1)
    assert list(out["EntityID"]) == [1]  # tiny far cell (area 0.25) filtered
