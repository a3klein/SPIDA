"""Tests for the shared mask->polygon conversion (backends/masks.py).

No cellpose/deepcell dependency, so these run in the standard preprocessing env.

``masks_to_geodataframe`` dissolves internally by (ID, z): tiling is an internal
detail, so the output is one polygon per cell per z-plane (no tile fragments leak).
"""

import numpy as np
import pytest

from spida.S.segmentation.backends import masks


def test_simple_two_cells_one_tile():
    m = np.zeros((32, 32), dtype=np.uint16)
    m[2:8, 2:8] = 1
    m[12:20, 12:20] = 2

    masks.max_cpu = 1
    gdf = masks.masks_to_geodataframe(m, tile_size=32, overlap=8)

    assert not gdf.empty
    assert set(gdf["ID"].unique()) == {1, 2}
    assert gdf.geometry.notna().all()
    # single tile, single plane -> exactly one polygon per cell
    assert len(gdf) == 2


def test_2d_cell_spanning_tiles_is_single_dissolved_polygon():
    """A cell straddling tile seams comes back as ONE polygon at the right place."""
    scale = 4
    m = np.zeros((1100, 1100), dtype=np.uint16)
    m[200:240, 220:260] = 1      # centroid col 239.5, row 219.5
    m[800:840, 900:940] = 2      # centroid col 919.5, row 819.5

    masks.max_cpu = 1
    gdf = masks.masks_to_geodataframe(m, tile_size=64, overlap=16)

    # dissolve is internal now: exactly one polygon per cell, no fragments
    assert len(gdf) == 2
    assert set(gdf["ID"]) == {1, 2}

    gdf.geometry = gdf.geometry.affine_transform([scale, 0, 0, scale, 0, 0])
    expected = {}
    for cid in (1, 2):
        rc = np.column_stack(np.where(m == cid)).mean(axis=0)   # (row, col)
        expected[cid] = np.array([rc[1] * scale, rc[0] * scale])  # (x=col, y=row)
    for _, row in gdf.iterrows():
        c = np.array([row["Geometry"].centroid.x, row["Geometry"].centroid.y])
        assert np.allclose(c, expected[row["ID"]], atol=scale)


def test_3d_preserves_z_planes():
    """A cell on multiple z-planes stays one polygon PER plane (planes not merged)."""
    m = np.zeros((3, 200, 200), dtype=np.uint16)
    m[0, 40:80, 40:80] = 1        # cell 1 on z=0
    m[1, 45:95, 45:85] = 1        # cell 1 on z=1 (different shape)
    m[0, 120:150, 120:150] = 2    # cell 2 only on z=0

    masks.max_cpu = 1
    gdf = masks.masks_to_geodataframe(m, tile_size=64, overlap=16)

    # cell 1 -> 2 rows (z=0, z=1); cell 2 -> 1 row (z=0)
    c1 = gdf[gdf["ID"] == 1]
    c2 = gdf[gdf["ID"] == 2]
    assert sorted(c1["z"].tolist()) == [0, 1]
    assert c2["z"].tolist() == [0]
    # each (cell, plane) is a single dissolved polygon
    assert len(c1) == 2 and len(c2) == 1
    assert c1.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).all()


def test_empty_mask_returns_empty():
    m = np.zeros((32, 32), dtype=np.uint16)
    masks.max_cpu = 1
    gdf = masks.masks_to_geodataframe(m, tile_size=32, overlap=8)
    assert gdf.empty
