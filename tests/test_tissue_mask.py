"""Unit tests for the KDE tissue mask (spida.P.tissue_mask), synthetic fixtures only."""
import numpy as np
import pytest

from spida.P.tissue_mask import (
    bin_transcripts,
    tissue_mask,
    mask_to_polygons,
    rasterize_pixel_mask,
    GB,
    REF_UM,
)


def _block_grid(n=100, lo=20, hi=80, count=20.0):
    tot = np.zeros((n, n))
    tot[lo:hi, lo:hi] = count
    return tot


def test_bin_transcripts_counts_and_origin():
    # points at known micron coords; gb=5 -> cell index = floor((v-min)/5)
    x = np.array([0.0, 0.0, 12.0, 100.0])
    y = np.array([0.0, 0.0, 3.0, 50.0])
    tot, x0, y0 = bin_transcripts(x, y, gb=5.0)
    assert (x0, y0) == (0.0, 0.0)
    assert tot[0, 0] == 2  # two points at origin cell
    assert tot[0, 2] == 1  # (x=12 -> col 2, y=3 -> row 0)
    assert tot.sum() == 4


def test_tissue_mask_dense_vs_sparse():
    tot = _block_grid()
    final = tissue_mask(tot)["final"]
    assert final[50, 50]                 # inside the dense block -> tissue
    assert not final[0, 0]               # sparse corner -> background
    # smoothing spreads the footprint beyond the raw block, erosion trims it back;
    # tissue is present but not the whole grid
    assert 3000 < final.sum() < tot.size


def test_tissue_mask_drops_small_components():
    # large block survives a_min; a FAR-SEPARATED tiny block is dropped (no halo bridge).
    tot = np.zeros((150, 150))
    tot[20:80, 20:80] = 20.0     # 60x60 -> footprint > a_min/gb^2 (=2000)
    tot[120:130, 120:130] = 20.0  # 10x10, isolated -> footprint << a_min -> dropped
    final = tissue_mask(tot, erode_sigma=0.0)["final"]
    assert final[50, 50]                          # large block kept
    assert not final[115:135, 115:135].any()      # small isolated component removed


def test_erosion_shrinks_boundary():
    tot = _block_grid()
    none = tissue_mask(tot, erode_sigma=0.0)["final"]
    one = tissue_mask(tot, erode_sigma=1.0)["final"]
    assert one.sum() < none.sum()
    # 1 sigma erosion ~ round(REF_UM/GB) cells removed from each side
    n_iter = round(REF_UM / GB)
    assert none.sum() > 0 and one.sum() > 0
    # eroded bounding box shrinks by ~n_iter on each side (allow slack for closing)
    rows_none = np.where(none.any(axis=1))[0]
    rows_one = np.where(one.any(axis=1))[0]
    assert rows_one[0] >= rows_none[0] + n_iter - 1
    assert rows_one[-1] <= rows_none[-1] - n_iter + 1


def test_mask_to_polygons_area_and_filtered():
    mask = np.zeros((100, 100), dtype=bool)
    mask[10:30, 10:40] = True     # 20 rows x 30 cols
    grid = mask_to_polygons(mask, x0=0.0, y0=0.0, gb=5.0)
    assert list(grid["filtered"].unique()) == [False]     # pass region flagged False
    expected_area = 20 * 30 * (5.0 ** 2)                   # px count x gb^2
    assert grid.geometry.area.sum() == pytest.approx(expected_area, rel=1e-6)


def test_rasterize_pixel_mask_roundtrip():
    mask = np.zeros((100, 100), dtype=bool)
    mask[10:30, 10:40] = True
    grid = mask_to_polygons(mask, x0=0.0, y0=0.0, gb=5.0)   # micron rect [50:200]x, [50:150]y
    m2m = np.eye(3)                                         # micron == pixel
    shape = (500, 500)
    px = rasterize_pixel_mask(grid, m2m, shape)
    assert px.shape == shape
    expected_frac = (150 * 100) / (500 * 500)              # micron area / total
    assert px.mean() == pytest.approx(expected_frac, abs=2e-3)


def test_rasterize_empty_polygons():
    empty = mask_to_polygons(np.zeros((10, 10), dtype=bool), 0.0, 0.0, gb=5.0)
    px = rasterize_pixel_mask(empty, np.eye(3), (50, 50))
    assert px.shape == (50, 50) and not px.any()
