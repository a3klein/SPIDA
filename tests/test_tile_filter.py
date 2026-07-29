"""Tests for the extracted tile-filter strategies (spida.S.filters.mask).

The otsu test is a GOLDEN regression: the extracted function must reproduce the old inline
decon formula (threshold_otsu(maxes)/2, keep if max>thr) exactly.
"""
import numpy as np
import pytest
import skimage as ski

from spida.S.filters.mask import (
    otsu_tile_filter,
    no_tile_filter,
    tissue_mask_tile_filter,
    select_tiles,
)


def _masked_tiles():
    """4 tiles over a 10x10 mask whose top-left 5x5 is tissue."""
    mask = np.zeros((10, 10), dtype=np.uint8)
    mask[0:5, 0:5] = 1
    tiles = [np.ones((5, 5), dtype=np.uint16) for _ in range(4)]
    tile_info = [
        {"tile_id": 0, "position": (0, 0), "end_position": (5, 5)},   # overlaps tissue -> keep
        {"tile_id": 1, "position": (5, 5), "end_position": (10, 10)},  # all background -> drop
        {"tile_id": 2, "position": (0, 5), "end_position": (5, 10)},   # all background -> drop
        {"tile_id": 3, "position": (3, 3), "end_position": (8, 8)},    # overlaps tissue -> keep
    ]
    return tiles, tile_info, mask


def _bimodal_tiles():
    # large-gap bimodal maxes so the liberal otsu/2 threshold drops the low ("background") tiles
    lows = [10, 12, 11, 9, 13, 8]
    highs = [5000, 5200, 4800, 5100, 4900, 5300]
    vals = lows + highs
    tiles = [np.full((3, 3), v, dtype=np.uint16) for v in vals]
    tile_info = [{"tile_id": i} for i in range(len(tiles))]
    return tiles, tile_info, vals


def test_otsu_tile_filter_matches_old_inline_formula():
    tiles, tile_info, vals = _bimodal_tiles()
    kept, kept_info, tile_maxes, thr = otsu_tile_filter(tiles, tile_info)

    # GOLDEN: recompute the OLD inline formula independently -> must match exactly
    old_maxes = [np.max(t) for t in tiles]
    old_thr = ski.filters.threshold_otsu(np.asarray(old_maxes)) / 2
    old_keep = [np.max(t) > old_thr for t in tiles]

    assert thr == old_thr
    assert list(tile_maxes) == old_maxes
    assert len(kept) == sum(old_keep)
    assert [int(np.max(t)) for t in kept] == [v for v, k in zip(vals, old_keep) if k]


def test_otsu_sets_tile_info_fields():
    tiles, tile_info, _ = _bimodal_tiles()
    kept, *_ = otsu_tile_filter(tiles, tile_info)
    assert all("max_intensity" in ti and "thresholded" in ti for ti in tile_info)
    assert sum(ti["thresholded"] for ti in tile_info) == len(kept)


def test_none_keeps_all_thr_zero():
    tiles, tile_info, vals = _bimodal_tiles()
    kept, kept_info, tile_maxes, thr = no_tile_filter(tiles, tile_info)
    assert len(kept) == len(tiles)
    assert thr == 0
    assert list(tile_maxes) == vals
    assert all("max_intensity" in ti for ti in tile_info)


def test_select_tiles_dispatch():
    tiles, tile_info, _ = _bimodal_tiles()
    k_otsu = select_tiles(tiles, [dict(ti) for ti in tile_info], method="otsu")
    k_ref = otsu_tile_filter(tiles, [dict(ti) for ti in tile_info])
    assert len(k_otsu[0]) == len(k_ref[0]) and k_otsu[3] == k_ref[3]

    k_none = select_tiles(tiles, [dict(ti) for ti in tile_info], method="none")
    assert len(k_none[0]) == len(tiles)

    assert select_tiles(tiles, [dict(ti) for ti in tile_info], method=None)[0]  # None -> keep all

    with pytest.raises(ValueError):
        select_tiles(tiles, tile_info, method="bogus")


def test_tissue_mask_tile_filter_keeps_overlapping_tiles():
    tiles, tile_info, mask = _masked_tiles()
    kept, kept_info, tile_maxes, thr = tissue_mask_tile_filter(tiles, tile_info, mask)
    assert thr == 0
    assert [ti["tile_id"] for ti in kept_info] == [0, 3]        # only tissue-overlapping tiles kept
    assert [ti["thresholded"] for ti in tile_info] == [True, False, False, True]
    assert all("max_intensity" in ti for ti in tile_info)       # set for plot_thr histogram


def test_select_tiles_tissue_mask_dispatch_and_requires_mask():
    tiles, tile_info, mask = _masked_tiles()
    kept, *_ = select_tiles(tiles, [dict(ti) for ti in tile_info], method="tissue_mask", mask_tif=mask)
    assert len(kept) == 2
    with pytest.raises(ValueError):
        select_tiles(tiles, tile_info, method="tissue_mask")   # missing mask_tif


def test_tissue_mask_tile_filter_reads_tif(tmp_path):
    import tifffile
    tiles, tile_info, mask = _masked_tiles()
    p = tmp_path / "tissue_mask.tif"
    tifffile.imwrite(p, mask * 255)
    kept, kept_info, _, _ = tissue_mask_tile_filter(tiles, tile_info, str(p))
    assert [ti["tile_id"] for ti in kept_info] == [0, 3]
