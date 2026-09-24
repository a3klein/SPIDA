"""Tests for the tile-selection strategies (spida.S.image_prep.tile_selection).

The otsu test is a GOLDEN regression: the extracted function must reproduce the old inline
decon formula (threshold_otsu(maxes)/2, keep if max>thr) exactly.
"""
import numpy as np
import pytest
import skimage as ski

from spida.S.image_prep.tile_selection import (
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


def _masked_tiles_3d():
    """The same 4 tiles with 3D (z, row, col) positions, as the 3D tiler emits.

    `tile_image_with_overlap` branches on `len(shape)`: 2D gives `(row, col)`, 3D gives
    `(d, row, col)`. `decon_image` passes the full (z, H, W) volume, so every 3D
    deconvolution takes the second branch. The filter used to unpack these as 2-tuples
    and raised `ValueError: too many values to unpack` on the first tile, which meant
    `--image-filter tissue_mask` had never worked for 3D decon. Keep/drop decisions must
    be identical to the 2D case, since the mask is a 2D mosaic-pixel raster and z is
    irrelevant to tissue overlap.
    """
    tiles, tile_info, mask = _masked_tiles()
    for ti in tile_info:
        ti["position"] = (0, *ti["position"])
        ti["end_position"] = (13, *ti["end_position"])
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


def test_tissue_mask_tile_filter_handles_3d_positions():
    """REGRESSION: 3D (z, row, col) positions must give the same answer as 2D.

    This is the case that broke job 2467 -- and that no existing fixture covered, so
    the suite passed while the only real caller (`decon_image` on a z-stack) crashed.
    """
    tiles, tile_info, mask = _masked_tiles_3d()
    kept, kept_info, tile_maxes, thr = tissue_mask_tile_filter(tiles, tile_info, mask)
    assert thr == 0
    assert [ti["tile_id"] for ti in kept_info] == [0, 3]
    assert [ti["thresholded"] for ti in tile_info] == [True, False, False, True]


def test_tissue_mask_2d_and_3d_agree():
    """The z component must not change which tiles are kept."""
    t2, i2, m = _masked_tiles()
    t3, i3, _ = _masked_tiles_3d()
    k2 = [ti["tile_id"] for ti in tissue_mask_tile_filter(t2, i2, m)[1]]
    k3 = [ti["tile_id"] for ti in tissue_mask_tile_filter(t3, i3, m)[1]]
    assert k2 == k3


def test_tile_rowcol_rejects_unexpected_shape():
    """A malformed position must name the shape, not mis-unpack silently."""
    from spida.S.image_prep.tile_selection import _tile_rowcol
    with pytest.raises(ValueError, match="2 .* or 3 .* elements"):
        _tile_rowcol({"position": (0, 1, 2, 3), "end_position": (4, 5, 6, 7)})


def test_tissue_mask_tile_filter_reads_tif(tmp_path):
    import tifffile
    tiles, tile_info, mask = _masked_tiles()
    p = tmp_path / "tissue_mask.tif"
    tifffile.imwrite(p, mask * 255)
    kept, kept_info, _, _ = tissue_mask_tile_filter(tiles, tile_info, str(p))
    assert [ti["tile_id"] for ti in kept_info] == [0, 3]
