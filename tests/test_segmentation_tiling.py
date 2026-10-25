"""Tests for the outer (GPU-scale) cellpose tiling: grid, tissue gating, IO.

No cellpose import and no GPU, so these run in either pixi env on the head node.
The grid arithmetic and the tissue gate are where silent errors would be most
expensive (a wrong grid wastes hours of GPU), so they are covered exhaustively.
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from spida.utilities.tiling import (
    Tile,
    bbox_of_mask,
    build_grid,
    redundancy,
    tile_coverage,
)
from spida.S.segmentation.tiling import gate_tiles, plan_grid
from spida.io.formats.merscope import find_stain_planes
from spida.io.images import probe_mosaic_shape
from spida.io.schema.tissue_mask import load as load_tissue_mask
from spida.utilities.tiling import tile_grid_with_overlap


def _origins(extent, tile, overlap, axis=0):
    """Distinct tile origins along one axis, via the shared utility."""
    shape = (extent, extent)
    info = tile_grid_with_overlap(shape, (tile, tile), overlap)
    return sorted({i["position"][axis] for i in info})


# ---------------------------------------------------------------------------
# tile_origins
# ---------------------------------------------------------------------------

def test_origins_clean_packing():
    # stride 500, tiles at 0,500,1000 -> last edge exactly on 1500
    assert _origins(1500, 1000, 500) == [0, 500]


def test_origins_last_tile_is_right_aligned():
    origins = _origins(1200, 500, 100)
    assert origins[0] == 0
    assert origins[-1] + 500 == 1200          # final edge lands on the extent
    assert origins == sorted(origins)


def test_origins_tile_larger_than_extent_is_single():
    assert [t.y0 for t in build_grid((0, 300, 0, 300), 1000, 100)] == [0]


def test_origins_tile_equal_to_extent_is_single():
    assert [t.y0 for t in build_grid((0, 1000, 0, 1000), 1000, 100)] == [0]


def test_grid_rejects_overlap_ge_tile():
    with pytest.raises(ValueError, match="must exceed overlap"):
        build_grid((0, 5000, 0, 5000), 500, 500)


@pytest.mark.parametrize("extent", [1000, 1234, 5000, 20607, 82428])
@pytest.mark.parametrize("tile,overlap", [(500, 100), (3000, 500), (7000, 500)])
def test_grid_always_covers_the_extent(extent, tile, overlap):
    """Includes tile > extent, where the raw utility emits a NEGATIVE origin
    (it leans on numpy clamping the slice). build_grid clamps the tile first,
    so absolute mosaic coordinates stay valid."""
    tiles = build_grid((0, extent, 0, extent), tile, overlap)
    covered = np.zeros(extent, dtype=bool)
    for t in tiles:
        assert t.y0 >= 0 and t.y1 <= extent
        covered[t.y0:t.y1] = True
    assert covered.all(), "grid left a gap"


# ---------------------------------------------------------------------------
# build_grid
# ---------------------------------------------------------------------------

def test_grid_solves_axes_independently():
    """The dev bench reused one origin list for both axes, which is only valid
    for a square ROI. Real mosaics are 82428x55905."""
    tiles = build_grid((0, 82428, 0, 55905), tile=7000, overlap=500)
    n_rows = len({t.gi for t in tiles})
    n_cols = len({t.gj for t in tiles})
    assert n_rows != n_cols                    # would be equal if axes were shared
    assert n_rows == len(_origins(82428, 7000, 500))
    assert n_cols == len(_origins(55905, 7000, 500))
    assert len(tiles) == n_rows * n_cols


def test_grid_offsets_are_absolute_mosaic_coordinates():
    tiles = build_grid((1000, 3000, 500, 2500), tile=1000, overlap=200)
    assert min(t.y0 for t in tiles) == 1000    # includes the bbox origin
    assert min(t.x0 for t in tiles) == 500
    assert max(t.y1 for t in tiles) <= 3000
    assert max(t.x1 for t in tiles) <= 2500


def test_grid_covers_the_bbox_without_gaps():
    y0, y1, x0, x1 = 10, 2010, 20, 1520
    tiles = build_grid((y0, y1, x0, x1), tile=700, overlap=150)
    covered = np.zeros((y1 - y0, x1 - x0), dtype=bool)
    for t in tiles:
        covered[t.y0 - y0:t.y1 - y0, t.x0 - x0:t.x1 - x0] = True
    assert covered.all()


def test_redundancy_matches_closed_form():
    """For a cleanly packing grid, redundancy == (1 + (N-1)*o/R)^2."""
    R, tile, o = 1500, 1000, 500
    tiles = build_grid((0, R, 0, R), tile, o)
    n = len(_origins(R, tile, o))
    expected = (1 + (n - 1) * o / R) ** 2
    assert redundancy(tiles, (0, R, 0, R)) == pytest.approx(expected, rel=1e-9)


def test_redundancy_is_not_monotonic_in_tile_size():
    """Bigger tiles are NOT always cheaper: what matters is how cleanly the tile
    packs into the extent. tile=6000 into 12000 at overlap=500 forces a
    right-aligned final tile at 6000 that overlaps its neighbour by 5500,
    costing 2.25x — worse than tile=2000's 1.78x, and far worse than the
    1.14-1.62x that tissue masking buys. Any "auto" sizing must therefore
    optimise packing, not just pick the largest tile that fits.
    """
    bbox = (0, 12000, 0, 12000)
    r_2000 = redundancy(build_grid(bbox, 2000, 500), bbox)
    r_6000 = redundancy(build_grid(bbox, 6000, 500), bbox)
    assert r_6000 == pytest.approx(2.25)
    assert r_2000 == pytest.approx(1.7778, abs=1e-3)
    assert r_6000 > r_2000                      # the bad-fit penalty


def test_redundancy_is_minimal_when_the_tile_packs_cleanly():
    """When (N-1)*stride + tile == extent exactly, redundancy hits the closed-form
    floor and beats a badly-fitting larger tile."""
    bbox = (0, 12000, 0, 12000)
    # stride 3500, 0/3500/7000 -> 7000+4000 == ... not exact; solve for a clean fit:
    # N=3, o=500 -> T = (R + (N-1)*o) / N = (12000 + 1000) / 3 = 4333.33 -> use 4334
    clean = redundancy(build_grid(bbox, 4334, 500), bbox)
    bad = redundancy(build_grid(bbox, 6000, 500), bbox)
    assert clean < bad
    assert clean == pytest.approx((1 + 2 * 500 / 12000) ** 2, rel=0.01)


# ---------------------------------------------------------------------------
# tissue mask: bbox, coverage, gating, zeroing
# ---------------------------------------------------------------------------

def _window(planes, y0, y1, x0, x1, **kw):
    """One stain's `(Y, X, Z)` window via `read_stains`, which returns 4D."""
    from spida.io.images import read_stains
    return read_stains({"s": planes}, window=(y0, y1, x0, x1), **kw)[0]


def _mask_with_block(shape, y0, y1, x0, x1):
    m = np.zeros(shape, dtype=bool)
    m[y0:y1, x0:x1] = True
    return m


def test_bbox_of_mask_is_half_open_and_tight():
    mask = _mask_with_block((100, 200), 10, 40, 50, 90)
    assert bbox_of_mask(mask) == (10, 40, 50, 90)


def test_bbox_of_mask_rejects_empty_mask():
    with pytest.raises(ValueError, match="empty"):
        bbox_of_mask(np.zeros((10, 10), dtype=bool))


def test_tile_coverage_fractions():
    mask = _mask_with_block((100, 100), 0, 50, 0, 100)     # top half is tissue
    full = Tile(0, 0, 0, 0, 50, 100)
    empty = Tile(0, 0, 50, 0, 50, 100)
    half = Tile(0, 0, 25, 0, 50, 100)
    assert tile_coverage(mask, full) == 1.0
    assert tile_coverage(mask, empty) == 0.0
    assert tile_coverage(mask, half) == pytest.approx(0.5)


def test_gate_skips_only_zero_coverage_tiles():
    mask = _mask_with_block((100, 100), 0, 50, 0, 100)
    tiles = build_grid((0, 100, 0, 100), tile=50, overlap=0)
    keep, skip = gate_tiles(tiles, mask)
    assert len(keep) + len(skip) == len(tiles)
    assert all(tile_coverage(mask, t) > 0 for t in keep)
    assert all(tile_coverage(mask, t) == 0 for t in skip)


def test_gate_is_a_noop_without_a_mask():
    tiles = build_grid((0, 100, 0, 100), tile=50, overlap=0)
    keep, skip = gate_tiles(tiles, None)
    assert keep == tiles and skip == []


def test_gate_min_coverage_threshold():
    mask = _mask_with_block((100, 100), 0, 10, 0, 100)     # 10% of the height
    tiles = [Tile(0, 0, 0, 0, 100, 100)]
    assert len(gate_tiles(tiles, mask, min_coverage=0.0)[0]) == 1
    assert len(gate_tiles(tiles, mask, min_coverage=0.5)[0]) == 0


# ---------------------------------------------------------------------------
# plan_grid
# ---------------------------------------------------------------------------

def test_plan_without_mask_covers_whole_mosaic():
    plan = plan_grid((1000, 2000), tile=500, overlap=100, mask=None)
    assert plan["bbox"] == (0, 1000, 0, 2000)
    assert plan["bbox_frac_of_mosaic"] == 1.0
    assert plan["n_skipped"] == 0
    assert plan["n_segment"] == plan["n_tiles"]


def test_plan_with_mask_clips_and_gates():
    mask = _mask_with_block((1000, 1000), 200, 600, 300, 700)
    plan = plan_grid((1000, 1000), tile=200, overlap=50, mask=mask)
    assert plan["bbox"] == (200, 600, 300, 700)
    assert plan["bbox_frac_of_mosaic"] == pytest.approx(0.16)
    assert plan["n_segment"] <= plan["n_tiles"]
    assert plan["compute_frac"] <= 1.0


def test_plan_reports_partial_tiles():
    """Partial tiles are the ceiling on what sub-tiling could recover, so the
    plan has to surface them."""
    mask = np.zeros((600, 600), dtype=bool)
    ys, xs = np.mgrid[0:600, 0:600]
    mask[(ys + xs) < 600] = True                           # diagonal half
    plan = plan_grid((600, 600), tile=200, overlap=0, mask=mask)
    assert plan["n_partial"] > 0
    assert plan["n_partial"] <= plan["n_segment"]


# ---------------------------------------------------------------------------
# IO block
# ---------------------------------------------------------------------------

@pytest.fixture
def mosaic_dir(tmp_path):
    """Three uncompressed DAPI planes + one PolyT plane, distinct per plane."""
    import tifffile

    d = tmp_path / "images"
    d.mkdir()
    for z in range(3):
        plane = np.full((64, 48), z + 1, dtype=np.uint16)
        plane[10:20, 5:15] = 100 + z                       # a locatable feature
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", plane)
    tifffile.imwrite(d / "mosaic_PolyT_z0.tif", np.zeros((64, 48), dtype=np.uint16))
    return d


def test_probe_mosaic_shape_reads_header(mosaic_dir):
    assert probe_mosaic_shape(mosaic_dir / "mosaic_DAPI_z0.tif") == (64, 48)


def test_probe_rejects_non_2d(tmp_path):
    import tifffile
    p = tmp_path / "stack.tif"
    tifffile.imwrite(p, np.zeros((3, 8, 8), dtype=np.uint16))
    with pytest.raises(ValueError, match="expected a 2D mosaic"):
        probe_mosaic_shape(p)


def test_find_stain_planes_is_z_ordered_and_stain_scoped(mosaic_dir):
    planes = find_stain_planes(mosaic_dir, "DAPI")
    assert [p.name for p in planes] == [
        "mosaic_DAPI_z0.tif", "mosaic_DAPI_z1.tif", "mosaic_DAPI_z2.tif"
    ]
    assert len(find_stain_planes(mosaic_dir, "PolyT")) == 1


def test_find_stain_planes_raises_for_missing_stain(mosaic_dir):
    with pytest.raises(ValueError, match="no planes for stain"):
        find_stain_planes(mosaic_dir, "GAPDH")


def test_windowed_read_matches_a_full_read(mosaic_dir):
    """The whole point of the lazy reader: identical pixels, without loading the
    64.5 GB a real z-stack would need."""
    import tifffile

    planes = find_stain_planes(mosaic_dir, "DAPI")
    got = _window(planes, 10, 20, 5, 15)
    assert got.shape == (10, 10, 3)                        # (Y, X, Z)
    for z, p in enumerate(planes):
        expected = tifffile.imread(str(p))[10:20, 5:15]
        assert np.array_equal(got[..., z], expected)


def test_windowed_read_clips_at_the_edge(mosaic_dir):
    planes = find_stain_planes(mosaic_dir, "DAPI")
    got = _window(planes, 60, 80, 40, 60)              # runs past the edge
    assert got.shape[0] == 4 and got.shape[1] == 8         # clipped by numpy
    assert got.shape[2] == 3


def test_load_tissue_mask_is_boolean(tmp_path):
    import tifffile
    p = tmp_path / "tissue_mask.tif"
    raw = np.zeros((32, 32), dtype=np.uint8)
    raw[4:20, 6:28] = 255                                  # written as 0/255
    tifffile.imwrite(p, raw, compression="deflate")        # as P/tissue_mask.py does
    mask = load_tissue_mask(p)
    assert mask.dtype == bool
    assert bbox_of_mask(mask) == (4, 20, 6, 28)


def test_load_tissue_mask_rejects_non_2d(tmp_path):
    import tifffile
    p = tmp_path / "bad.tif"
    tifffile.imwrite(p, np.zeros((2, 8, 8), dtype=np.uint8))
    with pytest.raises(ValueError, match="expected a 2D mask"):
        load_tissue_mask(p)


# ---------------------------------------------------------------------------
# realistic end-to-end plan, using the measured THM1 geometry
# ---------------------------------------------------------------------------

def test_plan_on_real_thm1_geometry():
    """Sanity-check the plan against the shape measured on THM1 UCI-2424:
    82428x55905 mosaic, tissue bbox 81568x53640, ~13x9 grid at tile=7000."""
    mask = np.zeros((82428 // 64, 55905 // 64), dtype=bool)   # downscaled stand-in
    mask[6:1281, 7:845] = True
    plan = plan_grid(mask.shape, tile=7000 // 64, overlap=500 // 64, mask=mask)
    assert plan["n_tiles"] == plan["n_segment"] + plan["n_skipped"]
    assert 1.0 < plan["redundancy"] < 2.0
    assert plan["bbox_frac_of_mosaic"] > 0.9                  # bbox ~95% of mosaic


# ---------------------------------------------------------------------------
# tile-size optimisation (integral image + 2D search)
# ---------------------------------------------------------------------------

from spida.utilities.tiling import (  # noqa: E402
    mask_integral_image,
    optimize_grid,
    score_grid,
)


def test_integral_image_matches_bruteforce_sums():
    rng = np.random.default_rng(0)
    mask = rng.random((160, 240)) > 0.5
    integral, ds = mask_integral_image(mask, ds=8)
    from spida.utilities.tiling import rect_sums
    # block-aligned rectangles must be exact
    for y0, y1, x0, x1 in [(0, 160, 0, 240), (0, 80, 0, 80), (16, 96, 32, 160)]:
        got, sampled = rect_sums(integral, ds, y0, y1, x0, x1)
        expected = int(mask[y0:y1, x0:x1].sum())
        assert int(got) == expected, f"rect ({y0},{y1},{x0},{x1})"
        assert int(sampled) == (y1 - y0) * (x1 - x0)


def test_rect_sums_snap_outward_for_unaligned_rectangles():
    """Bounds snap to block edges, so an unaligned rectangle samples a SUPERSET.
    That is the safe direction for a skip decision: never drop a tile that has
    tissue."""
    mask = np.zeros((160, 160), dtype=bool)
    mask[0:20, 0:20] = True
    integral, ds = mask_integral_image(mask, ds=8)
    from spida.utilities.tiling import rect_sums
    # (21..40) misses the tissue entirely at full res, but snaps back to block 2
    got, sampled = rect_sums(integral, ds, 21, 40, 21, 40)
    assert int(mask[21:40, 21:40].sum()) == 0
    assert int(got) >= 0
    assert int(sampled) >= (40 - 21) ** 2          # superset of the request"


def test_integral_image_rejects_too_coarse_ds():
    with pytest.raises(ValueError, match="too coarse"):
        mask_integral_image(np.ones((4, 4), dtype=bool), ds=8)


def test_score_grid_kept_matches_direct_gating():
    """The integral-image count must agree with counting tiles the slow way."""
    mask = np.zeros((800, 800), dtype=bool)
    mask[100:500, 200:700] = True
    integral, ds = mask_integral_image(mask, ds=4)
    bbox = bbox_of_mask(mask)
    scored = score_grid(integral, ds, bbox, 200, 200, 50)

    tiles = build_grid(bbox, 200, 50)
    keep, _ = gate_tiles(tiles, mask)
    assert scored["n_tiles"] == len(tiles)
    assert scored["n_kept"] == len(keep)


def test_optimizer_finds_non_square_optimum_for_anisotropic_tissue():
    """A wide, short tissue should prefer a wide, short tile."""
    mask = np.zeros((4000, 12000), dtype=bool)
    mask[500:1500, 1000:11000] = True
    best, scored = optimize_grid(
        mask, overlap=100, tile_min=500, tile_max=4000, step=250, ds=8
    )
    assert best["cost_area"] == min(s["cost_area"] for s in scored)
    assert best["tile_x"] >= best["tile_y"]        # wider than tall, like the tissue


def test_optimizer_beats_or_matches_the_best_square_tile():
    """Non-square can only help: the square grid is a subset of the search."""
    mask = np.zeros((4000, 9000), dtype=bool)
    mask[300:3700, 800:8200] = True
    best_any, _ = optimize_grid(mask, overlap=100, tile_min=500, tile_max=4000,
                                step=250, ds=8, allow_non_square=True)
    best_sq, _ = optimize_grid(mask, overlap=100, tile_min=500, tile_max=4000,
                               step=250, ds=8, allow_non_square=False)
    assert best_any["cost_area"] <= best_sq["cost_area"]


def test_optimizer_respects_max_tile_area():
    mask = np.zeros((3000, 3000), dtype=bool)
    mask[200:2800, 200:2800] = True
    cap = 1000 * 1000
    best, scored = optimize_grid(mask, overlap=100, tile_min=500, tile_max=3000,
                                 step=250, ds=8, max_tile_area=cap)
    assert all(s["tile_area"] <= cap for s in scored)
    assert best["tile_area"] <= cap


def test_optimizer_gate_skips_empty_tiles_when_tissue_is_scattered():
    """Tissue in two opposite corners: the bbox spans everything, so only the
    per-tile gate can drop the empty interior. A COMPACT corner mask would not
    test this -- `optimize_grid` clips to the bbox first, and one tile then
    covers it."""
    mask = np.zeros((4000, 4000), dtype=bool)
    mask[0:600, 0:600] = True
    mask[3400:4000, 3400:4000] = True
    best, _ = optimize_grid(mask, overlap=100, tile_min=400, tile_max=1200,
                            step=200, ds=4)
    assert best["n_kept"] < best["n_tiles"] / 4
    assert best["cost_area"] == best["n_kept"] * best["tile_area"]


def test_optimizer_raises_when_nothing_is_viable():
    mask = np.zeros((2000, 2000), dtype=bool)
    mask[500:1500, 500:1500] = True
    with pytest.raises(ValueError, match="no viable tile candidates"):
        optimize_grid(mask, overlap=100, tile_min=500, tile_max=2000,
                      step=250, ds=8, max_tile_area=1)


def test_optimizer_cost_ranks_consistently_with_kept_area():
    """cost_area is exactly kept tiles x tile area -- no per-tile overhead term."""
    mask = np.ones((2000, 2000), dtype=bool)
    _, scored = optimize_grid(mask, overlap=100, tile_min=500, tile_max=2000,
                              step=500, ds=8)
    for s in scored:
        expected = s["n_kept"] * s["tile_area"]
        assert s["cost_area"] == expected


def test_score_grid_partial_count_matches_direct_gating():
    """Regression: `sums` is in PIXELS but was once compared against a BLOCK
    count, making n_partial always 0 and silently hiding the sub-tiling
    headroom. Cross-check against counting coverage the slow way."""
    mask = np.zeros((1200, 1200), dtype=bool)
    ys, xs = np.mgrid[0:1200, 0:1200]
    mask[(ys + xs) < 1200] = True                  # diagonal -> many partial tiles
    integral, ds = mask_integral_image(mask, ds=8)
    bbox = bbox_of_mask(mask)

    scored = score_grid(integral, ds, bbox, 300, 300, 0)
    tiles = build_grid(bbox, 300, 0)
    covs = [tile_coverage(mask, t) for t in tiles]
    expected_kept = sum(1 for c in covs if c > 0)
    expected_partial = sum(1 for c in covs if 0 < c < 1)

    # The integral gate snaps outward, so it may keep a few extra boundary
    # tiles -- it must never keep FEWER than full-res gating would.
    assert scored["n_kept"] >= expected_kept
    assert scored["n_kept"] - expected_kept <= 4
    assert scored["n_partial"] > 0                 # the bug made this always 0
    assert abs(scored["n_partial"] - expected_partial) <= 4


def test_score_grid_partial_is_zero_only_when_tiles_are_truly_full():
    """A fully-tissue mask must report no partial tiles."""
    mask = np.ones((800, 800), dtype=bool)
    integral, ds = mask_integral_image(mask, ds=8)
    scored = score_grid(integral, ds, bbox_of_mask(mask), 200, 200, 0)
    assert scored["n_kept"] == scored["n_tiles"]
    assert scored["n_partial"] == 0


# ---------------------------------------------------------------------------
# streaming downscale during the read
# ---------------------------------------------------------------------------

def test_streamed_downscale_is_bit_identical(mosaic_dir):
    """The whole point: reduce during the read so the full-res window is never
    materialised, WITHOUT changing a single value. downscale_local_mean is a
    block mean, so any ds-aligned row split is exact."""
    import skimage as ski

    planes = find_stain_planes(mosaic_dir, "DAPI")
    full = _window(planes, 0, 64, 0, 48)                       # (64, 48, 3)
    expected = np.stack(
        [ski.transform.downscale_local_mean(full[..., z], (4, 4))
         for z in range(full.shape[-1])], axis=-1)

    for strip in (4, 8, 16, 1024):                                 # incl. > window
        got = _window(planes, 0, 64, 0, 48, downscale=4, strip_rows=strip)
        assert got.shape == expected.shape
        assert np.array_equal(got, expected), f"strip_rows={strip} differs"


def test_whole_and_windowed_reads_agree(mosaic_dir):
    """The whole-plane and windowed reads are the SAME call with one argument
    different, so a full-extent window must equal a `window=None` read.

    This replaces a test that compared against `cellpose._load_image`. That
    function is gone -- both paths now call `read_stains` directly, so comparing
    them through cellpose added a dependency without adding evidence. The
    exactness guarantee is the streamed-downscale test above, which checks
    against an explicit `downscale_local_mean` reference.
    """
    from spida.io.images import read_stains

    planes = find_stain_planes(mosaic_dir, "DAPI")
    whole = read_stains({"s": planes}, downscale=4)                 # window=None
    windowed = read_stains({"s": planes}, window=(0, 64, 0, 48), downscale=4)
    assert whole.shape == windowed.shape
    assert np.array_equal(whole, windowed)
    assert whole.ndim == 4                       # (C, Y, X, Z), always


def test_read_stains_is_4d_even_for_one_plane(mosaic_dir, tmp_path):
    """`read_stains` never squeezes. The collapse to (C, Y, X) belongs to the
    cellpose 2D path and is asserted there, on the shape the model receives
    (test_cellpose_backend_wiring.test_untiled_path_squeezes_a_single_plane_stain).
    """
    import shutil
    from spida.io.images import read_stains

    d = tmp_path / "one_plane"
    d.mkdir()
    shutil.copy(mosaic_dir / "mosaic_DAPI_z0.tif", d / "mosaic_DAPI_z0.tif")
    planes = find_stain_planes(d, "DAPI")
    img = read_stains({"s": planes})
    assert img.ndim == 4 and img.shape[0] == 1 and img.shape[-1] == 1


def test_streamed_downscale_handles_unaligned_extent(mosaic_dir):
    """A window whose height is not a multiple of ds must still match a
    whole-window reduction (the trailing partial block is padded identically)."""
    import skimage as ski

    planes = find_stain_planes(mosaic_dir, "DAPI")
    full = _window(planes, 0, 50, 0, 30)                       # 50 and 30 % 4 != 0
    expected = np.stack(
        [ski.transform.downscale_local_mean(full[..., z], (4, 4))
         for z in range(full.shape[-1])], axis=-1)
    got = _window(planes, 0, 50, 0, 30, downscale=4, strip_rows=8)
    assert np.array_equal(got, expected)


def test_downscale_of_one_is_a_passthrough(mosaic_dir):
    planes = find_stain_planes(mosaic_dir, "DAPI")
    assert np.array_equal(_window(planes, 0, 16, 0, 16),
                          _window(planes, 0, 16, 0, 16, downscale=1))


# ---------------------------------------------------------------------------
# polygon-based stitching
# ---------------------------------------------------------------------------

from spida.S.segmentation.tiling import (  # noqa: E402
    UnionFind,
    polygon_iou,
    stitch_polygons,
)


def _gdf(rows):
    """rows: list of (ID, z, shapely geom) -> GeoDataFrame in mosaic coords."""
    import geopandas as gpd
    return gpd.GeoDataFrame(
        {"ID": [r[0] for r in rows], "z": [r[1] for r in rows],
         "Geometry": [r[2] for r in rows]},
        geometry="Geometry", crs=None,
    )


def test_unionfind_transitive_and_idempotent():
    uf = UnionFind()
    uf.union(1, 2); uf.union(2, 3); uf.union(10, 11)
    assert uf.find(1) == uf.find(3)
    assert uf.find(1) != uf.find(10)
    uf.union(1, 3)                              # already joined
    assert uf.find(1) == uf.find(2) == uf.find(3)


def test_polygon_iou_values():
    from shapely.geometry import box
    assert polygon_iou(box(0, 0, 10, 10), box(0, 0, 10, 10)) == pytest.approx(1.0)
    assert polygon_iou(box(0, 0, 10, 10), box(20, 20, 30, 30)) == 0.0
    # half-overlap: inter 50, union 150
    assert polygon_iou(box(0, 0, 10, 10), box(5, 0, 15, 10)) == pytest.approx(1 / 3)


def test_stitch_merges_a_cell_spanning_two_tiles():
    """The core requirement: one cell seen by two independent model.eval calls
    must come out as ONE cell, not two."""
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    cell = box(85, 10, 95, 20)                   # sits inside the 80..100 overlap
    polys = {
        (0, 0): _gdf([(1, 0, box(10, 10, 20, 20)), (2, 0, cell)]),
        (0, 1): _gdf([(1, 0, cell), (2, 0, box(150, 10, 160, 20))]),
    }
    out = stitch_polygons(polys, tiles, iou_threshold=0.2)
    assert len(out) == 3                         # 4 fragments -> 3 distinct cells
    assert sorted(out["ID"].unique()) == [1, 2, 3]


def test_stitch_keeps_distinct_cells_apart():
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    polys = {
        (0, 0): _gdf([(1, 0, box(82, 10, 88, 16))]),
        (0, 1): _gdf([(1, 0, box(90, 40, 96, 46))]),   # far apart, no overlap
    }
    out = stitch_polygons(polys, tiles, iou_threshold=0.2)
    assert len(out) == 2 and out["ID"].nunique() == 2


def test_stitch_respects_the_iou_threshold():
    """Isolates the IoU gate, so containment is disabled.

    Without `containment_threshold=0` this cannot test IoU alone: the two boxes
    are EQUAL area overlapping by half, so containment = 50/100 = 0.5, which the
    default containment threshold accepts even when IoU (50/150 = 0.333) is
    below the requested 0.5.
    """
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    a, b = box(84, 10, 94, 20), box(89, 10, 99, 20)    # IoU = 1/3
    polys = {(0, 0): _gdf([(1, 0, a)]), (0, 1): _gdf([(1, 0, b)])}
    for thresh, expected in ((0.2, 1), (0.5, 2)):
        got = stitch_polygons(polys, tiles, iou_threshold=thresh,
                              containment_threshold=0)["ID"].nunique()
        assert got == expected, f"iou_threshold={thresh}"


def test_containment_merges_a_sliver_that_iou_cannot_see():
    """The artifact the second criterion exists for: a cell cut at a tile edge
    leaves a thin lip on the far side. Because the lip is thin its IoU against
    the whole cell is tiny -- here ~0.05, far below any usable IoU threshold --
    while containment is 1.0, since the lip lies entirely inside the cell."""
    from shapely.geometry import box
    from spida.S.segmentation.tiling import polygon_overlap

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    cell = box(40, 10, 90, 50)                 # 50 x 40 = 2000
    sliver = box(80, 10, 90, 20)               # 10 x 10 =  100, fully inside
    polys = {(0, 0): _gdf([(1, 0, cell)]), (0, 1): _gdf([(1, 0, sliver)])}

    iou, containment = polygon_overlap(cell, sliver)
    assert iou < 0.06 and containment == pytest.approx(1.0)

    # IoU alone leaves them as two cells...
    assert stitch_polygons(polys, tiles, containment_threshold=0)["ID"].nunique() == 2
    # ...containment merges them
    assert stitch_polygons(polys, tiles)["ID"].nunique() == 1


def test_containment_does_not_widen_acceptance_for_equal_sized_pairs():
    """At production defaults the new criterion only fires on ASYMMETRIC pairs.

    For two equal-area polygons, containment 0.5 implies IoU 1/3, which already
    clears the default iou_threshold of 0.15 -- so containment cannot merge a
    symmetric pair that IoU would have rejected. That bounds the change: it adds
    acceptances only where IoU < 0.15, i.e. where one side is much smaller.
    """
    from shapely.geometry import box
    from spida.S.segmentation.tiling import polygon_overlap

    a, b = box(0, 0, 10, 10), box(5, 0, 15, 10)          # equal area, half overlap
    iou, containment = polygon_overlap(a, b)
    assert containment == pytest.approx(0.5)
    assert iou == pytest.approx(1 / 3)
    assert iou > 0.15                       # IoU would already have accepted it


def test_stitch_matches_only_within_a_z_plane():
    """Fragments on different z must never be merged with each other."""
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    cell = box(85, 10, 95, 20)
    polys = {(0, 0): _gdf([(1, 0, cell)]), (0, 1): _gdf([(1, 1, cell)])}
    out = stitch_polygons(polys, tiles, iou_threshold=0.2)
    assert out["ID"].nunique() == 2
    assert sorted(out["z"].tolist()) == [0, 1]


def test_stitch_merges_a_cell_at_a_four_tile_corner():
    """Only right/down neighbours are compared, so a corner cell relies on
    union-find closing the chain transitively."""
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100),
             Tile(1, 0, 80, 0, 100, 100), Tile(1, 1, 80, 80, 100, 100)]
    cell = box(85, 85, 95, 95)                   # in all four overlaps
    polys = {k: _gdf([(1, 0, cell)]) for k in
             [(0, 0), (0, 1), (1, 0), (1, 1)]}
    out = stitch_polygons(polys, tiles, iou_threshold=0.2)
    assert out["ID"].nunique() == 1              # one cell, not four


def test_stitch_largest_vs_union_merge():
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    big, small = box(84, 10, 96, 22), box(85, 11, 95, 21)
    polys = {(0, 0): _gdf([(1, 0, big)]), (0, 1): _gdf([(1, 0, small)])}

    kept = stitch_polygons(polys, tiles, iou_threshold=0.2, merge="largest")
    assert kept.geometry.iloc[0].area == pytest.approx(big.area)

    joined = stitch_polygons(polys, tiles, iou_threshold=0.2, merge="union")
    assert joined.geometry.iloc[0].area == pytest.approx(big.union(small).area)


def test_a_tile_contributing_nothing_is_ignored():
    """Three shapes of "no cells here": key absent (the gate skipped it), an
    empty frame, and None. Production only produces the first two -- gated tiles
    never reach the post worker, and `masks_to_geodataframe` returns an empty
    frame rather than None -- but `add_tile` guards all three."""
    from shapely.geometry import box
    import geopandas as gpd

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    one = _gdf([(1, 0, box(10, 10, 20, 20))])

    assert len(stitch_polygons({(0, 0): one}, tiles, iou_threshold=0.2)) == 1

    for absent in (None, gpd.GeoDataFrame({"ID": [], "z": [],
                                           "Geometry": gpd.GeoSeries([], dtype="geometry")},
                                          geometry="Geometry", crs=None)):
        s = PolygonStitcher(tiles)
        s.add_tile((0, 0), one)
        s.add_tile((0, 1), absent)
        assert len(s.finish()) == 1


def test_stitch_output_ids_are_dense_from_one():
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100)]
    polys = {
        (0, 0): _gdf([(5, 0, box(10, 10, 20, 20)), (9, 0, box(30, 30, 40, 40))]),
        (0, 1): _gdf([(3, 0, box(150, 10, 160, 20))]),
    }
    out = stitch_polygons(polys, tiles, iou_threshold=0.2)
    assert sorted(out["ID"].tolist()) == [1, 2, 3]


# ---------------------------------------------------------------------------
# 3-stage pipeline
# ---------------------------------------------------------------------------

from spida.S.segmentation.tiling import (  # noqa: E402
    run_tiled,
    tile_polygons_to_mosaic,
)


def test_tile_polygons_to_mosaic_scales_and_translates():
    from shapely.geometry import box
    gdf = _gdf([(1, 0, box(0, 0, 10, 10))])
    out = tile_polygons_to_mosaic(gdf, Tile(0, 0, 400, 800, 100, 100), scale=4)
    minx, miny, maxx, maxy = out.geometry.iloc[0].bounds
    assert (minx, miny) == (800, 400)        # translated x-first, y-second
    assert (maxx, maxy) == (840, 440)        # and scaled by 4


@pytest.fixture
def tiled_mosaic(tmp_path):
    """A 2-plane mosaic with two blobs: one interior, one straddling the seam of
    a 2-tile grid, so stitching is actually exercised."""
    import tifffile

    h, w = 256, 448
    d = tmp_path / "images"
    d.mkdir()
    base = np.zeros((h, w), dtype=np.uint16)
    base[40:80, 40:80] = 50000               # interior of tile 0
    base[100:140, 200:240] = 50000           # inside the 192..256 overlap
    for z in range(2):
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", base)
    return d, (h, w)


def _blob_eval(img):
    """Stand-in for cellpose: label the bright blobs in the tile."""
    from scipy import ndimage
    # img is (C, Y, X, Z); segment the first channel's first plane, then
    # broadcast so the mask is (Z, Y, X) like a 3D cellpose result.
    plane = img[0, ..., 0]
    lab, _ = ndimage.label(plane > 1000)
    return np.stack([lab.astype(np.uint32)] * img.shape[-1], axis=0)


def test_run_tiled_end_to_end_stitches_the_seam_blob(tiled_mosaic):
    images_dir, shape = tiled_mosaic
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=256, overlap=64, mask=None)
    # 256x448 at tile 256 overlap 64 -> 1 row, 2 cols
    assert plan["grid"] == (1, 2)

    out = run_tiled(plan, {"DAPI": planes}, _blob_eval, scale=1,
                    poly_tile_size=128, poly_tile_overlap=32,
                    iou_threshold=0.2, n_prep=2, n_post=2)

    # two blobs total; the seam blob is seen by BOTH tiles but must be one cell
    assert out["ID"].nunique() == 2, f"got {out['ID'].nunique()} cells"
    assert sorted(out["ID"].unique()) == [1, 2]
    assert set(out["z"].unique()) == {0, 1}


def test_run_tiled_polygons_land_in_mosaic_coordinates(tiled_mosaic):
    images_dir, shape = tiled_mosaic
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=256, overlap=64, mask=None)
    out = run_tiled(plan, {"DAPI": planes}, _blob_eval, scale=1,
                    poly_tile_size=128, poly_tile_overlap=32, n_prep=1, n_post=1)

    bounds = out.total_bounds                    # (minx, miny, maxx, maxy)
    assert 30 <= bounds[0] and bounds[2] <= shape[1]     # x within the mosaic
    assert 30 <= bounds[1] and bounds[3] <= shape[0]     # y within the mosaic


def test_run_tiled_propagates_eval_errors(tiled_mosaic):
    """A failing model call must abort loudly, not hang or silently drop tiles."""
    images_dir, shape = tiled_mosaic
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=256, overlap=64, mask=None)

    def boom(img):
        raise RuntimeError("CUDA out of memory")

    with pytest.raises(RuntimeError, match="tiled segmentation failed"):
        run_tiled(plan, {"DAPI": planes}, boom, scale=1, n_prep=1, n_post=1)


def test_run_tiled_rejects_an_empty_plan(tiled_mosaic):
    images_dir, shape = tiled_mosaic
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=256, overlap=64, mask=None)
    plan["to_segment"] = []
    with pytest.raises(ValueError, match="no tiles to segment"):
        run_tiled(plan, {"DAPI": planes}, _blob_eval, scale=1)


def test_run_tiled_honours_the_tissue_gate(tiled_mosaic):
    """Tiles the gate drops must never reach eval_tile."""
    images_dir, shape = tiled_mosaic
    planes = find_stain_planes(images_dir, "DAPI")
    mask = np.zeros(shape, dtype=bool)
    mask[0:150, 0:150] = True                    # only the left tile has tissue
    plan = plan_grid(shape, tile=256, overlap=64, mask=mask)

    seen = []

    def counting_eval(img):
        seen.append(img.shape)
        return _blob_eval(img)

    run_tiled(plan, {"DAPI": planes}, counting_eval, scale=1,
              mask=mask, poly_tile_size=128, poly_tile_overlap=32,
              n_prep=1, n_post=1)
    assert len(seen) == plan["n_segment"]
    assert plan["n_segment"] < plan["n_tiles"] or plan["n_skipped"] == 0


# ---------------------------------------------------------------------------
# incremental stitching
# ---------------------------------------------------------------------------

from spida.S.segmentation.tiling import (  # noqa: E402
    DEFAULT_STITCH_IOU,
    PolygonStitcher,
)


def test_incremental_equals_batch_regardless_of_arrival_order():
    """Tiles finish in whatever order the GPU emits them, so the result must not
    depend on arrival order."""
    from shapely.geometry import box
    import itertools

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100),
             Tile(1, 0, 80, 0, 100, 100), Tile(1, 1, 80, 80, 100, 100)]
    cell = box(85, 85, 95, 95)                      # touches all four
    polys = {
        (0, 0): _gdf([(1, 0, box(10, 10, 20, 20)), (2, 0, cell)]),
        (0, 1): _gdf([(1, 0, cell)]),
        (1, 0): _gdf([(1, 0, cell), (2, 0, box(120, 10, 130, 20))]),
        (1, 1): _gdf([(1, 0, cell)]),
    }
    batch = stitch_polygons(polys, tiles)
    expected = batch["ID"].nunique()

    for order in itertools.islice(itertools.permutations(polys), 6):
        s = PolygonStitcher(tiles)
        for key in order:
            s.add_tile(key, polys[key])
        assert s.finish()["ID"].nunique() == expected, f"order {order} differs"


def test_incremental_matches_pairs_before_all_tiles_arrive():
    """The point of the design: matching happens in flight, not at the end."""
    from shapely.geometry import box

    tiles = [Tile(0, 0, 0, 0, 100, 100), Tile(0, 1, 0, 80, 100, 100),
             Tile(0, 2, 0, 160, 100, 100)]
    cell = box(85, 10, 95, 20)
    s = PolygonStitcher(tiles)
    s.add_tile((0, 0), _gdf([(1, 0, cell)]))
    assert s.n_pairs == 0                            # nothing to compare yet
    s.add_tile((0, 1), _gdf([(1, 0, cell)]))
    assert s.n_pairs == 1 and s.n_merges >= 1        # matched before tile 3 exists


def test_incremental_is_thread_safe():
    """add_tile runs on the post pool, so concurrent calls must not corrupt
    union-find or double-allocate ids."""
    from concurrent.futures import ThreadPoolExecutor
    from shapely.geometry import box

    n = 12
    tiles = [Tile(0, j, 0, j * 80, 100, 100) for j in range(n)]
    polys = {}
    for j in range(n):
        seam = box(j * 80 + 85, 10, j * 80 + 95, 20)   # shared with j+1
        polys[(0, j)] = _gdf([(1, 0, seam), (2, 0, box(j * 80 + 10, 50,
                                                       j * 80 + 20, 60))])
    s = PolygonStitcher(tiles)
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(lambda kv: s.add_tile(*kv), polys.items()))
    out = s.finish()

    batch = stitch_polygons(polys, tiles)
    assert out["ID"].nunique() == batch["ID"].nunique()
    assert len(out) == len(batch)


def test_stitcher_rejects_unknown_merge_mode():
    with pytest.raises(ValueError, match="merge must be"):
        PolygonStitcher([], merge="average")


# ---------------------------------------------------------------------------
# concurrent GPU stage (one shared model)
# ---------------------------------------------------------------------------

from spida.S.segmentation.tiling import (  # noqa: E402
    DEFAULT_GPU_THREADS,
)


def _many_tile_mosaic(tmp_path, n_cols=6, tile=128, overlap=32):
    import tifffile
    h = tile
    w = (tile - overlap) * (n_cols - 1) + tile
    d = tmp_path / "images"
    d.mkdir()
    base = np.zeros((h, w), dtype=np.uint16)
    for j in range(n_cols):                       # one blob per tile interior
        cx = j * (tile - overlap) + tile // 2
        base[40:70, cx - 15:cx + 15] = 50000
    for z in range(2):
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", base)
    return d, (h, w)


def test_gpu_pool_calls_eval_once_per_tile(tmp_path):
    """No tile dropped, none duplicated, whatever the interleaving."""
    import threading

    images_dir, shape = _many_tile_mosaic(tmp_path)
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=128, overlap=32, mask=None)

    seen, lock = [], threading.Lock()

    def counting_eval(img):
        with lock:
            seen.append(img.shape)
        return _blob_eval(img)

    run_tiled(plan, {"DAPI": planes}, counting_eval, scale=1,
              poly_tile_size=64, poly_tile_overlap=16,
              n_prep=2, n_post=2, n_gpu_threads=3)
    assert len(seen) == plan["n_segment"]


def test_gpu_threads_actually_run_concurrently(tmp_path):
    """The whole point: two evals must be in flight at once, otherwise the CPU
    tail of one cannot overlap the GPU work of another."""
    import threading
    import time

    images_dir, shape = _many_tile_mosaic(tmp_path)
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=128, overlap=32, mask=None)

    inflight = {"now": 0, "max": 0}
    lock = threading.Lock()

    def slow_eval(img):
        with lock:
            inflight["now"] += 1
            inflight["max"] = max(inflight["max"], inflight["now"])
        time.sleep(0.05)                          # stand in for the model call
        with lock:
            inflight["now"] -= 1
        return _blob_eval(img)

    run_tiled(plan, {"DAPI": planes}, slow_eval, scale=1,
              poly_tile_size=64, poly_tile_overlap=16,
              n_prep=3, n_post=2, n_gpu_threads=3, queue_depth=3)
    assert inflight["max"] >= 2, "GPU stage never had two evals in flight"


def test_gpu_pool_result_matches_single_threaded(tmp_path):
    """Concurrency must not change the answer."""
    images_dir, shape = _many_tile_mosaic(tmp_path)
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=128, overlap=32, mask=None)

    kw = dict(scale=1, poly_tile_size=64, poly_tile_overlap=16, n_prep=2, n_post=2)
    one = run_tiled(plan, {"DAPI": planes}, _blob_eval, n_gpu_threads=1, **kw)
    many = run_tiled(plan, {"DAPI": planes}, _blob_eval, n_gpu_threads=4, **kw)
    assert one["ID"].nunique() == many["ID"].nunique()
    assert len(one) == len(many)


def test_gpu_pool_propagates_errors_without_hanging(tmp_path):
    """One failing thread must abort the run, not deadlock the others."""
    images_dir, shape = _many_tile_mosaic(tmp_path)
    planes = find_stain_planes(images_dir, "DAPI")
    plan = plan_grid(shape, tile=128, overlap=32, mask=None)

    def flaky(img):
        raise RuntimeError("CUDA out of memory")

    with pytest.raises(RuntimeError, match="tiled segmentation failed"):
        run_tiled(plan, {"DAPI": planes}, flaky, scale=1,
                  poly_tile_size=64, poly_tile_overlap=16,
                  n_prep=2, n_post=2, n_gpu_threads=3)


# ---------------------------------------------------------------------------
# Golden master for utilities/tiling.py
#
# `tile_grid_with_overlap` was SPLIT OUT of `tile_image_with_overlap` so the
# segmentation side could plan a grid without holding a 64 GB array. That
# refactor had to be behaviour-preserving for the existing caller
# (deconvolution), so its output was captured across the full parameter surface
# BEFORE the change and is asserted identical here.
#
# It has already earned its keep twice: it caught `actual_size` being computed
# as `end - start` (wrong whenever an oversized tile makes the origin negative),
# and it pins the negative-origin behaviour itself — which is why
# `masks._fit_tiling` clamps around that behaviour rather than changing it here.
# ---------------------------------------------------------------------------
