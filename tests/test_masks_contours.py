"""`_contours_to_geometry` — every shape a single-z cross-section can take.

These exist because the previous implementation was `Polygon(contours[0])`: the
first contour only, no interior rings. `find_contours` orders by RASTER SCAN, not
by size, so "first" means topmost — a label whose cross-section was two blobs
kept only the upper one, and a ring became a filled disc.

Every case asserts the polygon area against the TRUE PIXEL COUNT, which is the
invariant that was being violated. Tolerance is ~0.5 px^2 per boundary: contours
trace the 0.5 iso-level, i.e. through pixel centres, so each part or hole costs
half a pixel. That offset predates this fix and is not what these guard.
"""
from __future__ import annotations

import numpy as np
import pytest
import skimage as ski

from spida.S.segmentation.backends.masks import _contours_to_geometry


def _contours(mask: np.ndarray):
    """Exactly what convert_mask feeds the helper."""
    padded = np.pad(mask, 1, mode="constant")
    return [c - 1 for c in ski.measure.find_contours(padded, 0.5)]


def _area(mask: np.ndarray) -> float:
    g = _contours_to_geometry(_contours(mask))
    assert g is not None, "no geometry produced"
    assert g.is_valid, f"invalid geometry: {g.geom_type}"
    return g.area


def _ring(size, cx, cy, r_out, r_in):
    yy, xx = np.mgrid[0:size, 0:size]
    r = np.hypot(yy - cy, xx - cx)
    m = np.zeros((size, size), dtype=bool)
    m[(r < r_out) & (r > r_in)] = True
    return m


# --- one boundary: area is true count minus ~0.5 -------------------------
def test_single_solid_blob():
    m = np.zeros((40, 40), bool); m[10:30, 10:30] = True
    assert _area(m) == pytest.approx(m.sum(), abs=0.6)


def test_blob_touching_the_array_edge():
    """convert_mask pads by 1 so an edge-touching region's contour closes."""
    m = np.zeros((30, 30), bool); m[0:12, 0:12] = True
    assert _area(m) == pytest.approx(m.sum(), abs=0.6)


@pytest.mark.parametrize("mask,label", [
    (np.pad(np.ones((1, 20), bool), ((10, 9), (5, 5))), "1-px line"),
    (np.pad(np.ones((1, 1), bool), ((10, 19), (10, 19))), "single pixel"),
])
def test_degenerate_shapes_do_not_raise(mask, label):
    """A 1-px line and a single pixel are near-zero area but must not crash or
    return an invalid geometry."""
    g = _contours_to_geometry(_contours(mask))
    assert g is None or g.is_valid, label


# --- THE REGRESSIONS: cases contours[0] got wrong ------------------------
def test_disconnected_parts_are_all_kept_regardless_of_order():
    """The headline bug. Two blobs of one label; `contours[0]` kept whichever the
    raster scan reached first, so a small blob ABOVE a large one discarded 95% of
    the cell. Both orderings must give the full area."""
    small_above = np.zeros((60, 60), bool)
    small_above[4:12, 4:12] = True; small_above[25:55, 10:50] = True

    large_above = np.zeros((60, 60), bool)
    large_above[5:35, 10:50] = True; large_above[45:53, 5:13] = True

    for m, name in ((small_above, "small above large"),
                    (large_above, "large above small")):
        assert _area(m) == pytest.approx(m.sum(), abs=1.1), name
    # and the ordering must not matter, which is what the old code got wrong
    assert _area(small_above) == pytest.approx(_area(large_above), abs=0.01)


def test_three_disconnected_parts():
    m = np.zeros((60, 60), bool)
    m[5:15, 5:15] = True; m[20:30, 20:30] = True; m[40:55, 5:25] = True
    assert _area(m) == pytest.approx(m.sum(), abs=1.6)


def test_ring_hole_is_subtracted_not_filled():
    """A ring became a filled disc: +33% area, and it manufactured phantom
    containment for any cell sitting in the hole."""
    m = _ring(60, 30, 30, 22, 11)
    assert _area(m) == pytest.approx(m.sum(), abs=1.1)


def test_blob_with_two_holes():
    m = np.zeros((60, 60), bool); m[8:52, 8:52] = True
    m[14:22, 14:22] = False; m[34:44, 34:44] = False
    assert _area(m) == pytest.approx(m.sum(), abs=1.6)


def test_part_inside_a_hole_is_added_back():
    """Nesting depth 2: the hole subtracts, the blob within it adds back. An
    orientation-plus-nesting implementation gets this wrong; symmetric
    difference gets it right for free."""
    m = _ring(60, 30, 30, 22, 11)
    m[27:33, 27:33] = True
    assert _area(m) == pytest.approx(m.sum(), abs=1.6)


def test_concentric_rings_nesting_depth_two():
    """Four contours alternating solid/hole/solid/hole."""
    m = _ring(80, 40, 40, 36, 28) | _ring(80, 40, 40, 18, 10)
    assert _area(m) == pytest.approx(m.sum(), abs=2.1)


def test_diagonally_touching_blobs_are_two_parts():
    """find_contours treats a diagonal-only touch as two regions; both count."""
    m = np.zeros((30, 30), bool); m[5:14, 5:14] = True; m[14:23, 14:23] = True
    assert _area(m) == pytest.approx(m.sum(), abs=1.1)


# --- shape of the return value ------------------------------------------
def test_multipart_returns_multipolygon_singlepart_returns_polygon():
    """MultiPolygon is NOT a new type downstream -- `dissolve` and the stitcher's
    union already produce it, and `ingest_polygons` coerces everything to
    MultiPolygon anyway. This just pins which case yields which."""
    solid = np.zeros((40, 40), bool); solid[10:30, 10:30] = True
    split = np.zeros((60, 60), bool); split[4:12, 4:12] = True; split[25:55, 10:50] = True
    assert _contours_to_geometry(_contours(solid)).geom_type == "Polygon"
    assert _contours_to_geometry(_contours(split)).geom_type == "MultiPolygon"


def test_no_usable_contour_returns_none():
    """Empty input, and a mask whose only contour has too few points, both give
    None rather than raising -- convert_mask skips the plane."""
    assert _contours_to_geometry([]) is None
    assert _contours_to_geometry([np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])]) is None


def test_degenerate_first_contour_no_longer_costs_the_plane():
    """The old guard tested `len(contours[0]) > 3` and skipped the ENTIRE plane
    when it failed. Filtering per contour means a real region behind a
    too-small first contour still produces geometry."""
    big = np.zeros((40, 40), bool); big[10:30, 10:30] = True
    tiny = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])       # 3 points
    g = _contours_to_geometry([tiny] + _contours(big))
    assert g is not None and g.area == pytest.approx(big.sum(), abs=0.6)
