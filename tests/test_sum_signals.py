"""Tests for the pure-Python VPT `sum-signals` reimplementation.

Two kinds of test:
  * synthetic (always runs): a uniform image with a known polygon -> known sums.
  * real-image (gated): runs on the committed 200-cell boundaries against the raw
    mosaic images if present, comparing to VPT output with a tolerance (not
    bit-exact: GDAL all_touched rasterization is version-dependent -- see the
    sum_signals docstring).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.segmentation_utils import sum_signals, _high_pass

DATA = Path(__file__).parent / "data" / "segmentation_utils"
SS_BOUNDARIES = DATA / "sum_signals_boundaries.parquet"
SS_EXPECTED = DATA / "expected_sum_signals.parquet"
REAL_IMAGES = Path(
    "/anvil/scratch/x-aklein2/BICAN/202506151211_BICAN-4x1-PU-01_VMSC31910/"
    "out/region_UCI4723/images"
)


# --------------------------------------------------------------------------- #
# synthetic (always runs)                                                     #
# --------------------------------------------------------------------------- #
def _write_synthetic(tmp_path, value=100, size=200):
    rasterio = pytest.importorskip("rasterio")
    img = np.full((size, size), value, dtype=np.uint16)
    img_path = tmp_path / "mosaic_FOO_z0.tif"
    with rasterio.open(
        img_path, "w", driver="GTiff", height=size, width=size, count=1,
        dtype="uint16",
    ) as dst:
        dst.write(img, 1)
    # identity micron->mosaic transform
    (tmp_path / "micron_to_mosaic_pixel_transform.csv").write_text(
        "1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n"
    )
    # one square cell, 50x50, at z0
    poly = geom.MultiPolygon([geom.box(20, 20, 70, 70)])
    gdf = gpd.GeoDataFrame(
        {"EntityID": [1], "ZIndex": [0], "Geometry": [poly]}, geometry="Geometry",
    )
    bpath = tmp_path / "boundaries.parquet"
    gdf.to_parquet(bpath)
    return bpath, tmp_path


def test_synthetic_uniform_image(tmp_path):
    bpath, imdir = _write_synthetic(tmp_path, value=100)
    df = sum_signals(
        bpath, imdir, imdir / "micron_to_mosaic_pixel_transform.csv", n_jobs=1
    )
    assert list(df.columns) == ["FOO_raw", "FOO_high_pass"]
    assert df.index.name == "EntityID"
    # uniform image => raw = value * (#pixels in mask); the 50x50 box masks ~2500
    # pixels (exact boundary count is GDAL-version dependent, hence the tolerance)
    pixel_count = df.loc[1, "FOO_raw"] / 100
    assert pixel_count == pytest.approx(2500, rel=0.05)
    assert pixel_count == round(pixel_count)  # integer count * uniform value
    # high-pass of a uniform region is ~0
    assert abs(df.loc[1, "FOO_high_pass"]) < 1.0


def test_high_pass_uniform_is_zero():
    window = np.full((40, 40), 500.0)
    hp = _high_pass(window)
    assert np.allclose(hp, 0.0, atol=1e-6)


def test_high_pass_nonnegative():
    rng = np.arange(64 * 64).reshape(64, 64).astype(float)
    hp = _high_pass(rng)
    assert (hp >= 0).all()


# --------------------------------------------------------------------------- #
# real-image equivalence vs VPT (gated on the raw images being present)       #
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(
    not (SS_BOUNDARIES.exists() and REAL_IMAGES.exists()),
    reason="sum-signals boundaries fixture or raw mosaic images not present",
)
def test_sum_signals_matches_vpt_within_tolerance(tmp_path):
    expected = pd.read_parquet(SS_EXPECTED)
    got = sum_signals(
        SS_BOUNDARIES, REAL_IMAGES,
        REAL_IMAGES / "micron_to_mosaic_pixel_transform.csv",
        n_jobs=4,
    )
    assert got.shape == expected.shape
    exp = expected.reindex(index=got.index, columns=got.columns)
    a = got.to_numpy(float).ravel()
    b = exp.to_numpy(float).ravel()
    # near-perfect correlation; small boundary-pixel differences accepted
    assert np.corrcoef(a, b)[0, 1] > 0.999
    rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-9)
    assert np.median(rel) < 0.02
