"""Light in-package regression for the dark-sectioning filter.

DELIBERATELY LIGHT. The authoritative numerical gate is the Octave comparison against
the upstream MATLAB (max|diff|/range = 2.98e-08, r = 1.000000000000, four optical
regimes), which lives in
`dev_notebooks/cellpose_3d_test/dark_staining/scripts/validate_against_matlab.py`.
Octave does not belong in this test suite, so these tests only check the things a
refactor here could plausibly break: shapes, dtypes, that background is actually
suppressed, that the injectable global terms are honoured, and that the pipeline
entry point writes a readable uint16 tile.
"""
import numpy as np
import pytest

from spida.S.image_prep.transforms import TileContext
from spida.S.image_prep.transforms import dark_sectioning as ds
from spida.S.image_prep.transforms.dark_sectioning import (
    DS_DEFAULTS,
    dark_sectioning_stack,
    get_atmosphere,
    get_dark_channel,
    run as dark_sectioning,
)

EM, NA, PX = 405.0, 1.4, 106.0


def _synthetic(nz=3, n=192, bg=300.0, seed=0):
    """Blobs on a haze ramp: something to sharpen and something to remove."""
    rng = np.random.default_rng(seed)
    yy, xx = np.ogrid[:n, :n]
    haze = bg + 120.0 * (xx / n)                      # left-right background ramp
    vol = np.empty((nz, n, n), np.float32)
    for z in range(nz):
        img = haze + rng.normal(0, 8, (n, n))
        for cy, cx in ((60, 60), (120, 130), (60, 140)):
            img = img + 3000.0 * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 14.0 ** 2))
        vol[z] = img
    return vol


def _kw(**over):
    kw = dict(DS_DEFAULTS)
    kw.pop("background", None)
    kw.update(over)
    return kw


def test_shape_dtype_and_finiteness_preserved():
    vol = _synthetic()
    out, meta = dark_sectioning_stack(vol, EM, NA, PX, **_kw())
    assert out.shape == vol.shape
    assert out.dtype == np.float32
    assert np.isfinite(out).all()
    assert meta["block_size"] > 0
    assert meta["dep"] == DS_DEFAULTS["dep"]


def test_single_slice_input_returns_2d():
    out, _ = dark_sectioning_stack(_synthetic(nz=1)[0], EM, NA, PX, **_kw())
    assert out.ndim == 2


def test_background_ramp_is_suppressed_relative_to_blobs():
    """The point of the method: flat haze down, structure kept.

    Scored as blob-to-background contrast rather than absolute values, since the
    output is rescaled to 0-65535 and absolute levels are therefore arbitrary.
    """
    vol = _synthetic()
    out, _ = dark_sectioning_stack(vol, EM, NA, PX, **_kw())
    yy, xx = np.ogrid[:vol.shape[1], :vol.shape[2]]
    blob = ((yy - 60) ** 2 + (xx - 60) ** 2) < 14 ** 2
    bg = ((yy - 160) ** 2 + (xx - 40) ** 2) < 14 ** 2   # empty corner

    def contrast(v):
        return float(v[:, blob].mean() / max(v[:, bg].mean(), 1e-9))

    assert contrast(out) > contrast(vol)


def test_dep_controls_correction_strength():
    """dep scales the atmosphere. The shipped 3 is above its own annotated 0.7-2
    range and caused the hollowing artifact, so this guards the knob's direction."""
    vol = _synthetic()
    lo, _ = dark_sectioning_stack(vol, EM, NA, PX, **_kw(dep=0.5))
    hi, _ = dark_sectioning_stack(vol, EM, NA, PX, **_kw(dep=3.0))
    assert not np.allclose(lo, hi)


@pytest.mark.parametrize("mode", ["abs", "pct", "otsu"])
def test_thres_modes_all_reachable(mode):
    thres = 50.0 if mode == "abs" else 30.0
    out, meta = dark_sectioning_stack(
        _synthetic(nz=1), EM, NA, PX, **_kw(thres=thres, thres_mode=mode))
    assert np.isfinite(out).all()
    assert meta["thres_mode"] == mode


def test_injected_input_range_changes_result_and_is_used():
    """Injecting whole-mosaic bounds is how per-tile seams are avoided: the input
    normalisation is an affine remap, so per-tile min/max means per-tile gain."""
    vol = _synthetic()
    base, m0 = dark_sectioning_stack(vol, EM, NA, PX, **_kw())
    wide, m1 = dark_sectioning_stack(
        vol, EM, NA, PX, **_kw(), input_range=(0.0, float(vol.max()) * 4))
    assert m0["input_min"] == pytest.approx(float(vol.min()))
    assert m1["input_min"] == 0.0
    assert not np.allclose(base, wide)


def test_injected_atmosphere_bounds_bypass_local_estimate(monkeypatch):
    calls = []
    real = ds.get_atmosphere
    monkeypatch.setattr(
        ds, "get_atmosphere",
        lambda img, dark: (calls.append(1), real(img, dark))[1])

    vol = _synthetic(nz=1)
    dark_sectioning_stack(vol, EM, NA, PX, **_kw())
    n_local = len(calls)
    assert n_local > 0, "baseline should estimate the atmosphere locally"

    calls.clear()
    dark_sectioning_stack(vol, EM, NA, PX, **_kw(), atmosphere_bounds=(1.0, 40.0))
    assert not calls, "injected bounds must skip the local estimate entirely"


def test_injected_output_max_scales_output():
    vol = _synthetic(nz=1)
    a, _ = dark_sectioning_stack(vol, EM, NA, PX, **_kw())
    b, _ = dark_sectioning_stack(vol, EM, NA, PX, **_kw(), output_max=1e9)
    assert b.max() < a.max()


def test_dark_channel_is_a_local_minimum():
    img = np.full((32, 32), 100.0)
    img[16, 16] = 5.0
    dc = get_dark_channel(img, 5)
    assert dc[16, 16] == 5.0
    assert dc[15, 15] == 5.0          # inside the window
    assert dc[0, 0] == 100.0          # outside it


def _ctx(channel="DAPI", colors=(405,)):
    return TileContext(channel=channel, colors=list(colors), z_step_um=1.5)


def test_pipeline_entry_writes_uint16_tile(tmp_path):
    """The tile loop is path-in/path-out and only reads the file back, so a writable
    uint16 with no negative wrap is the real contract."""
    import imageio.v3 as iio
    src = tmp_path / "tile_colDAPI.tif"
    iio.imwrite(src, _synthetic(nz=3).astype(np.uint16))
    dst = tmp_path / "tile_colDAPI.decon.tif"

    ret = dark_sectioning(src, dst, _ctx(), na=NA, pixel_size=PX, emwavelength=EM)
    assert str(ret) == str(dst) and dst.exists()
    out = iio.imread(dst)
    assert out.dtype == np.uint16
    assert out.shape == (3, 192, 192)
    assert out.max() > 0


def test_validate_rejects_another_transforms_params():
    """Params are no longer a shared dict every transform must tolerate."""
    with pytest.raises(ValueError, match="does not accept"):
        ds.validate({"iter": 50, "psigma": 1.5})


def test_pipeline_entry_rejects_4d(tmp_path):
    import imageio.v3 as iio
    src = tmp_path / "t_colDAPI.tif"
    iio.imwrite(src, np.zeros((2, 3, 32, 32), np.uint16))
    with pytest.raises(ValueError, match="expects"):
        dark_sectioning(src, tmp_path / "o.tif", _ctx(),
                        na=NA, pixel_size=PX, emwavelength=EM)


def test_registered_in_transform_dispatch():
    from spida.S.image_prep.transforms import TRANSFORMS, resolve_transform
    assert "dark_sectioning" in TRANSFORMS
    assert resolve_transform("dark_sectioning") is ds
    assert ds.USES_GPU is False
