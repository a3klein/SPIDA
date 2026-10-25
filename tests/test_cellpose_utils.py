from __future__ import annotations

from pathlib import Path

import numpy as np
import imageio.v3 as iio
import skimage as ski
import pytest
from cellpose import core

pytest.importorskip("cellpose")

from spida.S.segmentation.backends import cellpose
from spida.S.segmentation.backends.cellpose_config import CellposeConfig


def _write_tif(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(path, array)


def _load_region(images_dir: Path, **cfg_kw) -> np.ndarray:
    """Load a region's stains the way `_segment_whole` does.

    Three production steps: `_resolve_stains` -> `read_stains` -> the
    single-plane squeeze. The squeeze is repeated here rather than imported
    because it lives inline in `_segment_whole`, which also runs a model; the
    guard against the two drifting apart is
    `test_untiled_path_squeezes_a_single_plane_stain` in
    test_cellpose_backend_wiring.py, which asserts on the shape the MODEL sees.
    """
    from spida.io.images import read_stains

    cfg = CellposeConfig.from_kwargs(**cfg_kw)
    stains = cellpose._resolve_stains(cfg, images_dir)
    img = read_stains(stains, downscale=cfg.scale)
    return img[..., 0] if img.shape[-1] == 1 else img


def test_read_stains_downscale_2d_single_channel(tmp_path: Path):
    """`read_stains` is always 4D -- Z=1 here. The collapse to (C, Y, X) is the
    caller's, tested via `_load_region` below and on the model input in
    test_cellpose_backend_wiring.py."""
    from spida.io.images import read_stains

    img = np.arange(16, dtype=np.uint16).reshape(4, 4)
    _write_tif(tmp_path / "img_DAPI.tif", img)

    cfg = CellposeConfig.from_kwargs(nuc_stain_name="DAPI", cyto_stain_name=None,
                                     image_ext=".tif", scale=2)
    loaded = read_stains(cellpose._resolve_stains(cfg, tmp_path), downscale=2)

    expected = ski.transform.downscale_local_mean(img, (2, 2))
    assert loaded.shape == (1, 2, 2, 1)
    assert np.allclose(loaded[0, ..., 0], expected)
    # and the caller-side collapse the 2D path depends on
    assert _load_region(tmp_path, nuc_stain_name="DAPI", cyto_stain_name=None,
                        image_ext=".tif", scale=2).shape == (1, 2, 2)


def test_read_stains_downscale_3d_stack(tmp_path: Path):
    stack = [
        np.full((4, 4), fill_value=1, dtype=np.uint16),
        np.full((4, 4), fill_value=3, dtype=np.uint16),
    ]
    for idx, plane in enumerate(stack):
        _write_tif(tmp_path / f"img_DAPI_z{idx}.tif", plane)

    loaded = _load_region(tmp_path, nuc_stain_name="DAPI", cyto_stain_name=None,
                          image_ext=".tif", scale=2)

    expected = np.stack(
        [ski.transform.downscale_local_mean(p, (2, 2)) for p in stack], axis=-1
    )
    assert loaded.shape == (1, 2, 2, 2)      # multi-plane: no squeeze
    assert np.allclose(loaded[0], expected)


def test_read_stains_downscale_with_cyto(tmp_path: Path):
    nuc = np.arange(16, dtype=np.uint16).reshape(4, 4)
    cyto = np.arange(16, dtype=np.uint16).reshape(4, 4) * 2
    _write_tif(tmp_path / "img_DAPI.tif", nuc)
    _write_tif(tmp_path / "img_PolyT.tif", cyto)

    loaded = _load_region(tmp_path, nuc_stain_name="DAPI", cyto_stain_name="PolyT",
                          image_ext=".tif", scale=2)

    expected_nuc = ski.transform.downscale_local_mean(nuc, (2, 2))
    expected_cyto = ski.transform.downscale_local_mean(cyto, (2, 2))
    assert loaded.shape == (2, 2, 2)
    assert np.allclose(loaded[0], expected_nuc)
    assert np.allclose(loaded[1], expected_cyto)


# mask->polygon tests moved to tests/test_masks.py (backends/masks.py has no
# cellpose dependency, so they run in the standard preprocessing env).


class _DummyModel:
    def __init__(self):
        self.img = None
        self.kwargs = None

    def eval(self, img, **kwargs):
        self.img = img
        self.kwargs = kwargs
        return "masks", "flows", "styles"


class _DummyResolution:
    """Stand-in for the ModelResolution that _load_model now returns."""

    name = "cpsam"
    pretrained_model = "/models/cpsam"
    model_dir = Path("/models")
    exists_locally = True
    will_download = False
    size_bytes = 1

    def log(self):
        pass


def _patch_load_model(monkeypatch, dummy):
    monkeypatch.setattr(cellpose, "_load_model",
                        lambda cfg: (dummy, _DummyResolution()))


def _cfg(**kwargs):
    """Build a CellposeConfig; seg_mode replaces the old do_3D/stitch flags."""
    return CellposeConfig.from_kwargs(**kwargs)


def test_cellpose_wrapper_3d_reorders_channels(monkeypatch):
    dummy = _DummyModel()
    _patch_load_model(monkeypatch, dummy)

    img = np.zeros((2, 3, 5, 7), dtype=np.float32)
    cellpose._cellpose_wrapper(img, _cfg(seg_mode="3d_true"))

    assert dummy.img.shape == (7, 3, 5, 2)
    assert dummy.kwargs["channel_axis"] == -1
    assert dummy.kwargs["z_axis"] == 0


def test_cellpose_wrapper_2d_stack_reorders_channels(monkeypatch):
    dummy = _DummyModel()
    _patch_load_model(monkeypatch, dummy)

    img = np.zeros((2, 3, 5, 7), dtype=np.float32)
    cellpose._cellpose_wrapper(img, _cfg(seg_mode="3d_stitched"))

    assert dummy.img.shape == (7, 3, 5, 2)
    assert dummy.kwargs["channel_axis"] == -1
    assert dummy.kwargs["z_axis"] == 0


def _gaussian_2d(shape, center, sigma, amplitude=1.0):
    yy, xx = np.indices(shape)
    cy, cx = center
    return amplitude * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma**2)))


def _synthetic_stack(z_slices=7, shape=(256, 256), diameters=(30, 60, 90)):
    rng = np.random.default_rng(42)
    base_centers = [
        (shape[0] * 0.3, shape[1] * 0.3),
        (shape[0] * 0.6, shape[1] * 0.5),
        (shape[0] * 0.4, shape[1] * 0.75),
    ]
    dapi = np.zeros((z_slices, *shape), dtype=np.float32)
    polyt = np.zeros((z_slices, *shape), dtype=np.float32)
    centers_per_z: list[list[tuple[float, float]]] = []

    for z in range(z_slices):
        z_factor = 1.0 - abs(z - (z_slices - 1) / 2) / z_slices
        centers_this_z = []
        for (cy, cx), diameter in zip(base_centers, diameters):
            sigma = diameter / 2.355
            jitter = rng.normal(scale=1.5, size=2)
            center = (cy + jitter[0], cx + jitter[1])
            centers_this_z.append(center)
            dapi[z] += _gaussian_2d(shape, center, sigma, amplitude=1.0 * z_factor)
            polyt[z] += _gaussian_2d(shape, center, sigma * 1.3, amplitude=0.6 * z_factor)
        centers_per_z.append(centers_this_z)

    dapi = (dapi / dapi.max()).astype(np.float32)
    polyt = (polyt / polyt.max()).astype(np.float32)
    return np.stack([dapi, polyt], axis=0), centers_per_z


def _mask_contains_center(masks: np.ndarray, center: tuple[float, float], radius: int = 3):
    y, x = int(round(center[0])), int(round(center[1]))
    y0 = max(0, y - radius)
    y1 = min(masks.shape[-2], y + radius + 1)
    x0 = max(0, x - radius)
    x1 = min(masks.shape[-1], x + radius + 1)
    if y0 >= y1 or x0 >= x1:
        return False
    window = masks[y0:y1, x0:x1]
    return np.any(window > 0)


@pytest.mark.slow
def test_cellpose_3d_vs_stitched_2d_on_synthetic_stack():
    if not core.use_gpu():
        pytest.skip("GPU not available for Cellpose model")

    img, centers_per_z = _synthetic_stack(z_slices=7, shape=(256, 256))

    masks_3d, _, _ = cellpose._cellpose_wrapper(
        img,
        do_3D=True,
        flow_threshold=0.0,
        cellprob_threshold=-4,
        min_size=100,
    )

    masks_2d, _, _ = cellpose._cellpose_wrapper(
        img,
        do_3D=False,
        flow_threshold=0.0,
        cellprob_threshold=-4,
        min_size=100,
        stitch_threshold=0.25,
    )

    assert masks_3d.ndim == 3
    assert masks_2d.ndim == 3
    assert masks_3d.shape == masks_2d.shape
    assert np.unique(masks_3d).size > 1
    assert np.unique(masks_2d).size > 1

    for z, centers in enumerate(centers_per_z):
        for center in centers:
            assert _mask_contains_center(masks_3d[z], center)
            assert _mask_contains_center(masks_2d[z], center)


# --- flow3D_smooth reaches model.eval, and only on the do_3D path ---------
def test_flow3d_smooth_passed_only_for_do_3d():
    """Z-ONLY by construction: the ring artifact is a z-axis flow-noise effect,
    and blurring the high-resolution XY flow risks merging touching nuclei. The
    config exposes a scalar sigma_z so XY can never be set by accident."""
    img = np.zeros((1, 8, 8, 4), dtype=np.uint16)          # (C, Y, X, Z)

    _, kw = cellpose.build_eval_call(
        _cfg(seg_mode="3d_true", flow3d_smooth_z=2), img, anisotropy=3.5)
    assert kw["flow3D_smooth"] == [2.0, 0.0, 0.0]
    assert kw["do_3D"] is True

    # absent when disabled
    _, kw = cellpose.build_eval_call(
        _cfg(seg_mode="3d_true", flow3d_smooth_z=0), img, anisotropy=3.5)
    assert "flow3D_smooth" not in kw

    # absent on the stitched path even if set
    _, kw = cellpose.build_eval_call(_cfg(seg_mode="3d_stitched",
                                          flow3d_smooth_z=2), img)
    assert "flow3D_smooth" not in kw


def test_flow3d_smooth_is_a_real_cellpose_eval_kwarg():
    """Guards against the kwarg being silently swallowed by **kwargs if a future
    cellpose renames or drops it -- which would revert the ring fix invisibly,
    since cellpose logs nothing when it applies flow3D_smooth."""
    import inspect
    from cellpose import models
    assert "flow3D_smooth" in inspect.signature(models.CellposeModel.eval).parameters
