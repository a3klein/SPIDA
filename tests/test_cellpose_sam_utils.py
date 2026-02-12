from __future__ import annotations

from pathlib import Path

import numpy as np
import imageio.v3 as iio
import skimage as ski
import pytest
from cellpose import core

pytest.importorskip("cellpose")

from spida.S.segmentation import cellposeSAM


def _write_tif(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(path, array)


def test_load_image_downscale_2d_single_channel(tmp_path: Path):
    img = np.arange(16, dtype=np.uint16).reshape(4, 4)
    _write_tif(tmp_path / "img_DAPI.tif", img)

    loaded = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name=None,
        downscale=2,
    )

    expected = ski.transform.downscale_local_mean(img, (2, 2))
    assert loaded.shape == (1, 2, 2)
    assert np.allclose(loaded[0], expected)


def test_load_image_downscale_3d_stack(tmp_path: Path):
    stack = [
        np.full((4, 4), fill_value=1, dtype=np.uint16),
        np.full((4, 4), fill_value=3, dtype=np.uint16),
    ]
    for idx, plane in enumerate(stack):
        _write_tif(tmp_path / f"img_DAPI_z{idx}.tif", plane)

    loaded = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name=None,
        downscale=2,
    )

    expected = np.stack(
        [ski.transform.downscale_local_mean(p, (2, 2)) for p in stack], axis=-1
    )
    assert loaded.shape == (1, 2, 2, 2)
    assert np.allclose(loaded[0], expected)


def test_load_image_downscale_with_cyto(tmp_path: Path):
    nuc = np.arange(16, dtype=np.uint16).reshape(4, 4)
    cyto = np.arange(16, dtype=np.uint16).reshape(4, 4) * 2
    _write_tif(tmp_path / "img_DAPI.tif", nuc)
    _write_tif(tmp_path / "img_PolyT.tif", cyto)

    loaded = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name="PolyT",
        downscale=2,
    )

    expected_nuc = ski.transform.downscale_local_mean(nuc, (2, 2))
    expected_cyto = ski.transform.downscale_local_mean(cyto, (2, 2))
    assert loaded.shape == (2, 2, 2)
    assert np.allclose(loaded[0], expected_nuc)
    assert np.allclose(loaded[1], expected_cyto)


def test_masks_to_geodataframe_simple(tmp_path: Path):
    masks = np.zeros((32, 32), dtype=np.uint16)
    masks[2:8, 2:8] = 1
    masks[12:20, 12:20] = 2

    cellposeSAM.max_cpu = 1
    gdf = cellposeSAM.masks_to_geodataframe(masks, tile_size=32, overlap=8)

    assert not gdf.empty
    assert set(gdf["ID"].unique()) == {1, 2}
    assert gdf.geometry.notna().all()


def test_masks_to_geodataframe_upscale_preserves_location(tmp_path: Path):
    scale = 4
    masks = np.zeros((1100, 1100), dtype=np.uint16)
    masks[200:240, 220:260] = 1
    masks[800:840, 900:940] = 2

    cellposeSAM.max_cpu = 1
    gdf = cellposeSAM.masks_to_geodataframe(masks, tile_size=64, overlap=16)
    gdf.geometry = gdf.geometry.affine_transform([scale, 0, 0, scale, 0, 0])

    assert set(gdf["ID"].unique()) == {1, 2}

    expected_centroids = {}
    for mask_id in (1, 2):
        coords = np.column_stack(np.where(masks == mask_id))
        # coords are (row, col); geometry uses (x=col, y=row)
        centroid_rc = coords.mean(axis=0)
        expected_centroids[mask_id] = np.array(
            [centroid_rc[1] * scale, centroid_rc[0] * scale]
        )

    for _, row in gdf.iterrows():
        geom = row["Geometry"]
        centroid = np.array([geom.centroid.x, geom.centroid.y])
        assert np.allclose(centroid, expected_centroids[row["ID"]], atol=scale)


class _DummyModel:
    def __init__(self):
        self.img = None
        self.kwargs = None

    def eval(self, img, **kwargs):
        self.img = img
        self.kwargs = kwargs
        return "masks", "flows", "styles"


def test_cellpose_wrapper_3d_reorders_channels(monkeypatch):
    dummy = _DummyModel()
    monkeypatch.setattr(cellposeSAM, "_load_model", lambda **_: dummy)

    img = np.zeros((2, 3, 5, 7), dtype=np.float32)
    cellposeSAM._cellpose_wrapper(img, do_3D=True)

    assert dummy.img.shape == (7, 3, 5, 2)
    assert dummy.kwargs["channel_axis"] == -1
    assert dummy.kwargs["z_axis"] == 0


def test_cellpose_wrapper_2d_stack_reorders_channels(monkeypatch):
    dummy = _DummyModel()
    monkeypatch.setattr(cellposeSAM, "_load_model", lambda **_: dummy)

    img = np.zeros((2, 3, 5, 7), dtype=np.float32)
    cellposeSAM._cellpose_wrapper(img, do_3D=False)

    assert dummy.img.shape == (7, 3, 5, 2)
    assert dummy.kwargs["channel_axis"] == -1
    assert dummy.kwargs["z_axis"] == 0


def test_pipeline_2d_sets_channel_axis_and_no_z_axis(tmp_path: Path, monkeypatch):
    dummy = _DummyModel()
    monkeypatch.setattr(cellposeSAM, "_load_model", lambda **_: dummy)

    nuc = np.arange(16, dtype=np.uint16).reshape(4, 4)
    cyto = np.arange(16, dtype=np.uint16).reshape(4, 4) * 2
    _write_tif(tmp_path / "img_DAPI.tif", nuc)
    _write_tif(tmp_path / "img_PolyT.tif", cyto)

    img = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name="PolyT",
        downscale=None,
    )

    cellposeSAM._cellpose_wrapper(img, do_3D=False)

    assert dummy.img.shape == (4, 4, 2)
    assert dummy.kwargs["channel_axis"] == -1
    assert "z_axis" not in dummy.kwargs


def test_pipeline_3d_stitched_sets_channel_and_z_axis(tmp_path: Path, monkeypatch):
    dummy = _DummyModel()
    monkeypatch.setattr(cellposeSAM, "_load_model", lambda **_: dummy)

    stack = [
        np.full((4, 4), fill_value=1, dtype=np.uint16),
        np.full((4, 4), fill_value=3, dtype=np.uint16),
    ]
    for idx, plane in enumerate(stack):
        _write_tif(tmp_path / f"img_DAPI_z{idx}.tif", plane)

    img = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name=None,
        downscale=None,
    )

    cellposeSAM._cellpose_wrapper(img, do_3D=False, stitch_threshold=0.4)

    assert dummy.img.shape == (2, 4, 4, 1)
    assert dummy.kwargs["channel_axis"] == -1
    assert dummy.kwargs["z_axis"] == 0


def test_pipeline_3d_do3d_sets_channel_and_z_axis(tmp_path: Path, monkeypatch):
    dummy = _DummyModel()
    monkeypatch.setattr(cellposeSAM, "_load_model", lambda **_: dummy)

    stack = [
        np.full((4, 4), fill_value=1, dtype=np.uint16),
        np.full((4, 4), fill_value=3, dtype=np.uint16),
    ]
    for idx, plane in enumerate(stack):
        _write_tif(tmp_path / f"img_DAPI_z{idx}.tif", plane)

    img = cellposeSAM._load_image(
        tmp_path,
        image_ext=".tif",
        nuc_stain_name="DAPI",
        cyto_stain_name=None,
        downscale=None,
    )

    cellposeSAM._cellpose_wrapper(img, do_3D=True)

    assert dummy.img.shape == (2, 4, 4, 1)
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

    masks_3d, _, _ = cellposeSAM._cellpose_wrapper(
        img,
        do_3D=True,
        flow_threshold=0.0,
        cellprob_threshold=-4,
        min_size=100,
    )

    masks_2d, _, _ = cellposeSAM._cellpose_wrapper(
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
