"""Tests for the CellposeConfig -> cellpose backend wiring.

GPU-free: ``model.eval`` and the model load are mocked, so these run on the head
node. What they cover is the part that has historically gone wrong silently —
which parameters actually reach ``model.eval``, and which mode the seg_mode
selects.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("cellpose")

from spida.S.segmentation.backends import cellpose as cpb
from spida.S.segmentation.backends import masks as masks_mod
from spida.S.segmentation.backends.cellpose_config import CellposeConfig
from spida.S.segmentation.backends.config import ConfigError


def _cfg(**kw) -> CellposeConfig:
    return CellposeConfig.from_kwargs(**kw)


# ---------------------------------------------------------------------------
# build_eval_call: mode -> model.eval kwargs
# ---------------------------------------------------------------------------

def test_true_3d_sets_do_3d_and_axes():
    img = np.zeros((2, 8, 8, 7), dtype=np.uint16)          # (C, Y, X, Z)
    out, kw = cpb.build_eval_call(_cfg(seg_mode="3d_true"), img, anisotropy=3.5)
    assert out.shape == (7, 8, 8, 2)                       # -> (Z, Y, X, C)
    assert kw["do_3D"] is True
    assert kw["z_axis"] == 0 and kw["channel_axis"] == -1
    assert kw["anisotropy"] == 3.5
    assert "stitch_threshold" not in kw                    # not applicable in 3D


def test_true_3d_omits_anisotropy_when_none():
    img = np.zeros((1, 8, 8, 7), dtype=np.uint16)
    _, kw = cpb.build_eval_call(_cfg(seg_mode="3d_true"), img, anisotropy=None)
    assert "anisotropy" not in kw


def test_true_3d_rejects_single_plane():
    img = np.zeros((1, 8, 8), dtype=np.uint16)
    with pytest.raises(ConfigError, match="needs a multi-plane z-stack"):
        cpb.build_eval_call(_cfg(seg_mode="3d_true"), img)


def test_stitched_sets_stitch_threshold():
    img = np.zeros((2, 8, 8, 7), dtype=np.uint16)
    out, kw = cpb.build_eval_call(_cfg(seg_mode="3d_stitched"), img)
    assert out.shape == (7, 8, 8, 2)
    assert kw["do_3D"] is False
    assert kw["stitch_threshold"] == 0.2
    assert kw["z_axis"] == 0


def test_stitched_with_single_plane_falls_back_to_2d(caplog):
    """--force_2d_segmentation selects one z-slice file, so the default
    seg_mode must degrade gracefully rather than raise."""
    img = np.zeros((1, 8, 8), dtype=np.uint16)
    with caplog.at_level("WARNING"):
        out, kw = cpb.build_eval_call(_cfg(seg_mode="3d_stitched"), img)
    assert out.shape == (8, 8, 1)                          # (C,Y,X) -> (Y,X,C)
    assert kw["stitch_threshold"] == 0
    assert "z_axis" not in kw
    assert "single z-plane" in caplog.text


def test_2d_mode_kwargs():
    img = np.zeros((2, 8, 8), dtype=np.uint16)
    out, kw = cpb.build_eval_call(_cfg(seg_mode="2d"), img)
    assert out.shape == (8, 8, 2)
    assert kw["do_3D"] is False and kw["stitch_threshold"] == 0


def test_2d_mode_rejects_unreduced_stack():
    img = np.zeros((1, 8, 8, 7), dtype=np.uint16)
    with pytest.raises(ConfigError, match="needs a single z-plane"):
        cpb.build_eval_call(_cfg(seg_mode="2d"), img)


def test_config_values_reach_eval_kwargs():
    """The parameters that used to be swallowed by _load_model's **kwargs."""
    img = np.zeros((1, 8, 8, 7), dtype=np.uint16)
    cfg = _cfg(seg_mode="3d_true", batch_size=32, flow_threshold=0.4,
               cellprob_threshold=-1.5, tile_norm_blocksize=2960,
               diameter=30.0, cp_min_size=250)
    _, kw = cpb.build_eval_call(cfg, img, anisotropy=3.5)
    assert kw["batch_size"] == 32                          # used to raise TypeError
    assert kw["flow_threshold"] == 0.4
    assert kw["cellprob_threshold"] == -1.5
    assert kw["normalize"] == {"tile_norm_blocksize": 2960}
    assert kw["diameter"] == 30.0
    assert kw["min_size"] == 250                           # cellpose-internal


# ---------------------------------------------------------------------------
# z reduction
# ---------------------------------------------------------------------------

def test_reduce_z_middle_picks_centre_plane():
    img = np.stack([np.full((4, 4), i, dtype=np.uint16) for i in range(7)], axis=-1)
    img = img[None]                                        # (1, 4, 4, 7)
    out = cpb.reduce_z(img, "middle")
    assert out.shape == (1, 4, 4)
    assert np.all(out[0] == 3)                             # plane 7 // 2 == 3


def test_reduce_z_mip_takes_maximum():
    img = np.stack([np.full((4, 4), i, dtype=np.uint16) for i in range(7)], axis=-1)
    out = cpb.reduce_z(img[None], "mip")
    assert out.shape == (1, 4, 4)
    assert np.all(out[0] == 6)


def test_reduce_z_passes_through_single_plane():
    img = np.zeros((1, 4, 4), dtype=np.uint16)
    assert cpb.reduce_z(img, "middle").shape == (1, 4, 4)


# ---------------------------------------------------------------------------
# end-to-end run_cellpose with a mocked model
# ---------------------------------------------------------------------------

class _FakeModel:
    """Records the eval kwargs and returns a deterministic label array."""

    def __init__(self, masks):
        self._masks = masks
        self.eval_kwargs = None
        self.eval_shape = None

    def eval(self, img, **kwargs):
        self.eval_kwargs = kwargs
        self.eval_shape = img.shape
        return self._masks, None, None


class _FakeResolution:
    name = "cpsam"
    pretrained_model = "/models/cpsam"
    model_dir = Path("/models")
    exists_locally = True
    will_download = False
    size_bytes = 123

    def log(self):
        pass


@pytest.fixture
def images_dir(tmp_path) -> Path:
    """A minimal region images dir: one DAPI plane + a micron->mosaic transform."""
    import imageio.v3 as iio

    d = tmp_path / "images"
    d.mkdir()
    iio.imwrite(d / "mosaic_DAPI_z3.tif", np.zeros((16, 16), dtype=np.uint16))
    (d / "micron_to_mosaic_pixel_transform.csv").write_text(
        "9.433962 0.0 0.0\n0.0 9.433962 0.0\n0.0 0.0 1.0\n"
    )
    return d


def _two_cell_mask_2d() -> np.ndarray:
    m = np.zeros((4, 4), dtype=np.uint32)
    m[0:2, 0:2] = 1
    m[2:4, 2:4] = 2
    return m


def _patch_model(monkeypatch, masks):
    fake = _FakeModel(masks)
    monkeypatch.setattr(cpb, "_load_model", lambda cfg: (fake, _FakeResolution()))
    # keep polygon extraction single-process so the mask fixtures stay cheap
    monkeypatch.setattr(masks_mod, "max_cpu", 1, raising=False)
    return fake


def test_run_cellpose_2d_writes_outputs_and_meta(monkeypatch, images_dir, tmp_path):
    fake = _patch_model(monkeypatch, _two_cell_mask_2d())
    out = tmp_path / "out"

    is_3d = cpb.run_cellpose(str(images_dir), str(out), "region_001",
                             seg_mode="2d", scale=1, min_size=0, min_z=1)

    assert is_3d is False
    assert (out / "polygons.parquet").is_file()
    meta = json.loads((out / "meta.json").read_text())
    assert meta["region"] == "region_001"
    assert meta["config"]["seg_mode"] == "2d"
    assert meta["model"]["name"] == "cpsam"
    assert meta["masks"]["n_labels"] == 2
    assert meta["polygon_counts"]["final"] >= 1
    assert fake.eval_kwargs["do_3D"] is False


def test_run_cellpose_records_bfloat16_default(monkeypatch, images_dir, tmp_path):
    _patch_model(monkeypatch, _two_cell_mask_2d())
    out = tmp_path / "out"
    cpb.run_cellpose(str(images_dir), str(out), "r", seg_mode="2d", scale=1,
                     min_size=0, min_z=1)
    meta = json.loads((out / "meta.json").read_text())
    assert meta["config"]["use_bfloat16"] is True          # deliberate default


def test_run_cellpose_batch_size_reaches_eval(monkeypatch, images_dir, tmp_path):
    """Regression: run_cellpose used to hardcode batch_size=16 positionally, so
    passing --batch_size raised TypeError."""
    fake = _patch_model(monkeypatch, _two_cell_mask_2d())
    cpb.run_cellpose(str(images_dir), str(tmp_path / "o"), "r", seg_mode="2d",
                     scale=1, min_size=0, min_z=1, batch_size=32)
    assert fake.eval_kwargs["batch_size"] == 32


def test_run_cellpose_rejects_unknown_parameter(monkeypatch, images_dir, tmp_path):
    _patch_model(monkeypatch, _two_cell_mask_2d())
    with pytest.raises(ConfigError, match="unknown parameter"):
        cpb.run_cellpose(str(images_dir), str(tmp_path / "o"), "r",
                         seg_mode="2d", anisotrpy=3.5)


def test_run_cellpose_writes_nothing_outside_output_dir(monkeypatch, images_dir, tmp_path):
    """The caller resolves the layout; the backend must not append {region}."""
    _patch_model(monkeypatch, _two_cell_mask_2d())
    out = tmp_path / "out"
    cpb.run_cellpose(str(images_dir), str(out), "region_001", seg_mode="2d",
                     scale=1, min_size=0, min_z=1)
    assert not (out / "region_001").exists()
    assert sorted(p.name for p in out.iterdir()) == ["meta.json", "polygons.parquet"]


def test_run_cellpose_legacy_project_3d_to_2d(monkeypatch, images_dir, tmp_path):
    fake = _patch_model(monkeypatch, _two_cell_mask_2d())
    with pytest.warns(DeprecationWarning):
        cpb.run_cellpose(str(images_dir), str(tmp_path / "o"), "r",
                         project_3d_to_2d=True, scale=1, min_size=0, min_z=1)
    meta = json.loads((tmp_path / "o" / "meta.json").read_text())
    assert meta["config"]["seg_mode"] == "2d"
    assert meta["config"]["z_reduce"] == "mip"             # legacy did a projection
    assert fake.eval_kwargs["do_3D"] is False


def test_run_cellpose_3d_applies_min_z_filter(monkeypatch, images_dir, tmp_path):
    """A cell present on only one plane is dropped when min_z=2."""
    masks = np.zeros((3, 4, 4), dtype=np.uint32)
    masks[:, 0:2, 0:2] = 1          # spans all 3 planes -> kept
    masks[0, 2:4, 2:4] = 2          # one plane only     -> dropped
    _patch_model(monkeypatch, masks)
    out = tmp_path / "out"

    is_3d = cpb.run_cellpose(str(images_dir), str(out), "r", seg_mode="2d",
                             scale=1, min_size=0, min_z=2)

    assert is_3d is True            # keyed off the MASK, not the input image
    meta = json.loads((out / "meta.json").read_text())
    import geopandas as gpd
    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert set(gdf["ID"]) == {1}
    assert meta["polygon_counts"]["final"] < meta["polygon_counts"]["post_area"]


def test_run_cellpose_handles_zero_cells(monkeypatch, images_dir, tmp_path, caplog):
    """A blank mask makes masks_to_geodataframe return a frame with only tiling
    internals (no ID/z). Must not crash, and must still write a readable
    parquet that ingest_polygons can consume."""
    _patch_model(monkeypatch, np.zeros((3, 4, 4), dtype=np.uint32))
    out = tmp_path / "out"

    with caplog.at_level("WARNING"):
        cpb.run_cellpose(str(images_dir), str(out), "region_001", seg_mode="2d",
                         scale=1, min_size=0, min_z=2)

    assert "NO cells" in caplog.text
    import geopandas as gpd
    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert len(gdf) == 0
    assert {"ID", "z"}.issubset(gdf.columns)          # readable downstream
    meta = json.loads((out / "meta.json").read_text())
    assert meta["masks"]["n_labels"] == 0
    assert meta["polygon_counts"]["final"] == 0


def test_run_cellpose_anisotropy_auto_from_transform(monkeypatch, tmp_path):
    """seg_mode=3d_true derives anisotropy from the region's transform CSV."""
    import imageio.v3 as iio

    d = tmp_path / "images"
    d.mkdir()
    for z in range(3):
        iio.imwrite(d / f"mosaic_DAPI_z{z}.tif", np.zeros((8, 8), dtype=np.uint16))
    (d / "micron_to_mosaic_pixel_transform.csv").write_text(
        "9.433962 0.0 0.0\n0.0 9.433962 0.0\n0.0 0.0 1.0\n"
    )
    fake = _patch_model(monkeypatch, np.zeros((3, 8, 8), dtype=np.uint32))

    cpb.run_cellpose(str(d), str(tmp_path / "o"), "r", seg_mode="3d_true",
                     scale=4, min_size=0, min_z=1)

    # z=1.5um / (0.106um * 4) ~= 3.54
    assert fake.eval_kwargs["anisotropy"] == pytest.approx(3.538, abs=0.01)
    meta = json.loads((tmp_path / "o" / "meta.json").read_text())
    assert meta["anisotropy_resolved"] == pytest.approx(3.538, abs=0.01)
