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


# ---------------------------------------------------------------------------
# GPU-tiled path: run_cellpose dispatch
#
# The whole point of the tiled path is that a cell straddling a tile edge is
# segmented twice, by two INDEPENDENT model.eval calls, and must come back as
# one cell. So the fake model here is not a canned array: it labels whatever it
# is handed, which lets a single painted square genuinely appear in four tiles.
# ---------------------------------------------------------------------------

class _LabellingModel:
    """Connected-component labels whatever it is given. Counts its calls."""

    def __init__(self):
        self.n_eval = 0
        self.eval_shapes = []
        self.last_img = None          # for tests that assert on the PIXELS given

    def eval(self, img, **kwargs):
        import scipy.ndimage as ndi

        self.n_eval += 1
        self.eval_shapes.append(tuple(img.shape))
        self.last_img = np.asarray(img).copy()
        arr = np.asarray(img)
        if arr.ndim >= 3:
            # (Y, X, C) or (Z, Y, X, C) — build_eval_call always puts channels
            # last, so dropping that axis leaves the spatial array. Labelling a
            # (Z, Y, X) stack in 3D mimics what stitch_threshold does: one label
            # per object across planes.
            arr = arr[..., 0]
        labels, _n = ndi.label(arr > 0)
        return labels.astype(np.uint32), None, None


@pytest.fixture
def mosaic_dir(tmp_path) -> Path:
    """A 64x64 mosaic, 3 z-planes, with three bright squares.

    Placement is deliberate for a tile=40 / overlap=8 grid (tiles at 0 and 24 on
    each axis, so four tiles):
      * centre square [28:38] is inside ALL FOUR tiles -> must stitch to 1 cell
      * corner squares sit inside exactly one tile each
    """
    import tifffile

    d = tmp_path / "images"
    d.mkdir()
    plane = np.zeros((64, 64), dtype=np.uint16)
    plane[28:38, 28:38] = 400                # seen by every tile
    plane[4:14, 4:14] = 400                  # tile (0,0) only
    plane[50:60, 50:60] = 400                # tile (1,1) only
    for z in range(3):
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", plane)
    return d


def _write_mask(images_dir: Path, mask: np.ndarray) -> Path:
    import tifffile

    path = images_dir / "tissue_mask.tif"
    tifffile.imwrite(path, (mask.astype(np.uint8) * 255), compression="deflate")
    return path


def _patch_labelling_model(monkeypatch):
    model = _LabellingModel()
    calls = {"n_load": 0}

    def fake_load(cfg):
        calls["n_load"] += 1
        return model, _FakeResolution()

    monkeypatch.setattr(cpb, "_load_model", fake_load)
    monkeypatch.setattr(masks_mod, "max_cpu", 1, raising=False)
    return model, calls


def _run_tiled(images_dir, out, **kw):
    # poly_tile_size must stay at or below gpu_tile_size / scale — that is the
    # min_gpu_tile() floor, and at these toy sizes the default 1000 would trip it.
    params = dict(seg_mode="2d", scale=1, min_size=0, min_z=1,
                  gpu_tile_size=40, gpu_tile_overlap=8, n_gpu_threads=2,
                  poly_tile_size=20, poly_tile_overlap=5)
    params.update(kw)
    return cpb.run_cellpose(str(images_dir), str(out), "region_001", **params)


def test_tiled_path_stitches_a_cell_seen_by_every_tile(monkeypatch, mosaic_dir, tmp_path):
    import geopandas as gpd

    model, _calls = _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"
    _run_tiled(mosaic_dir, out)

    assert model.n_eval == 4                       # 2x2 grid, all tissue
    gdf = gpd.read_parquet(out / "polygons.parquet")
    # 3 squares, not 6: the centre one appears in 4 tiles and must collapse.
    assert gdf["ID"].nunique() == 3


def test_tiled_path_loads_the_model_exactly_once(monkeypatch, mosaic_dir, tmp_path):
    """Every GPU thread shares one model; per-tile loading would dominate."""
    model, calls = _patch_labelling_model(monkeypatch)
    _run_tiled(mosaic_dir, tmp_path / "out")
    assert calls["n_load"] == 1
    assert model.n_eval == 4


def test_tiled_polygons_land_in_full_resolution_mosaic_coordinates(
        monkeypatch, mosaic_dir, tmp_path):
    import geopandas as gpd

    _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"
    _run_tiled(mosaic_dir, out)

    gdf = gpd.read_parquet(out / "polygons.parquet")
    # the corner square was painted at [50:60, 50:60] of the MOSAIC
    far = gdf[gdf.geometry.bounds["minx"] > 45]
    assert len(far) == 1
    bounds = far.geometry.bounds.iloc[0]
    assert bounds["minx"] == pytest.approx(50, abs=1.5)
    assert bounds["miny"] == pytest.approx(50, abs=1.5)


def test_tiled_meta_records_the_plan(monkeypatch, mosaic_dir, tmp_path):
    _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"
    _run_tiled(mosaic_dir, out)

    meta = json.loads((out / "meta.json").read_text())
    assert meta["path"] == "tiled"
    assert meta["config"]["gpu_tile_size"] == 40
    tiling = meta["tiling"]
    assert tiling["n_segment"] == 4
    assert tiling["overlap"] == 8
    assert tiling["tissue_mask_path"] is None
    assert tiling["iou_threshold"] == 0.15
    assert len(tiling["segment_tiles"]) == 4
    assert meta["masks"]["n_labels"] == 3


def test_untiled_path_is_still_the_default(monkeypatch, mosaic_dir, tmp_path):
    """gpu_tile_size defaults to None, so nothing changes for existing runs."""
    model, _calls = _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"
    cpb.run_cellpose(str(mosaic_dir), str(out), "r", seg_mode="2d", scale=1,
                     min_size=0, min_z=1, poly_tile_size=20, poly_tile_overlap=5)
    assert model.n_eval == 1
    meta = json.loads((out / "meta.json").read_text())
    assert meta["path"] == "single"
    assert "tiling" not in meta


def test_tissue_mask_restricts_the_grid_and_the_result(monkeypatch, mosaic_dir, tmp_path):
    """The mask sets the bbox, so the square outside it is never segmented."""
    import geopandas as gpd

    _patch_labelling_model(monkeypatch)
    mask = np.zeros((64, 64), dtype=bool)
    mask[:45, :45] = True                          # excludes the [50:60] square
    mask_path = _write_mask(mosaic_dir, mask)
    out = tmp_path / "out"

    _run_tiled(mosaic_dir, out, tissue_mask_path=mask_path)

    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert gdf["ID"].nunique() == 2
    meta = json.loads((out / "meta.json").read_text())
    assert meta["tiling"]["bbox"] == [0, 45, 0, 45]
    assert meta["tiling"]["tissue_mask_path"] == str(mask_path)


def test_tissue_mask_zeroes_signal_outside_it(monkeypatch, mosaic_dir, tmp_path):
    """Out-of-mask pixels must be blanked, not merely un-gridded."""
    import geopandas as gpd

    _patch_labelling_model(monkeypatch)
    mask = np.ones((64, 64), dtype=bool)
    mask[45:, 45:] = False                         # punches out the far square
    mask_path = _write_mask(mosaic_dir, mask)
    out = tmp_path / "out"

    _run_tiled(mosaic_dir, out, tissue_mask_path=mask_path)

    gdf = gpd.read_parquet(out / "polygons.parquet")
    # bbox still spans the mosaic, so that tile IS segmented — but sees nothing.
    assert gdf["ID"].nunique() == 2
    assert json.loads((out / "meta.json").read_text())["tiling"]["bbox"] == [0, 64, 0, 64]


def test_mask_shape_mismatch_raises(monkeypatch, mosaic_dir, tmp_path):
    _patch_labelling_model(monkeypatch)
    mask_path = _write_mask(mosaic_dir, np.ones((32, 32), dtype=bool))
    with pytest.raises(ConfigError, match="pixel-aligned"):
        _run_tiled(mosaic_dir, tmp_path / "out", tissue_mask_path=mask_path)


def test_auto_tile_size_without_a_mask_raises(monkeypatch, mosaic_dir, tmp_path):
    """'auto' optimises AGAINST the tissue; with no mask it must not guess."""
    _patch_labelling_model(monkeypatch)
    with pytest.raises(ConfigError, match="tissue-mask"):
        _run_tiled(mosaic_dir, tmp_path / "out", gpu_tile_size="auto")


def test_auto_tile_size_uses_the_mask(monkeypatch, mosaic_dir, tmp_path):
    _patch_labelling_model(monkeypatch)
    mask = np.zeros((64, 64), dtype=bool)
    mask[:45, :45] = True
    mask_path = _write_mask(mosaic_dir, mask)
    out = tmp_path / "out"

    _run_tiled(mosaic_dir, out, gpu_tile_size="auto", tissue_mask_path=mask_path,
               gpu_tile_overlap=8)

    tiling = json.loads((out / "meta.json").read_text())["tiling"]
    assert tiling["bbox"] == [0, 45, 0, 45]
    assert tiling["n_segment"] >= 1


def test_untiled_path_applies_the_tissue_mask(monkeypatch, mosaic_dir,
                                              tmp_path, caplog):
    """The untiled path used to WARN that tissue_mask_path was ignored. Now that
    both paths read through `S.io.mosaic`, masking is a reader argument and comes
    for free, so the mask must actually reach the pixels handed to the model."""
    model, _calls = _patch_labelling_model(monkeypatch)
    mask = np.zeros((64, 64), dtype=bool)
    mask[:32] = True                              # top half only
    mask_path = _write_mask(mosaic_dir, mask)

    with caplog.at_level("WARNING"):
        cpb.run_cellpose(str(mosaic_dir), str(tmp_path / "o"), "r", seg_mode="2d",
                         scale=1, min_size=0, min_z=1, poly_tile_size=20,
                         poly_tile_overlap=5, tissue_mask_path=mask_path)

    assert "ignored without gpu_tile_size" not in caplog.text
    arr = model.last_img                          # (Y, X, C) at scale=1
    assert arr is not None
    assert arr[32:].max() == 0, "out-of-tissue rows were not zeroed"
    assert arr[:32].max() > 0, "in-tissue rows were zeroed too"


def test_meta_records_the_tissue_mask_on_both_paths(monkeypatch, mosaic_dir,
                                                    tmp_path):
    """`meta.json` must say whether a mask was applied, WHICHEVER path ran.

    The untiled path can mask now, so provenance that only existed for the tiled
    path is not enough -- and it silently was: `_segment_whole` set the key in
    its `info` dict, but `run_cellpose` copies only NAMED keys into meta, so it
    was never written. The tiled path's copy survived only because it is nested
    inside `info["tiling"]`, which is copied wholesale.
    """
    import json

    _patch_labelling_model(monkeypatch)
    mask_path = _write_mask(mosaic_dir, np.ones((64, 64), dtype=bool))

    for label, extra in (("single", {}), ("tiled", dict(gpu_tile_size=40,
                                                        gpu_tile_overlap=8))):
        out = tmp_path / f"out_{label}"
        cpb.run_cellpose(str(mosaic_dir), str(out), "r", seg_mode="2d",
                         scale=1, min_size=0, min_z=1, poly_tile_size=20,
                         poly_tile_overlap=5, tissue_mask_path=mask_path, **extra)
        meta = json.loads((out / "meta.json").read_text())
        assert meta["path"] == label
        assert meta["tissue_mask_path"] == str(mask_path), (
            f"{label} path did not record tissue_mask_path in meta.json"
        )


def test_untiled_path_squeezes_a_single_plane_stain(monkeypatch, mosaic_dir,
                                                   tmp_path):
    """The single-plane collapse, asserted on what the MODEL actually receives.

    `read_stains` always returns (C, Y, X, Z); `_segment_whole` collapses Z=1 to
    (C, Y, X), and `build_eval_call` dispatches its 2D branch on that ndim. This
    is the --force_2d_segmentation shape, where IMAGE_EXT names one z-slice file.
    Asserting on the model's input rather than on a loader's return value means
    the test still holds whatever the load path is refactored into -- it is what
    keeps `_load_region` in test_cellpose_utils.py from drifting from production.
    """
    import shutil
    # NOT tmp_path/"images" -- the mosaic_dir fixture already owns that name.
    d = tmp_path / "one_plane_images"
    d.mkdir()
    shutil.copy(mosaic_dir / "mosaic_DAPI_z0.tif", d / "mosaic_DAPI_z0.tif")

    model, _calls = _patch_labelling_model(monkeypatch)
    cpb.run_cellpose(str(d), str(tmp_path / "o"), "r", seg_mode="2d",
                     scale=1, min_size=0, min_z=1, poly_tile_size=20,
                     poly_tile_overlap=5)

    # build_eval_call moves channels last for 2D: (C, Y, X) -> (Y, X, C).
    # A 4-D eval shape here would mean the squeeze did not happen.
    assert model.eval_shapes, "model was never called"
    assert len(model.eval_shapes[0]) == 3, (
        f"expected a 2D (Y, X, C) eval, got {model.eval_shapes[0]} -- the "
        f"single-plane squeeze in _segment_whole did not fire"
    )


def test_untiled_path_rejects_a_misaligned_tissue_mask(monkeypatch, mosaic_dir,
                                                       tmp_path):
    """A mask on a different pixel grid would silently mask the wrong region."""
    from spida.S.segmentation.backends.config import ConfigError

    _patch_labelling_model(monkeypatch)
    mask_path = _write_mask(mosaic_dir, np.ones((32, 32), dtype=bool))
    with pytest.raises(ConfigError, match="pixel-aligned"):
        cpb.run_cellpose(str(mosaic_dir), str(tmp_path / "o"), "r", seg_mode="2d",
                         scale=1, min_size=0, min_z=1, poly_tile_size=20,
                         poly_tile_overlap=5, tissue_mask_path=mask_path)


def test_tiled_path_handles_zero_cells(monkeypatch, tmp_path, caplog):
    import geopandas as gpd
    import tifffile

    d = tmp_path / "images"
    d.mkdir()
    for z in range(3):
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", np.zeros((64, 64), dtype=np.uint16))
    _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"

    with caplog.at_level("WARNING"):
        _run_tiled(d, out)

    assert "NO cells" in caplog.text
    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert len(gdf) == 0
    assert {"ID", "z"}.issubset(gdf.columns)       # readable by ingest_polygons


def test_tiled_3d_stitched_keeps_z_and_applies_min_z(monkeypatch, mosaic_dir, tmp_path):
    """The 3D mask path through tiling: z survives to the min_z filter."""
    import geopandas as gpd

    _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"
    _run_tiled(mosaic_dir, out, seg_mode="3d_stitched", min_z=3, n_z_planes=3)

    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert set(gdf["z"]) == {0, 1, 2}              # all three planes are painted
    assert gdf["ID"].nunique() == 3
    assert json.loads((out / "meta.json").read_text())["masks"]["n_labels"] == 3


def test_tiled_min_z_drops_a_shallow_cell(monkeypatch, tmp_path):
    """A cell on one plane only is dropped, exactly as on the untiled path."""
    import geopandas as gpd
    import tifffile

    d = tmp_path / "images"
    d.mkdir()
    for z in range(3):
        plane = np.zeros((64, 64), dtype=np.uint16)
        plane[4:14, 4:14] = 400                    # every plane -> kept
        if z == 0:
            plane[50:60, 50:60] = 400              # one plane   -> dropped
        tifffile.imwrite(d / f"mosaic_DAPI_z{z}.tif", plane)
    _patch_labelling_model(monkeypatch)
    out = tmp_path / "out"

    _run_tiled(d, out, seg_mode="3d_stitched", min_z=3, n_z_planes=3)

    gdf = gpd.read_parquet(out / "polygons.parquet")
    assert gdf["ID"].nunique() == 1
