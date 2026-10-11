import json
from pathlib import Path

import pytest

from spida.config import SpidaConfig
from spida.io.layout import Dataset, TreeLayout, dataset_from, element_keys, read_namespace
from spida.utilities.sd_utils import _gen_keys

EXP, REG = "202506151211_BICAN-4x1-PU-01_VMSC31910", "region_UCI5224"


def _cfg(tmp_path):
    return SpidaConfig(zarr_store=tmp_path / "z", root_path=tmp_path / "raw",
                       segmentation_store=tmp_path / "seg", image_store=tmp_path / "img")


@pytest.mark.parametrize("prefix", ["default", "cellpose_cell", "proseg_cell"])
def test_named_keys_match_gen_keys(prefix):
    assert element_keys(prefix, f"{EXP}_{REG}") == _gen_keys(prefix, EXP, REG)


def test_nameless_keys():
    k = element_keys("cellpose_cell", None)
    assert list(k.values()) == ["cellpose_cell_polygons", "cellpose_cell_transcripts", "cellpose_cell_table", "default_z3"]


def test_tree_paths(tmp_path):
    d = dataset_from(_cfg(tmp_path), EXP, REG)
    assert d.zarr == tmp_path / "z" / EXP / REG
    assert d.raw_images() == tmp_path / "raw" / EXP / "out" / REG / "images"
    assert d.seg_dir("cellpose_cell") == tmp_path / "seg" / EXP / REG / "cellpose_cell"
    assert d.figures_dir("cellpose_cell") == tmp_path / "img" / EXP / REG / "cellpose_cell"
    assert d.keys("default") == _gen_keys("default", EXP, REG)


def test_paths_override_location_names_keep_namespace(tmp_path):
    d = dataset_from(_cfg(tmp_path), EXP, REG, zarr="/elsewhere.zarr", seg="/segs", figures="/figs")
    assert d.zarr == Path("/elsewhere.zarr")
    assert d.seg_dir("cellpose_cell") == Path("/segs/cellpose_cell")
    assert d.figures_dir("x") == Path("/figs/x")
    assert d.namespace == f"{EXP}_{REG}"


def test_argument_errors(tmp_path):
    with pytest.raises(ValueError, match="both"):
        dataset_from(_cfg(tmp_path), exp=EXP)
    with pytest.raises(ValueError, match="--zarr_path"):
        dataset_from(_cfg(tmp_path))


def _store(tmp_path, attrs=None, shapes=()):
    z = tmp_path / "s.zarr"
    (z / "shapes").mkdir(parents=True)
    for s in shapes:
        (z / "shapes" / s).mkdir()
    (z / "zarr.json").write_text(json.dumps({"attributes": attrs or {}, "zarr_format": 3}))
    return z


def test_namespace_from_attrs(tmp_path):
    z = _store(tmp_path, {"spida": {"namespace": "e_r"}}, ["default_x_y_polygons"])
    assert read_namespace(z) == "e_r"
    assert dataset_from(SpidaConfig(), zarr=z).keys("default")["shapes_key"] == "default_e_r_polygons"


def test_namespace_discovered(tmp_path):
    assert read_namespace(_store(tmp_path, shapes=[f"default_{EXP}_{REG}_polygons", "qc_mask"])) == f"{EXP}_{REG}"


def test_namespace_nameless_store(tmp_path):
    assert read_namespace(_store(tmp_path, shapes=["default_polygons"])) == ""


def test_namespace_ambiguous(tmp_path):
    with pytest.raises(ValueError, match="candidates"):
        read_namespace(_store(tmp_path, shapes=["default_a_polygons", "default_b_polygons"]))


def test_units(tmp_path):
    cfg = _cfg(tmp_path)
    for r in ("region_2", "region_1"):
        (cfg.zarr_store / EXP / r).mkdir(parents=True)
        (cfg.root_path / EXP / "out" / r).mkdir(parents=True)
    (cfg.segmentation_store / EXP / "region_1" / "cellpose_cell").mkdir(parents=True)
    (cfg.segmentation_store / "legacy" / "cellpose_cell" / "region_3").mkdir(parents=True)
    lay = TreeLayout(cfg)
    assert lay.units(EXP) == ["region_1", "region_2"]
    assert lay.units(EXP, under="raw") == ["region_1", "region_2"]
    assert lay.units(EXP, under="seg", label="cellpose_cell") == ["region_1"]
    assert lay.units("legacy", under="seg", label="cellpose_cell") == ["region_3"]


def test_missing_locations_raise():
    d = Dataset(zarr=Path("/z"))
    with pytest.raises(ValueError, match="seg_path"):
        d.seg_dir("cellpose_cell")
    with pytest.raises(ValueError, match="raw"):
        d.raw_images()


def test_seg_legacy_fallback_on_read(tmp_path):
    cfg = _cfg(tmp_path)
    legacy = cfg.segmentation_store / EXP / "cellpose_cell" / REG
    legacy.mkdir(parents=True)
    d = dataset_from(cfg, EXP, REG)
    assert d.seg_dir("cellpose_cell", must_exist=True) == legacy
    assert d.seg_dir("cellpose_cell") == cfg.segmentation_store / EXP / REG / "cellpose_cell"
    assert dataset_from(cfg, EXP, REG, seg="/segs").seg_dir("cellpose_cell", must_exist=True) == Path("/segs/cellpose_cell")


def test_names_checked_against_store(tmp_path):
    z = _store(tmp_path, {"spida": {"namespace": f"{EXP}_{REG}"}})
    assert dataset_from(SpidaConfig(), EXP, REG, zarr=z).namespace == f"{EXP}_{REG}"
    with pytest.raises(ValueError, match="ingested as"):
        dataset_from(SpidaConfig(), EXP, "region_other", zarr=z).namespace


def test_names_on_nameless_store_raise(tmp_path):
    z = _store(tmp_path, shapes=["default_polygons"])
    with pytest.raises(ValueError, match="ingested as ''"):
        dataset_from(SpidaConfig(), EXP, REG, zarr=z).namespace


def test_names_on_new_store_ok(tmp_path):
    assert dataset_from(SpidaConfig(), EXP, REG, zarr=tmp_path / "new.zarr").namespace == f"{EXP}_{REG}"
