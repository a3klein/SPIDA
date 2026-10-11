import json

import pytest
from click.testing import CliRunner

from spida.P import cli as p_cli


@pytest.fixture
def cfg(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({k: str(tmp_path / v) for k, v in {
        "ZARR_STORAGE_PATH": "z", "PROCESSED_ROOT_PATH": "raw",
        "IMAGE_STORE_PATH": "img", "ANNDATA_STORE_PATH": "ad"}.items()}))
    return p


def _run(monkeypatch, target, cfg, *args):
    calls = []
    monkeypatch.setattr(target, lambda *a, **k: calls.append((a, k)))
    res = CliRunner().invoke(p_cli.cli, ["--config", str(cfg), *args])
    return res, calls


def test_filter_region_names(monkeypatch, cfg, tmp_path):
    res, calls = _run(monkeypatch, "spida.P.main.filter_cells_region", cfg,
                      "filter-cells-region", "cellpose_cell", "--exp", "E", "--reg", "region_1")
    assert res.exit_code == 0, res.output
    (d, prefix), _ = calls[0]
    assert prefix == "cellpose_cell" and d.zarr == tmp_path / "z" / "E" / "region_1"
    assert d.figures == tmp_path / "img" / "E" / "region_1"


def test_setup_region_path_mode(monkeypatch, cfg, tmp_path):
    res, calls = _run(monkeypatch, "spida.P.main.setup_adata_region", cfg,
                      "setup-adata-region", "cellpose_cell", "--zarr_path", str(tmp_path / "x.zarr"),
                      "--figures_path", str(tmp_path / "f"))
    assert res.exit_code == 0, res.output
    (d, _), _ = calls[0]
    assert not d.named and d.zarr == tmp_path / "x.zarr" and d.figures == tmp_path / "f"


def test_filter_all_honours_zarr_store(monkeypatch, cfg, tmp_path):
    mine = tmp_path / "mine"
    for r in ("region_a", "region_b"):
        (mine / "E" / r).mkdir(parents=True)
    (tmp_path / "z" / "E" / "region_config").mkdir(parents=True)
    res, calls = _run(monkeypatch, "spida.P.main.filter_cells_region", cfg,
                      "filter-cells-all", "cellpose_cell", "--exp", "E", "--zarr_store", str(mine))
    assert res.exit_code == 0, res.output
    assert [a[0].zarr for a, _ in calls] == [mine / "E" / "region_a", mine / "E" / "region_b"]


def test_call_region_tz_positionals(monkeypatch, cfg):
    res, calls = _run(monkeypatch, "spida.P.main.call_region_tz", cfg,
                      "call-region-tz", "default", "wm_regions", "--exp", "E", "--reg", "region_1")
    assert res.exit_code == 0, res.output
    (d,), k = calls[0]
    assert k["prefix_name"] == "default" and k["geoms_name"] == "wm_regions" and d.named


def test_write_adata_defaults_to_anndata_store(monkeypatch, cfg, tmp_path):
    seen = []
    monkeypatch.setattr("spida.P.main._write_adata", lambda zarr, key, out: seen.append((zarr, key, out)))
    res = CliRunner().invoke(p_cli.cli, ["--config", str(cfg), "write-adata", "--exp", "E", "--reg", "region_1",
                                         "--prefix_names", "cellpose_cell", "--suffix", "_filt"])
    assert res.exit_code == 0, res.output
    zarr, key, out = seen[0]
    assert zarr == tmp_path / "z" / "E" / "region_1"
    assert key == "cellpose_cell_table_filt"
    assert out == tmp_path / "ad" / "E" / "region_1" / "cellpose_cell_table_filt.h5ad"


def test_tissue_mask_output_path_override(monkeypatch, cfg, tmp_path):
    res, calls = _run(monkeypatch, "spida.P.main.tissue_mask_region", cfg,
                      "tissue-mask", "--exp", "E", "--reg", "region_1", "--tissue_mask_output_path", str(tmp_path / "masks"))
    assert res.exit_code == 0, res.output
    (d,), k = calls[0]
    assert k["mask_dir"] == str(tmp_path / "masks") and d.raw == tmp_path / "raw" / "E" / "out" / "region_1"
