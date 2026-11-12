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
    (tmp_path / "z" / "E" / "region_1" / "tables" / "cellpose_cell_table_filt").mkdir(parents=True)
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


# ---- filter-cells --overwrite (phase 2) ----

def test_filter_cells_guard(monkeypatch, tmp_path):
    import json
    from spida.P import main as m
    from spida.io.layout import Dataset

    d = Dataset(zarr=tmp_path / "s.zarr")
    obs = d.zarr / "tables" / "cellpose_cell_table" / "obs"
    obs.mkdir(parents=True)
    (obs / "zarr.json").write_text(json.dumps({"attributes": {"column-order": ["CELL_ID", "pass_qc"]}}))
    (d.zarr / "zarr.json").write_text(json.dumps({"attributes": {"spida": {"namespace": "", "exp": None, "reg": None}}}))
    read = []
    monkeypatch.setattr(m, "read_table", lambda *a: read.append(a) or (_ for _ in ()).throw(RuntimeError("reached")))
    with pytest.raises(FileExistsError, match="pass_qc"):
        m.filter_cells_region(d, "cellpose_cell")
    assert not read  # stops before reading the table
    with pytest.raises(RuntimeError, match="reached"):
        m.filter_cells_region(d, "cellpose_cell", overwrite=True)


def test_filter_cells_all_skips_filtered(monkeypatch, cfg, tmp_path):
    import json
    from spida.P import main as m
    from spida.io.layout import element_keys
    from spida._constants import TABLE_KEY

    for r in ("region_1", "region_2"):
        (tmp_path / "z" / "E" / r).mkdir(parents=True)
    obs = tmp_path / "z" / "E" / "region_1" / "tables" / element_keys("cellpose_cell", "E_region_1")[TABLE_KEY] / "obs"
    obs.mkdir(parents=True)
    (obs / "zarr.json").write_text(json.dumps({"attributes": {"column-order": ["pass_qc"]}}))
    seen = []
    monkeypatch.setattr(m, "filter_cells_region", lambda d, p, **k: seen.append((d.reg, k["overwrite"])))
    res = CliRunner().invoke(p_cli.cli, ["--config", str(cfg), "filter-cells-all", "cellpose_cell", "--exp", "E"])
    assert res.exit_code == 0, res.output
    assert seen == [("region_2", False)]


@pytest.mark.parametrize("cmd", [["filter-cells-region", "cellpose_cell", "--exp", "E", "--reg", "region_1"],
                                 ["filter-cells-all", "cellpose_cell", "--exp", "E"]])
def test_seg_fam_accepted_ignored_and_warned(monkeypatch, cfg, tmp_path, caplog, cmd):
    # --seg_fam is deprecated: still parses (old job scripts), never reaches P.main, logs a warning
    (tmp_path / "z" / "E" / "region_1").mkdir(parents=True)
    res, calls = _run(monkeypatch, "spida.P.main." + cmd[0].replace("-", "_"), cfg, *cmd, "--seg_fam", "proseg")
    assert res.exit_code == 0, res.output
    assert calls and all("seg_fam" not in k for _, k in calls)
    assert "--seg_fam is ignored" in caplog.text and "process-custom-segmentation" in caplog.text
