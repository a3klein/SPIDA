import json

import pytest
from click.testing import CliRunner

from spida import site


@pytest.fixture
def cfg(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({"ZARR_STORAGE_PATH": str(tmp_path / "z"), "SPIDA_SITE_DIR": str(tmp_path / "site")}))
    return p


def test_seg_qc_figs_names_and_site_dir(monkeypatch, cfg, tmp_path):
    seen = {}

    def fake_load(dataset, site_dir, **k):
        seen.update(d=dataset, site=site_dir, **k)
        raise SystemExit(0)

    monkeypatch.setattr(site, "_load_metadata", fake_load)
    res = CliRunner().invoke(site.cli, ["--config", str(cfg), "generate-seg-qc-figs", "cellpose_cell",
                                        "--exp", "E", "--reg", "region_1"])
    assert res.exit_code == 0, res.output
    assert seen["d"].zarr == tmp_path / "z" / "E" / "region_1" and seen["site"] == tmp_path / "site"
    assert seen["prefix_name"] == "cellpose_cell"


def test_load_metadata_needs_names(tmp_path):
    from spida.io.layout import Dataset
    from spida.utilities.site_utils import _load_metadata

    with pytest.raises(ValueError, match="--exp and --reg"):
        _load_metadata(Dataset(zarr=tmp_path / "z"), tmp_path / "site")
