from pathlib import Path

from click.testing import CliRunner

from spida.S import cli as cli_module


def _invoke(monkeypatch, tmp_path, *args):
    captured = {}

    def _fake(dataset, **kwargs):
        captured.update(kwargs, dataset=dataset)

    monkeypatch.setattr("spida.S.io.main.load_segmentation_region", _fake)
    cfg = tmp_path / "cfg.json"
    cfg.write_text("{}")
    result = CliRunner().invoke(cli_module.cli, ["--config", str(cfg), "load-segmentation-region", *args])
    return result, captured


def test_load_segmentation_region_transcript_qc_flag(monkeypatch, tmp_path: Path):
    seg = tmp_path / "seg"
    result, captured = _invoke(monkeypatch, tmp_path, "--exp", "EXP", "--reg", "region_001",
                               "--zarr_path", str(tmp_path / "z"), "--seg_path", str(seg), "--transcript-qc")
    assert result.exit_code == 0, result.output
    d = captured["dataset"]
    assert (d.exp, d.reg, d.seg) == ("EXP", "region_001", seg)
    # default is the tissue-mask pass region, not the legacy hex QC shapes
    assert captured["qc_shapes_key"] == "qc_mask"


def test_load_segmentation_region_legacy_qc_shapes_key(monkeypatch, tmp_path: Path):
    """The legacy transcript-QC hex mask stays reachable by naming it explicitly."""
    result, captured = _invoke(monkeypatch, tmp_path, "--exp", "EXP", "--reg", "region_001",
                               "--zarr_path", str(tmp_path / "z"),
                               "--transcript-qc", "--qc_shapes_key", "transcript_qc_shapes")
    assert result.exit_code == 0, result.output
    assert captured["qc_shapes_key"] == "transcript_qc_shapes"


def test_load_segmentation_region_rejects_positional_names(monkeypatch, tmp_path: Path):
    """Names are options now; positional EXP/REG must fail to parse, not be silently ignored."""
    result, _ = _invoke(monkeypatch, tmp_path, "EXP", "region_001")
    assert result.exit_code != 0
    assert "unexpected extra argument" in result.output.lower()
