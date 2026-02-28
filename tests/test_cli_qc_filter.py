from pathlib import Path

from click.testing import CliRunner

from spida.S import cli as cli_module


def test_load_segmentation_region_transcript_qc_flag(monkeypatch, tmp_path: Path):
    captured = {}

    def _fake_load_segmentation_region(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(
        "spida.S.io.main.load_segmentation_region",
        _fake_load_segmentation_region,
    )

    cfg = tmp_path / ".env"
    cfg.write_text("\n")
    seg_dir = tmp_path / "seg"
    seg_dir.mkdir(parents=True, exist_ok=True)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "--config",
            str(cfg),
            "load-segmentation-region",
            "EXP",
            "region_001",
            str(seg_dir),
            "--transcript-qc",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["exp_name"] == "EXP"
    assert captured["reg_name"] == "region_001"
    assert captured["qc_shapes_key"] == "transcript_qc_shapes"
