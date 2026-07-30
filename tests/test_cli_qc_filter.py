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
            # seg_dir became an OPTION when the storage layout moved to
            # {store}/{exp}/{reg}/{label}; passing it positionally is a parse error.
            "--seg_dir",
            str(seg_dir),
            "--transcript-qc",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["exp_name"] == "EXP"
    assert captured["reg_name"] == "region_001"
    assert captured["seg_dir"] == str(seg_dir)
    # default is now the tissue-mask pass region, not the legacy hex QC shapes
    assert captured["qc_shapes_key"] == "qc_mask"


def test_load_segmentation_region_legacy_qc_shapes_key(monkeypatch, tmp_path: Path):
    """The legacy transcript-QC hex mask stays reachable by naming it explicitly."""
    captured = {}

    monkeypatch.setattr(
        "spida.S.io.main.load_segmentation_region",
        lambda **kwargs: captured.update(kwargs),
    )

    cfg = tmp_path / ".env"
    cfg.write_text("\n")

    result = CliRunner().invoke(
        cli_module.cli,
        [
            "--config",
            str(cfg),
            "load-segmentation-region",
            "EXP",
            "region_001",
            "--transcript-qc",
            "--qc_shapes_key",
            "transcript_qc_shapes",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["qc_shapes_key"] == "transcript_qc_shapes"


def test_load_segmentation_region_rejects_positional_seg_dir(monkeypatch, tmp_path: Path):
    """Regression guard: the SLURM templates passed seg_dir positionally, which
    now fails to parse. Pin the behaviour so the fix isn't silently undone."""
    monkeypatch.setattr(
        "spida.S.io.main.load_segmentation_region",
        lambda **kwargs: None,
    )

    cfg = tmp_path / ".env"
    cfg.write_text("\n")

    result = CliRunner().invoke(
        cli_module.cli,
        [
            "--config", str(cfg),
            "load-segmentation-region",
            "EXP", "region_001", str(tmp_path / "seg"),
        ],
    )

    assert result.exit_code != 0
    assert "unexpected extra argument" in result.output.lower()
