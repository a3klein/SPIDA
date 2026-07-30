"""Tests for `segment-region <method> --list-params` and the config type check.

Method-specific run parameters are deliberately NOT Click options (three
backends would make the shared command unusable), so `--list-params` is the only
way a user discovers them. It therefore needs to work without a region, without
a GPU, and from any pixi env.
"""

from __future__ import annotations

from click.testing import CliRunner
import pytest

from spida.S import cli as cli_module
from spida.S.segmentation.backends import (
    describe_method_params,
    get_config_class,
)
from spida.S.segmentation.backends.cellpose_config import CellposeConfig
from spida.S.segmentation.backends.config import ConfigError


# ---------------------------------------------------------------------------
# dispatcher
# ---------------------------------------------------------------------------

def test_get_config_class_for_cellpose():
    assert get_config_class("cellpose") is CellposeConfig


def test_get_config_class_none_for_methods_without_one():
    assert get_config_class("proseg") is None


def test_describe_cellpose_lists_parameters():
    text = describe_method_params("cellpose")
    assert "CellposeConfig" in text
    for name in ("seg_mode", "model", "use_bfloat16", "anisotropy", "poly_tile_size"):
        assert name in text


def test_describe_documents_boolean_convention():
    """parse_click_kwargs needs --flag=true; a bare --flag eats the next token."""
    text = describe_method_params("cellpose")
    assert "--apply_clahe=true" in text


def test_describe_known_method_without_config_is_graceful():
    text = describe_method_params("proseg")
    assert "does not have a run-configuration class yet" in text
    assert "cellpose" in text                      # points at what does


def test_describe_unknown_method_lists_known_ones():
    text = describe_method_params("cellpose2")
    assert "Unknown segmentation method" in text
    assert "cellpose" in text and "proseg" in text


# ---------------------------------------------------------------------------
# CLI wiring
# ---------------------------------------------------------------------------

def _invoke(args, tmp_path):
    cfg = tmp_path / ".env"
    cfg.write_text("\n")
    return CliRunner().invoke(cli_module.cli, ["--config", str(cfg), *args])


def test_list_params_needs_no_region(tmp_path):
    result = _invoke(["segment-region", "cellpose", "--list-params"], tmp_path)
    assert result.exit_code == 0, result.output
    assert "seg_mode" in result.output


def test_list_params_underscore_spelling(tmp_path):
    result = _invoke(["segment-region", "cellpose", "--list_params"], tmp_path)
    assert result.exit_code == 0, result.output
    assert "seg_mode" in result.output


def test_list_params_does_not_run_segmentation(monkeypatch, tmp_path):
    called = []
    monkeypatch.setattr(
        "spida.S.segmentation.main.segment_region",
        lambda **kw: called.append(kw),
    )
    result = _invoke(["segment-region", "cellpose", "--list-params"], tmp_path)
    assert result.exit_code == 0
    assert called == []


def test_missing_region_without_list_params_is_a_usage_error(tmp_path):
    result = _invoke(["segment-region", "cellpose", "EXP"], tmp_path)
    assert result.exit_code != 0
    assert "REG_NAME are required" in result.output


def test_normal_invocation_still_works(monkeypatch, tmp_path):
    """--list-params must not disturb the ordinary path."""
    called = {}
    monkeypatch.setattr(
        "spida.S.segmentation.main.segment_region",
        lambda **kw: called.update(kw),
    )
    result = _invoke(
        ["segment-region", "cellpose", "EXP", "region_001", "--scale=4"],
        tmp_path,
    )
    assert result.exit_code == 0, result.output
    assert called["method"] == "cellpose"
    assert called["exp_name"] == "EXP"
    assert called["reg_name"] == "region_001"
    assert called["scale"] == 4


# ---------------------------------------------------------------------------
# type checking (the reason the boolean convention matters)
# ---------------------------------------------------------------------------

def test_bool_field_rejects_a_string():
    """`--apply_clahe --scale=4` makes parse_click_kwargs hand us the string
    '--scale=4' as apply_clahe's value; that must not sail through as truthy."""
    with pytest.raises(ConfigError) as exc:
        CellposeConfig.from_kwargs(apply_clahe="--scale=4")
    msg = str(exc.value)
    assert "apply_clahe" in msg
    assert "--apply_clahe=true" in msg               # actionable hint


def test_bool_field_accepts_a_bool():
    assert CellposeConfig.from_kwargs(apply_clahe=True).apply_clahe is True


def test_int_field_rejects_a_bool():
    with pytest.raises(ConfigError, match="is a bool"):
        CellposeConfig.from_kwargs(scale=True)


def test_int_field_rejects_a_non_numeric_string():
    with pytest.raises(ConfigError, match="expected int"):
        CellposeConfig.from_kwargs(scale="lots")


def test_float_field_accepts_an_int():
    cfg = CellposeConfig.from_kwargs(micron_per_z=2)
    assert cfg.micron_per_z == 2


def test_optional_field_accepts_none():
    assert CellposeConfig.from_kwargs(diameter=None).diameter is None


def test_str_field_rejects_a_number():
    with pytest.raises(ConfigError, match="expected str"):
        CellposeConfig.from_kwargs(nuc_stain_name=3)


def test_literal_union_field_accepts_both_forms():
    assert CellposeConfig.from_kwargs(seg_mode="3d_true", anisotropy="auto").anisotropy == "auto"
    assert CellposeConfig.from_kwargs(seg_mode="3d_true", anisotropy=3.5).anisotropy == 3.5


def test_choices_field_keeps_its_specific_message():
    """CHOICES fields are skipped by the type check so validate() reports them."""
    with pytest.raises(ConfigError, match="is not one of"):
        CellposeConfig.from_kwargs(seg_mode="3d")
