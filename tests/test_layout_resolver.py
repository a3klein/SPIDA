"""Tests for resolve_region_dir: the single source of truth for the segmentation
output layout. Current layout is {store}/{exp}/{reg}/{label}; a single legacy
fallback {store}/{exp}/{label}/{reg} keeps pre-reorder data readable.
"""
from spida.S.segmentation.backends.base import resolve_region_dir

EXP, REG, LABEL = "exp1", "region_001", "cellpose_cell"


def test_write_returns_current_layout(tmp_path):
    # writes (must_exist=False) always target the current {exp}/{reg}/{label} layout
    out = resolve_region_dir(tmp_path, EXP, REG, LABEL)
    assert out == tmp_path / EXP / REG / LABEL


def test_read_prefers_current(tmp_path):
    current = tmp_path / EXP / REG / LABEL
    legacy = tmp_path / EXP / LABEL / REG
    current.mkdir(parents=True)
    legacy.mkdir(parents=True)
    assert resolve_region_dir(tmp_path, EXP, REG, LABEL, must_exist=True) == current


def test_read_falls_back_to_legacy(tmp_path):
    legacy = tmp_path / EXP / LABEL / REG
    legacy.mkdir(parents=True)                       # only the legacy layout exists
    assert resolve_region_dir(tmp_path, EXP, REG, LABEL, must_exist=True) == legacy


def test_read_defaults_to_current_when_absent(tmp_path):
    # neither exists -> return current (error surfaces downstream, not here)
    out = resolve_region_dir(tmp_path, EXP, REG, LABEL, must_exist=True)
    assert out == tmp_path / EXP / REG / LABEL


def test_label_distinguishes_runs(tmp_path):
    # cellpose_nuc vs cellpose_cell land in separate dirs under the same region
    nuc = resolve_region_dir(tmp_path, EXP, REG, "cellpose_nuc")
    cell = resolve_region_dir(tmp_path, EXP, REG, "cellpose_cell")
    assert nuc.parent == cell.parent == tmp_path / EXP / REG
    assert nuc != cell
