"""Tests for ``tissue_mask_path`` on ``segment_region``.

The tissue mask is a property of the REGION, not of the segmenter, so it is a
parameter of ``segment_region`` rather than a field of any ``{Method}Config``.
These tests pin that wiring: where the mask is discovered, when its absence is
tolerated, and when a bad path must fail loudly instead of silently segmenting
the whole mosaic.

Free of cellpose and of any GPU: the backend call is stubbed out.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from spida.S.segmentation import main as seg_main


@pytest.fixture(autouse=True)
def _any_pixi_env(monkeypatch):
    """Drop the env guard: this suite tests wiring, not which env it runs in.

    ``spec.require_env`` no-ops when ``PIXI_ENVIRONMENT_NAME`` is unset, so the
    same tests run in either env.
    """
    monkeypatch.delenv("PIXI_ENVIRONMENT_NAME", raising=False)


@pytest.fixture
def region(tmp_path) -> tuple[Path, Path]:
    """``(root, images_dir)`` for a minimal EXP/region layout."""
    images = tmp_path / "root" / "EXP" / "out" / "REG" / "images"
    images.mkdir(parents=True)
    (images / "mosaic_DAPI_z0.tif").write_bytes(b"")                     # inputs segment-region requires
    (images.parent / "detected_transcripts.csv").write_text("")
    return tmp_path / "root", images


# ---------------------------------------------------------------------------
# resolve_tissue_mask_path
# ---------------------------------------------------------------------------

def test_auto_discovers_the_conventional_filename(region):
    _root, images = region
    mask = images / seg_main.TISSUE_MASK_FILENAME
    mask.write_bytes(b"")
    assert seg_main.resolve_tissue_mask_path(None, images) == mask


def test_absent_mask_is_not_an_error(region, caplog):
    """A region that has never had `tissue-mask` run must still segment."""
    _root, images = region
    with caplog.at_level("INFO"):
        assert seg_main.resolve_tissue_mask_path(None, images) is None
    assert "no tissue mask" in caplog.text


def test_explicit_missing_path_raises(region):
    """A typo must not degrade into segmenting the entire mosaic."""
    _root, images = region
    with pytest.raises(FileNotFoundError) as exc:
        seg_main.resolve_tissue_mask_path(images / "typo.tif", images)
    assert "tissue-mask" in str(exc.value)


def test_explicit_directory_gets_the_filename_appended(region):
    _root, images = region
    other = images.parent / "elsewhere"
    other.mkdir()
    mask = other / seg_main.TISSUE_MASK_FILENAME
    mask.write_bytes(b"")
    assert seg_main.resolve_tissue_mask_path(other, images) == mask


def test_explicit_path_wins_over_the_conventional_one(region):
    _root, images = region
    (images / seg_main.TISSUE_MASK_FILENAME).write_bytes(b"")
    custom = images.parent / "custom_mask.tif"
    custom.write_bytes(b"")
    assert seg_main.resolve_tissue_mask_path(custom, images) == custom


def test_user_expansion_is_applied(region, monkeypatch):
    _root, images = region
    mask = images / seg_main.TISSUE_MASK_FILENAME
    mask.write_bytes(b"")
    monkeypatch.setenv("HOME", str(images))
    assert seg_main.resolve_tissue_mask_path(
        f"~/{seg_main.TISSUE_MASK_FILENAME}", images) == mask


# ---------------------------------------------------------------------------
# segment_region wiring
# ---------------------------------------------------------------------------

@pytest.fixture
def captured(monkeypatch) -> dict:
    """Stub the cellpose backend and capture exactly what it was handed.

    Injected as a module rather than monkeypatched on the real one: importing
    ``backends.cellpose`` pulls in cellpose + torch, which the preprocessing env
    does not have and which this wiring does not need.
    """
    import sys
    import types

    calls: dict = {}

    def fake_run(root_dir, output_dir, region, **kwargs):
        calls.update(root_dir=root_dir, output_dir=output_dir, region=region,
                     kwargs=kwargs)
        return False

    stub = types.ModuleType("spida.S.segmentation.backends.cellpose")
    stub.run_cellpose = fake_run
    monkeypatch.setitem(sys.modules, stub.__name__, stub)
    return calls


def _dataset(root):
    from spida.config import SpidaConfig
    from spida.io.layout import dataset_from
    return dataset_from(SpidaConfig(root_path=Path(root), segmentation_store=Path(root).parent / "seg"), "EXP", "REG")


def _segment(root, **kw):
    return seg_main.segment_region("cellpose", _dataset(root), **kw)


def test_segment_region_forwards_the_discovered_mask(region, captured):
    root, images = region
    mask = images / seg_main.TISSUE_MASK_FILENAME
    mask.write_bytes(b"")
    _segment(root)
    assert captured["kwargs"]["tissue_mask_path"] == mask


def test_segment_region_forwards_none_when_absent(region, captured):
    root, _images = region
    _segment(root)
    assert captured["kwargs"]["tissue_mask_path"] is None


def test_segment_region_forwards_an_explicit_mask(region, captured):
    root, images = region
    custom = images.parent / "custom.tif"
    custom.write_bytes(b"")
    _segment(root, tissue_mask_path=str(custom))
    assert captured["kwargs"]["tissue_mask_path"] == custom


def test_segment_region_reads_images_from_the_region_dir(region, captured):
    root, images = region
    _segment(root)
    assert Path(captured["root_dir"]) == images


def test_method_kwargs_still_pass_through(region, captured):
    root, _images = region
    _segment(root, gpu_tile_size=9000, seg_mode="2d")
    assert captured["kwargs"]["gpu_tile_size"] == 9000
    assert captured["kwargs"]["seg_mode"] == "2d"


def test_mask_rejected_for_a_backend_that_cannot_use_it(region, monkeypatch):
    """Silently dropping it would mean paying for a mask that does nothing."""
    root, images = region
    (images / seg_main.TISSUE_MASK_FILENAME).write_bytes(b"")
    with pytest.raises(ValueError, match="does not consume it"):
        seg_main.segment_region(
            "proseg", _dataset(root),
            tissue_mask_path=str(images / seg_main.TISSUE_MASK_FILENAME),
        )


def test_proseg_unaffected_when_no_mask_is_requested(region, monkeypatch):
    """Auto-discovery must not leak into backends that do not accept the kwarg."""
    from spida.S.segmentation.backends import proseg as pb

    root, images = region
    (images / seg_main.TISSUE_MASK_FILENAME).write_bytes(b"")
    seen: dict = {}
    monkeypatch.setattr(pb, "run_proseg",
                        lambda *a, **k: seen.update(args=a, kwargs=k))
    seg_main.segment_region("proseg", _dataset(root))
    assert "tissue_mask_path" not in seen["kwargs"]
