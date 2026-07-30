"""Tests for CellposeConfig and cellpose model resolution.

Deliberately free of any cellpose import and of any GPU requirement, so this
runs on the head node in either pixi env.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from spida.S.segmentation.backends.config import ConfigError
from spida.S.segmentation.backends import cellpose_config as cc
from spida.S.segmentation.backends.cellpose_config import (
    CellposeConfig,
    KNOWN_MODELS,
    SPIDA_ENV_VAR,
    CELLPOSE_ENV_VAR,
    derive_anisotropy,
    read_native_xy_um,
    resolve_model,
    resolve_model_dir,
)


# ---------------------------------------------------------------------------
# construction + defaults
# ---------------------------------------------------------------------------

def test_defaults_match_production():
    cfg = CellposeConfig.from_kwargs()
    assert cfg.seg_mode == "3d_stitched"      # current SPIDA 3a behaviour
    assert cfg.model == "cpsam"
    assert cfg.scale == 4
    assert cfg.stitch_threshold == 0.2
    assert cfg.min_size == 100
    assert cfg.min_z == 3
    assert cfg.use_bfloat16 is True           # deliberate change from fp32


def test_production_template_kwargs_are_accepted():
    """The exact parameter set the 3a SLURM template passes today."""
    cfg = CellposeConfig.from_kwargs(
        scale=4,
        image_ext="_z3.decon.tif",
        nuc_stain_name="DAPI",
        cyto_stain_name="PolyT",
        flow_threshold=0,
        cellprob_threshold=-2,
        tile_norm_blocksize=2960,
    )
    assert cfg.cyto_stain_name == "PolyT"
    assert cfg.tile_norm_blocksize == 2960
    assert cfg.flow_threshold == 0


# ---------------------------------------------------------------------------
# unknown keys and aliases
# ---------------------------------------------------------------------------

def test_unknown_parameter_raises_with_suggestion():
    with pytest.raises(ConfigError) as exc:
        CellposeConfig.from_kwargs(anisotrpy=3.5)
    msg = str(exc.value)
    assert "unknown parameter 'anisotrpy'" in msg
    assert "anisotropy" in msg          # close-match suggestion
    assert "--list-params" in msg


def test_unknown_parameter_without_close_match_still_raises():
    with pytest.raises(ConfigError, match="unknown parameter"):
        CellposeConfig.from_kwargs(completely_made_up_thing=1)


@pytest.mark.parametrize(
    "alias,target,value",
    [
        ("model_name", "model", "cpsam"),
        ("tile_size", "poly_tile_size", 512),
        ("overlap", "poly_tile_overlap", 64),
        ("clahe_clipLimit", "clahe_clip_limit", 5.0),
        ("clahe_tileGridSize", "clahe_tile_grid_size", 16),
    ],
)
def test_deprecated_aliases_map_and_warn(alias, target, value):
    with pytest.warns(DeprecationWarning):
        cfg = CellposeConfig.from_kwargs(**{alias: value})
    assert getattr(cfg, target) == value


def test_alias_and_canonical_together_raises():
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ConfigError, match="supplied twice"):
            CellposeConfig.from_kwargs(poly_tile_size=512, tile_size=256)


# ---------------------------------------------------------------------------
# legacy mode flags -> seg_mode
# ---------------------------------------------------------------------------

def test_do_3d_true_maps_to_3d_true():
    with pytest.warns(DeprecationWarning, match="do_3D"):
        cfg = CellposeConfig.from_kwargs(do_3D=True)
    assert cfg.seg_mode == "3d_true"


def test_do_3d_false_maps_to_stitched():
    with pytest.warns(DeprecationWarning, match="do_3D"):
        cfg = CellposeConfig.from_kwargs(do_3D=False)
    assert cfg.seg_mode == "3d_stitched"


def test_project_3d_to_2d_maps_to_2d_mip():
    with pytest.warns(DeprecationWarning, match="project_3d_to_2d"):
        cfg = CellposeConfig.from_kwargs(project_3d_to_2d=True)
    assert cfg.seg_mode == "2d"
    assert cfg.z_reduce == "mip"        # legacy flag did a max projection


def test_legacy_flag_with_seg_mode_raises():
    with pytest.raises(ConfigError, match="cannot be combined"):
        CellposeConfig.from_kwargs(seg_mode="3d_true", do_3D=True)


# ---------------------------------------------------------------------------
# validation rules
# ---------------------------------------------------------------------------

def test_invalid_seg_mode_rejected():
    with pytest.raises(ConfigError, match="seg_mode"):
        CellposeConfig.from_kwargs(seg_mode="3d")


def test_invalid_model_rejected():
    with pytest.raises(ConfigError, match="model"):
        CellposeConfig.from_kwargs(model="cellpose3")


def test_anisotropy_rejected_outside_3d_true():
    """Silently ignoring a supplied anisotropy is exactly the failure mode
    this config layer exists to prevent."""
    with pytest.raises(ConfigError, match="only used when seg_mode='3d_true'"):
        CellposeConfig.from_kwargs(seg_mode="3d_stitched", anisotropy=3.5)


def test_anisotropy_allowed_for_3d_true():
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true", anisotropy=3.5)
    assert cfg.anisotropy == 3.5


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (dict(scale=0), "scale must be >= 1"),
        (dict(n_z_planes=0), "n_z_planes must be >= 1"),
        (dict(micron_per_z=0), "micron_per_z must be > 0"),
        (dict(poly_tile_size=100, poly_tile_overlap=100), "must be smaller"),
        (dict(min_z=9, n_z_planes=7), "exceeds n_z_planes"),
        (dict(seg_mode="3d_true", anisotropy=-1), "anisotropy must be > 0"),
    ],
)
def test_validation_rules(kwargs, match):
    with pytest.raises(ConfigError, match=match):
        CellposeConfig.from_kwargs(**kwargs)


def test_string_values_are_coerced():
    cfg = CellposeConfig.from_kwargs(scale="4", cellprob_threshold="-2.5")
    assert cfg.scale == 4 and isinstance(cfg.scale, int)
    assert cfg.cellprob_threshold == -2.5


# ---------------------------------------------------------------------------
# presentation / provenance
# ---------------------------------------------------------------------------

def test_describe_lists_every_field_and_is_readable():
    text = CellposeConfig.describe()
    for name in CellposeConfig.field_names():
        assert name in text
    assert "seg_mode" in text and "3d_stitched" in text
    assert "Deprecated aliases" in text
    assert "tile_size -> poly_tile_size" in text


def test_to_meta_is_json_serialisable():
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true")
    meta = cfg.to_meta()
    json.dumps(meta)                       # must not raise
    assert meta["seg_mode"] == "3d_true"
    assert meta["model"] == "cpsam"


# ---------------------------------------------------------------------------
# model directory resolution
# ---------------------------------------------------------------------------

def test_model_dir_precedence(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    spida = tmp_path / "spida"
    cellpose = tmp_path / "cellpose"
    monkeypatch.setenv(SPIDA_ENV_VAR, str(spida))
    monkeypatch.setenv(CELLPOSE_ENV_VAR, str(cellpose))

    assert resolve_model_dir(explicit) == explicit          # explicit wins
    assert resolve_model_dir(None) == spida                 # then SPIDA's key

    monkeypatch.delenv(SPIDA_ENV_VAR)
    assert resolve_model_dir(None) == cellpose              # then cellpose's

    monkeypatch.delenv(CELLPOSE_ENV_VAR)
    assert resolve_model_dir(None) == Path.home() / ".cellpose" / "models"


def test_spida_env_var_is_actually_honoured(tmp_path, monkeypatch):
    """Regression guard for the original bug: CELLPOSE_MODEL_PATH was defined
    in every config JSON but never reached the model loader."""
    monkeypatch.setenv(SPIDA_ENV_VAR, str(tmp_path))
    (tmp_path / "cpsam").write_bytes(b"weights")
    res = resolve_model("cpsam")
    assert res.exists_locally
    assert res.pretrained_model == str(tmp_path / "cpsam")


# ---------------------------------------------------------------------------
# model resolution
# ---------------------------------------------------------------------------

def test_resolve_model_uses_absolute_path_when_present(tmp_path):
    (tmp_path / "cpsam_v2").write_bytes(b"x" * 1234)
    res = resolve_model("cpsam_v2", model_dir=tmp_path)
    assert res.exists_locally and not res.will_download
    assert res.pretrained_model == str(tmp_path / "cpsam_v2")
    assert res.size_bytes == 1234
    res.log()                               # must not raise


def test_resolve_model_unknown_name_raises(tmp_path):
    with pytest.raises(ConfigError, match="unknown cellpose model"):
        resolve_model("not_a_model", model_dir=tmp_path)


def test_resolve_model_missing_but_fetchable_downloads(tmp_path, monkeypatch):
    monkeypatch.setattr(cc, "_cellpose_known_models", lambda: {"cpsam", "cpsam_v2"})
    res = resolve_model("cpsam_v2", model_dir=tmp_path)
    assert res.will_download and not res.exists_locally
    assert res.pretrained_model == "cpsam_v2"    # bare name, so cellpose fetches
    res.log()


def test_resolve_model_missing_and_unfetchable_raises(tmp_path, monkeypatch):
    """The dangerous case: on cellpose 4.0.9, MODEL_NAMES == ['cpsam'], so a
    request for cpsam_v2 would silently load cpsam. We must refuse instead."""
    monkeypatch.setattr(cc, "_cellpose_known_models", lambda: {"cpsam"})
    with pytest.raises(ConfigError) as exc:
        resolve_model("cpsam_v2", model_dir=tmp_path)
    msg = str(exc.value)
    assert "does not know how to fetch it" in msg
    assert "4.2.1.1" in msg                       # actionable upgrade hint


def test_resolve_model_permits_download_when_cellpose_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(cc, "_cellpose_known_models", lambda: set())
    res = resolve_model("cpdino", model_dir=tmp_path)
    assert res.will_download


def test_every_known_model_has_a_description():
    assert set(KNOWN_MODELS) == set(CellposeConfig.CHOICES["model"])
    assert all(desc.strip() for desc in KNOWN_MODELS.values())


# ---------------------------------------------------------------------------
# anisotropy
# ---------------------------------------------------------------------------

def test_derive_anisotropy_matches_measured_value():
    """xy = 106 nm, z = 1500 nm, scale 4 -> the dev bench's 3.5."""
    assert derive_anisotropy(1.5, 0.106, 4) == pytest.approx(3.538, abs=0.01)


def test_derive_anisotropy_is_coupled_to_scale():
    at4 = derive_anisotropy(1.5, 0.106, 4)
    at2 = derive_anisotropy(1.5, 0.106, 2)
    assert at2 == pytest.approx(at4 * 2, rel=1e-6)


def _write_transform(directory: Path, px_per_um: float) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "micron_to_mosaic_pixel_transform.csv").write_text(
        f"{px_per_um} 0.0 0.0\n0.0 {px_per_um} 0.0\n0.0 0.0 1.0\n"
    )


def test_read_native_xy_um(tmp_path):
    _write_transform(tmp_path, 1 / 0.106)
    assert read_native_xy_um(tmp_path) == pytest.approx(0.106, rel=1e-6)


def test_read_native_xy_um_missing_file_raises(tmp_path):
    with pytest.raises(ConfigError, match="cannot derive anisotropy"):
        read_native_xy_um(tmp_path)


def test_resolved_anisotropy_auto(tmp_path):
    _write_transform(tmp_path, 1 / 0.106)
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true")
    assert cfg.resolved_anisotropy(tmp_path) == pytest.approx(3.538, abs=0.01)


def test_resolved_anisotropy_none_for_other_modes(tmp_path):
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_stitched")
    assert cfg.resolved_anisotropy(tmp_path) is None


def test_resolved_anisotropy_explicit_skips_transform(tmp_path):
    """No transform file present: an explicit value must not need one."""
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true", anisotropy=3.5)
    assert cfg.resolved_anisotropy(tmp_path) == 3.5


# ---------------------------------------------------------------------------
# micron_to_mosaic_path override
# ---------------------------------------------------------------------------

def test_read_native_xy_um_accepts_an_explicit_file(tmp_path):
    _write_transform(tmp_path, 1 / 0.106)
    explicit = tmp_path / "micron_to_mosaic_pixel_transform.csv"
    assert read_native_xy_um(explicit) == pytest.approx(0.106, rel=1e-6)


def test_micron_to_mosaic_path_overrides_auto_discovery(tmp_path):
    """The transform lives somewhere other than the images dir."""
    images = tmp_path / "images"
    images.mkdir()                                  # deliberately has NO transform
    elsewhere = tmp_path / "elsewhere"
    _write_transform(elsewhere, 1 / 0.212)          # 2x coarser pixels

    cfg = CellposeConfig.from_kwargs(
        seg_mode="3d_true",
        micron_to_mosaic_path=str(elsewhere / "micron_to_mosaic_pixel_transform.csv"),
    )
    # 1.5 / (0.212 * 4) ~= 1.77
    assert cfg.resolved_anisotropy(images) == pytest.approx(1.769, abs=0.01)


def test_micron_to_mosaic_path_accepts_a_directory(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    _write_transform(elsewhere, 1 / 0.106)
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true",
                                     micron_to_mosaic_path=str(elsewhere))
    assert cfg.resolved_anisotropy(tmp_path) == pytest.approx(3.538, abs=0.01)


def test_auto_discovery_is_still_the_default(tmp_path):
    """Default None must keep finding <images_dir>/micron_to_mosaic...csv."""
    _write_transform(tmp_path, 1 / 0.106)
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true")
    assert cfg.micron_to_mosaic_path is None
    assert cfg.resolved_anisotropy(tmp_path) == pytest.approx(3.538, abs=0.01)


def test_missing_transform_error_mentions_both_escape_hatches(tmp_path):
    cfg = CellposeConfig.from_kwargs(seg_mode="3d_true")
    with pytest.raises(ConfigError) as exc:
        cfg.resolved_anisotropy(tmp_path)
    msg = str(exc.value)
    assert "--micron_to_mosaic_path" in msg
    assert "--anisotropy" in msg


def test_unused_micron_to_mosaic_path_warns_but_does_not_raise(tmp_path, caplog):
    with caplog.at_level("WARNING"):
        cfg = CellposeConfig.from_kwargs(
            seg_mode="3d_stitched", micron_to_mosaic_path=str(tmp_path),
        )
    assert cfg.micron_to_mosaic_path == str(tmp_path)
    assert "unused" in caplog.text
