"""Tests for the narrow config -> environment bridge.

Regression coverage for a real bug: ``CELLPOSE_MODEL_PATH`` was present in every
per-region config JSON but never reached the model loader, because
``ConfigDefaultGroup`` only installs config values as Click *option* defaults
and no command declares an option for that key. The bridge fixes that for the
handful of keys consumed via ``os.getenv`` — and must NOT touch the rest.
"""

from __future__ import annotations

import json
import os

import pytest

from spida.config import (
    DEFAULT_CONFIG,
    CONFIG_DESCRIPTIONS,
    ENV_ONLY_CONFIG_KEYS,
    RENAME_CONFIG_KEYS,
    THIRD_PARTY_ENV_ALIASES,
    export_env_only_config,
    resolve_config,
)


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for var in ("CELLPOSE_MODEL_PATH", "CELLPOSE_LOCAL_MODELS_PATH"):
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------------------
# registration
# ---------------------------------------------------------------------------

def test_cellpose_model_path_is_a_first_class_key():
    assert "CELLPOSE_MODEL_PATH" in DEFAULT_CONFIG
    assert "CELLPOSE_MODEL_PATH" in CONFIG_DESCRIPTIONS


def test_default_is_usable_not_a_placeholder():
    """The key is optional, so its default must be a real location rather than
    one of the `/path/to/...` placeholders used for mandatory keys."""
    default = DEFAULT_CONFIG["CELLPOSE_MODEL_PATH"]
    assert not default.startswith("/path/to")
    assert default == "~/.cellpose/models"


def test_env_only_keys_have_no_click_option():
    """Anything with a Click option must not be exported — that is precisely the
    blanket-export behaviour this bridge avoids."""
    for key in ENV_ONLY_CONFIG_KEYS:
        assert key not in RENAME_CONFIG_KEYS


# ---------------------------------------------------------------------------
# export behaviour
# ---------------------------------------------------------------------------

def test_export_sets_the_key_and_the_third_party_alias():
    exported = export_env_only_config({"CELLPOSE_MODEL_PATH": "/models"})
    assert os.environ["CELLPOSE_MODEL_PATH"] == "/models"
    assert os.environ["CELLPOSE_LOCAL_MODELS_PATH"] == "/models"
    assert exported["CELLPOSE_MODEL_PATH"] == "/models"


def test_export_does_not_touch_unrelated_keys():
    """The whole point of the narrow bridge: storage paths keep coming from
    .env, not from the per-region JSON."""
    before = dict(os.environ)
    export_env_only_config({
        "CELLPOSE_MODEL_PATH": "/models",
        "ZARR_STORAGE_PATH": "/scratch/data/zarr_store",
        "PROCESSED_ROOT_PATH": "/scratch",
        "VPT_BIN_PATH": "/path/to/vpt/bin",
    })
    changed = {k for k in set(before) | set(os.environ)
               if before.get(k) != os.environ.get(k)}
    assert changed == {"CELLPOSE_MODEL_PATH", "CELLPOSE_LOCAL_MODELS_PATH"}


def test_export_skips_missing_keys():
    assert export_env_only_config({}) == {}
    assert "CELLPOSE_MODEL_PATH" not in os.environ


def test_existing_third_party_alias_is_not_clobbered(monkeypatch, capsys):
    monkeypatch.setenv("CELLPOSE_LOCAL_MODELS_PATH", "/user/choice")
    export_env_only_config({"CELLPOSE_MODEL_PATH": "/models"})
    assert os.environ["CELLPOSE_LOCAL_MODELS_PATH"] == "/user/choice"
    assert os.environ["CELLPOSE_MODEL_PATH"] == "/models"
    assert "differs from" in capsys.readouterr().out       # conflict surfaced


def test_matching_alias_produces_no_warning(monkeypatch, capsys):
    monkeypatch.setenv("CELLPOSE_LOCAL_MODELS_PATH", "/models")
    export_env_only_config({"CELLPOSE_MODEL_PATH": "/models"})
    assert "differs from" not in capsys.readouterr().out


def test_non_string_values_are_coerced():
    """A blanket os.environ[k] = v would raise TypeError here."""
    from pathlib import Path
    export_env_only_config({"CELLPOSE_MODEL_PATH": Path("/models")})
    assert os.environ["CELLPOSE_MODEL_PATH"] == "/models"


def test_every_alias_target_is_an_exported_key():
    assert set(THIRD_PARTY_ENV_ALIASES).issubset(set(ENV_ONLY_CONFIG_KEYS))


# ---------------------------------------------------------------------------
# end-to-end through resolve_config
# ---------------------------------------------------------------------------

def test_json_config_reaches_the_environment(tmp_path):
    """The end-to-end path that was broken: a value in a per-region config JSON
    must be visible to the cellpose model loader."""
    cfg = tmp_path / "BR_config.json"
    cfg.write_text(json.dumps({
        "ZARR_STORAGE_PATH": "/scratch/data/zarr_store",
        "CELLPOSE_MODEL_PATH": "/shared/weights",
    }))

    resolved = resolve_config(config_path=str(cfg), defaults=DEFAULT_CONFIG)
    export_env_only_config(resolved)

    from spida.S.segmentation.backends.cellpose_config import resolve_model_dir
    assert str(resolve_model_dir()) == "/shared/weights"
    # unrelated key stayed out of the environment
    assert os.environ.get("ZARR_STORAGE_PATH") != "/scratch/data/zarr_store"


def test_absent_key_falls_back_to_the_usable_default(tmp_path):
    cfg = tmp_path / "BR_config.json"
    cfg.write_text(json.dumps({"ZARR_STORAGE_PATH": "/scratch/data/zarr_store"}))

    resolved = resolve_config(config_path=str(cfg), defaults=DEFAULT_CONFIG)
    export_env_only_config(resolved)

    from pathlib import Path
    from spida.S.segmentation.backends.cellpose_config import resolve_model_dir
    assert resolve_model_dir() == Path.home() / ".cellpose" / "models"


def test_env_var_overlays_when_key_absent_from_json(tmp_path, monkeypatch):
    """Now that the key is in DEFAULT_CONFIG, `_config_from_env` picks it up, so
    a .env value works too."""
    monkeypatch.setenv("CELLPOSE_MODEL_PATH", "/from/env")
    cfg = tmp_path / "BR_config.json"
    cfg.write_text(json.dumps({"ZARR_STORAGE_PATH": "/scratch"}))
    resolved = resolve_config(config_path=str(cfg), defaults=DEFAULT_CONFIG)
    assert resolved["CELLPOSE_MODEL_PATH"] == "/from/env"


def test_json_beats_env(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLPOSE_MODEL_PATH", "/from/env")
    cfg = tmp_path / "BR_config.json"
    cfg.write_text(json.dumps({"CELLPOSE_MODEL_PATH": "/from/json"}))
    resolved = resolve_config(config_path=str(cfg), defaults=DEFAULT_CONFIG)
    assert resolved["CELLPOSE_MODEL_PATH"] == "/from/json"


def test_model_dir_argument_beats_everything(tmp_path, monkeypatch):
    """`--model_dir=` on segment-region is the per-run override."""
    monkeypatch.setenv("CELLPOSE_MODEL_PATH", "/from/config")
    from spida.S.segmentation.backends.cellpose_config import resolve_model_dir
    assert str(resolve_model_dir("/from/cli")) == "/from/cli"


# ---------------------------------------------------------------------------
# the real config files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["config", "CTX", "CB", "HIPP", "THM", "BG"])
def test_real_configs_point_at_an_existing_model_dir(name):
    """Guards the value that was wrong: the configured directory must exist and
    contain cpsam, otherwise cellpose silently substitutes a default model."""
    from pathlib import Path
    path = Path(f"/home/ubuntu/aklein/spida_config/{name}_config.json"
                if name != "config" else
                "/home/ubuntu/aklein/spida_config/config.json")
    if not path.is_file():
        pytest.skip(f"{path} not present on this machine")
    value = json.loads(path.read_text()).get("CELLPOSE_MODEL_PATH")
    assert value, f"{path.name} has no CELLPOSE_MODEL_PATH"
    model_dir = Path(value).expanduser()
    assert model_dir.is_dir(), f"{path.name}: {model_dir} does not exist"
    assert (model_dir / "cpsam").is_file(), f"{path.name}: cpsam missing in {model_dir}"
