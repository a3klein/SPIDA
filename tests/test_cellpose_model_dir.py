"""Where cellpose weights come from and where downloads land.

One directory per run: ``--model_dir`` > ``CELLPOSE_MODEL_PATH`` (SPIDA config, filled in by
``segment-region``) > ``CELLPOSE_LOCAL_MODELS_PATH`` > ``~/.cellpose/models``. cellpose fixes its
cache (``models.MODEL_DIR``) from ``CELLPOSE_LOCAL_MODELS_PATH`` when ``cellpose.models`` is first
imported, so ``_load_model`` points cellpose at the resolved directory before that import.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from spida.config import DEFAULT_CONFIG, load
from spida.S.segmentation.backends.cellpose_config import (
    CELLPOSE_ENV_VAR,
    CellposeConfig,
    point_cellpose_at,
)

SPIDA_ROOT = Path(__file__).resolve().parents[1]
CELLPOSE_PY = SPIDA_ROOT / ".pixi" / "envs" / "cellpose" / "bin" / "python"


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for var in ("CELLPOSE_MODEL_PATH", CELLPOSE_ENV_VAR):
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------

def test_config_key_is_unset_by_default():
    assert DEFAULT_CONFIG["CELLPOSE_MODEL_PATH"] == ""
    assert load().cellpose_model_path is None


def test_config_key_from_json_and_env(tmp_path, monkeypatch):
    cfg = tmp_path / "c.json"
    cfg.write_text(json.dumps({"CELLPOSE_MODEL_PATH": "~/weights"}))
    assert load(cfg).cellpose_model_path == Path.home() / "weights"
    monkeypatch.setenv("CELLPOSE_MODEL_PATH", "/from/env")
    assert load().cellpose_model_path == Path("/from/env")


def test_cli_group_does_not_touch_the_environment(tmp_path):
    """The old group exported CELLPOSE_* on every call; nothing should now."""
    from spida.S import cli as s_cli

    cfg = tmp_path / "c.json"
    cfg.write_text(json.dumps({"CELLPOSE_MODEL_PATH": "/models"}))
    before = dict(os.environ)
    s_cli.cli.make_context("cli", ["--config", str(cfg), "segment-region", "cellpose", "--list-params"])
    assert dict(os.environ) == before


# ---------------------------------------------------------------------------
# point_cellpose_at
# ---------------------------------------------------------------------------

def test_point_sets_expanded_env_before_import(monkeypatch):
    monkeypatch.delitem(sys.modules, "cellpose.models", raising=False)
    point_cellpose_at("~/cp_weights")
    assert os.environ[CELLPOSE_ENV_VAR] == str(Path.home() / "cp_weights")
    assert "cellpose.models" not in sys.modules          # does not import cellpose itself


def test_point_updates_an_already_imported_cellpose(tmp_path, monkeypatch):
    import types
    fake = types.ModuleType("cellpose.models")
    fake.MODEL_DIR = Path("/old")
    monkeypatch.setitem(sys.modules, "cellpose.models", fake)
    point_cellpose_at(tmp_path)
    assert fake.MODEL_DIR == tmp_path


# ---------------------------------------------------------------------------
# _load_model import order, against a fake cellpose that behaves like the real one
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_cellpose(tmp_path, monkeypatch):
    """A `cellpose` package whose `models` captures the env var at import, as cellpose 4.x does."""
    pkg = tmp_path / "fakepkgs" / "cellpose"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "core.py").write_text("def use_gpu():\n    return True\n")
    (pkg / "models.py").write_text(textwrap.dedent('''
        import os
        from pathlib import Path
        _ENV = os.environ.get("CELLPOSE_LOCAL_MODELS_PATH")
        MODEL_DIR = Path(_ENV) if _ENV else Path.home().joinpath(".cellpose", "models")
        DIR_AT_IMPORT = MODEL_DIR
        MODEL_NAMES = ["cpsam", "cpsam_v2"]
        class CellposeModel:
            def __init__(self, pretrained_model=None, gpu=False, use_bfloat16=False):
                self.pretrained_model = pretrained_model
    '''))
    for m in ("cellpose", "cellpose.models", "cellpose.core"):
        monkeypatch.delitem(sys.modules, m, raising=False)
    monkeypatch.syspath_prepend(str(pkg.parent))
    yield
    for m in ("cellpose", "cellpose.models", "cellpose.core"):
        sys.modules.pop(m, None)


def test_load_model_points_cellpose_before_import(fake_cellpose, tmp_path):
    from spida.S.segmentation.backends import cellpose as cpb

    weights = tmp_path / "w"
    weights.mkdir()
    (weights / "cpsam").write_bytes(b"x")
    model, res = cpb._load_model(CellposeConfig.from_kwargs(model="cpsam", model_dir=str(weights)))

    import cellpose.models as m
    assert m.DIR_AT_IMPORT == weights                   # env was set before the first import
    assert model.pretrained_model == str(weights / "cpsam") and res.exists_locally


def test_download_lands_where_spida_looks(fake_cellpose, tmp_path):
    """Weights absent: bare name goes to cellpose, whose cache is the resolved dir."""
    from spida.S.segmentation.backends import cellpose as cpb

    model, res = cpb._load_model(CellposeConfig.from_kwargs(model="cpsam_v2", model_dir=str(tmp_path / "dl")))
    import cellpose.models as m
    assert res.will_download and model.pretrained_model == "cpsam_v2"
    assert m.MODEL_DIR == tmp_path / "dl" == res.model_dir


def test_second_run_in_the_same_process_repoints(fake_cellpose, tmp_path):
    from spida.S.segmentation.backends import cellpose as cpb

    cpb._load_model(CellposeConfig.from_kwargs(model="cpsam_v2", model_dir=str(tmp_path / "a")))
    cpb._load_model(CellposeConfig.from_kwargs(model="cpsam_v2", model_dir=str(tmp_path / "b")))
    import cellpose.models as m
    assert m.DIR_AT_IMPORT == tmp_path / "a" and m.MODEL_DIR == tmp_path / "b"


# ---------------------------------------------------------------------------
# segment-region fills model_dir from the config
# ---------------------------------------------------------------------------

def _segment(monkeypatch, tmp_path, method, *extra, cfg_values=None):
    from spida.S import cli as s_cli

    seen = {}
    monkeypatch.setattr("spida.S.segmentation.main.segment_region", lambda *a, **k: seen.update(k))
    args = []
    if cfg_values is not None:
        cfg = tmp_path / "c.json"
        cfg.write_text(json.dumps(cfg_values))
        args = ["--config", str(cfg)]
    res = CliRunner().invoke(s_cli.cli, [*args, "segment-region", method, "--zarr_path", str(tmp_path / "z"),
                                         "--raw_path", str(tmp_path / "r"), "--seg_path", str(tmp_path / "s"), *extra])
    assert res.exit_code == 0, res.output
    return seen


def test_config_dir_reaches_cellpose(monkeypatch, tmp_path):
    k = _segment(monkeypatch, tmp_path, "cellpose", cfg_values={"CELLPOSE_MODEL_PATH": "/shared/weights"})
    assert k["model_dir"] == "/shared/weights"


def test_cli_model_dir_beats_config(monkeypatch, tmp_path):
    k = _segment(monkeypatch, tmp_path, "cellpose", "--model_dir=/cli/weights",
                 cfg_values={"CELLPOSE_MODEL_PATH": "/shared/weights"})
    assert k["model_dir"] == "/cli/weights"


def test_env_config_key_reaches_cellpose(monkeypatch, tmp_path):
    monkeypatch.setenv("CELLPOSE_MODEL_PATH", "/env/weights")
    assert _segment(monkeypatch, tmp_path, "cellpose")["model_dir"] == "/env/weights"


def test_no_config_leaves_resolution_to_cellpose_defaults(monkeypatch, tmp_path):
    assert "model_dir" not in _segment(monkeypatch, tmp_path, "cellpose")


def test_other_backends_get_no_model_dir(monkeypatch, tmp_path):
    k = _segment(monkeypatch, tmp_path, "proseg", cfg_values={"CELLPOSE_MODEL_PATH": "/shared/weights"})
    assert "model_dir" not in k


# ---------------------------------------------------------------------------
# the real cellpose (cellpose env interpreter, read-only; no pixi sync)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not CELLPOSE_PY.exists(), reason="cellpose pixi env not installed")
def test_real_cellpose_uses_the_pointed_dir(tmp_path):
    script = textwrap.dedent(f'''
        from pathlib import Path
        from spida.S.segmentation.backends.cellpose_config import point_cellpose_at
        a, b = Path({str(tmp_path / "a")!r}), Path({str(tmp_path / "b")!r})
        point_cellpose_at(a)
        from cellpose import models
        assert models.MODEL_DIR == a, models.MODEL_DIR
        b.mkdir(); (b / "fake_weights").write_bytes(b"x")
        point_cellpose_at(b)
        assert models.MODEL_DIR == b, models.MODEL_DIR
        assert models.cache_model_path("fake_weights") == str(b / "fake_weights")
        print("ok")
    ''')
    env = {k: v for k, v in os.environ.items() if k not in ("CELLPOSE_MODEL_PATH", CELLPOSE_ENV_VAR)}
    env["PYTHONPATH"] = str(SPIDA_ROOT / "src")
    r = subprocess.run([str(CELLPOSE_PY), "-c", script], env=env, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr[-3000:]
