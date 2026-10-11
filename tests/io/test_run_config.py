"""Config reaches commands only through run_config (no click default_map); every command runs with no config.

Real click parsing (group make_context included) and real files; nothing mocked.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from spida.P import cli as p_cli

EXP, REG = "E", "region_1"
KEYS = ("ZARR_STORAGE_PATH", "ANNDATA_STORE_PATH")


def _store(root: Path) -> Path:
    zarr = root / EXP / REG
    (zarr / "tables").mkdir(parents=True)
    ad.AnnData(np.ones((3, 2), dtype=np.float32), obs=pd.DataFrame(index=["a", "b", "c"])).write_zarr(
        zarr / "tables" / "p_table")
    return root


def _cfg(path: Path, **kv) -> Path:
    path.write_text(json.dumps({k: str(v) for k, v in kv.items()}))
    return path


def _out(root: Path) -> Path:
    return root / EXP / REG / "p_table.h5ad"


def _invoke(*args, name="write-adata"):
    return CliRunner().invoke(p_cli.cli, [*args[:2], name, "--exp", EXP, "--reg", REG, "--prefix_names", "p", *args[2:]])


@pytest.fixture
def no_env(monkeypatch):
    for k in KEYS:
        monkeypatch.delenv(k, raising=False)


def test_no_command_gets_a_default_map(tmp_path):
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=tmp_path / "z")
    ctx = p_cli.cli.make_context("cli", ["--config", str(cfg), "write-adata", "--exp", EXP, "--prefix_names", "p"])
    assert not ctx.default_map


def test_config_file_supplies_roots(tmp_path, no_env):
    z = _store(tmp_path / "z")
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=z, ANNDATA_STORE_PATH=tmp_path / "ad")
    res = _invoke("--config", str(cfg))
    assert res.exit_code == 0, res.output
    assert ad.read_h5ad(_out(tmp_path / "ad")).n_obs == 3


def test_cli_option_beats_config_file(tmp_path, no_env):
    z = _store(tmp_path / "z")
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=tmp_path / "wrong", ANNDATA_STORE_PATH=tmp_path / "ad")
    res = _invoke("--config", str(cfg), "--zarr_store", str(z), "--anndata_store", str(tmp_path / "mine"))
    assert res.exit_code == 0, res.output
    assert _out(tmp_path / "mine").exists() and not _out(tmp_path / "ad").exists()


def test_env_without_config_file(tmp_path, monkeypatch):
    z = _store(tmp_path / "z")
    monkeypatch.setenv("ZARR_STORAGE_PATH", str(z))
    monkeypatch.setenv("ANNDATA_STORE_PATH", str(tmp_path / "env_ad"))
    res = CliRunner().invoke(p_cli.cli, ["write-adata", "--exp", EXP, "--reg", REG, "--prefix_names", "p"])
    assert res.exit_code == 0, res.output
    assert _out(tmp_path / "env_ad").exists()


def test_config_file_beats_env(tmp_path, monkeypatch):
    z = _store(tmp_path / "z")
    monkeypatch.setenv("ANNDATA_STORE_PATH", str(tmp_path / "env_ad"))
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=z, ANNDATA_STORE_PATH=tmp_path / "file_ad")
    res = _invoke("--config", str(cfg))
    assert res.exit_code == 0, res.output
    assert _out(tmp_path / "file_ad").exists() and not _out(tmp_path / "env_ad").exists()


def test_alias_behaves_the_same(tmp_path, no_env):
    z = _store(tmp_path / "z")
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=z, ANNDATA_STORE_PATH=tmp_path / "ad")
    res = _invoke("--config", str(cfg), name="write_adata")
    assert res.exit_code == 0, res.output
    assert _out(tmp_path / "ad").exists()


def test_unset_root_fails_clearly_not_placeholder(tmp_path, no_env, monkeypatch):
    """Unset ANNDATA_STORE_PATH resolves to None, never to the '/path/to/…' placeholder default."""
    monkeypatch.chdir(tmp_path)
    z = _store(tmp_path / "z")
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=z)
    res = _invoke("--config", str(cfg))
    assert res.exit_code != 0
    assert "anndata" in str(res.exception).lower()
    assert not Path("/path/to/anndata/store").exists()


def test_real_subprocess(tmp_path):
    """Outside CliRunner: the module entry point, a different cwd, config file only."""
    z = _store(tmp_path / "z")
    cfg = _cfg(tmp_path / "c.json", ZARR_STORAGE_PATH=z, ANNDATA_STORE_PATH=tmp_path / "ad")
    env = {k: v for k, v in os.environ.items() if k not in KEYS}
    r = subprocess.run([sys.executable, "-m", "spida.P.cli", "--config", str(cfg), "write-adata", "--exp", EXP,
                        "--reg", REG, "--prefix_names", "p"], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    assert _out(tmp_path / "ad").exists()


CLIS = ["spida.S.cli", "spida.P.cli", "spida.I.cli", "spida.site"]


@pytest.mark.parametrize("mod", CLIS)
def test_every_help_runs_without_config_or_env(mod, tmp_path, monkeypatch):
    """No --config, no .env in the cwd, no SPIDA env vars: --help works for the group and every command."""
    import importlib
    from spida.config import DEFAULT_CONFIG

    monkeypatch.chdir(tmp_path)
    for k in DEFAULT_CONFIG:
        monkeypatch.delenv(k, raising=False)
    cli = importlib.import_module(mod).cli
    assert CliRunner().invoke(cli, ["--help"]).exit_code == 0
    for name in cli.commands:
        res = CliRunner().invoke(cli, [name, "--help"])
        assert res.exit_code == 0, f"{mod} {name}: {res.output[-400:]}"


def test_real_command_with_cli_paths_only(tmp_path, monkeypatch):
    """No config file and no env: the command runs when the CLI gives what it needs."""
    from spida.config import DEFAULT_CONFIG

    monkeypatch.chdir(tmp_path)
    for k in DEFAULT_CONFIG:
        monkeypatch.delenv(k, raising=False)
    z = _store(tmp_path / "z")
    res = CliRunner().invoke(p_cli.cli, ["write-adata", "--exp", EXP, "--reg", REG, "--prefix_names", "p",
                                         "--zarr_store", str(z), "--anndata_store", str(tmp_path / "ad")])
    assert res.exit_code == 0, res.output
    assert _out(tmp_path / "ad").exists()


def test_env_config_file_is_rejected_clearly(tmp_path):
    env = tmp_path / "old.env"
    env.write_text("ZARR_STORAGE_PATH=/x\n")
    res = CliRunner().invoke(p_cli.cli, ["--config", str(env), "write-adata", "--exp", EXP, "--prefix_names", "p"])
    assert res.exit_code != 0 and "JSON" in str(res.exception)

