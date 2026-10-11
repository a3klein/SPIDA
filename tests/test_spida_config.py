import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from spida.config import DEFAULT_CONFIG, ConfigDefaultGroup, SpidaConfig, load


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for k in (*DEFAULT_CONFIG, "ANNOTATIONS_STORE_PATH"):
        monkeypatch.delenv(k, raising=False)


def _json(tmp_path, d):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(d))
    return p


def test_unset_is_none():
    cfg = load()
    assert cfg.zarr_store is None and cfg.mmc_store_path is None
    assert cfg.cellpose_model_path is None  # unset: cellpose resolves its own default


def test_precedence_override_over_file_over_env(tmp_path, monkeypatch):
    monkeypatch.setenv("ZARR_STORAGE_PATH", "/env")
    monkeypatch.setenv("IMAGE_STORE_PATH", "/env_img")
    cfg_path = _json(tmp_path, {"ZARR_STORAGE_PATH": "/file"})
    assert load().zarr_store == Path("/env")
    assert load(cfg_path).zarr_store == Path("/file")
    assert load(cfg_path).image_store == Path("/env_img")
    assert load(cfg_path, zarr_store="/cli").zarr_store == Path("/cli")
    assert load(cfg_path, zarr_store=None).zarr_store == Path("/file")


def test_field_names_follow_cli_options(tmp_path):
    cfg = load(_json(tmp_path, {"PROCESSED_ROOT_PATH": "/raw", "MMC_DIR": "/mmc", "RUST_BIN_PATH": "/rust"}))
    assert cfg.root_path == Path("/raw") and cfg.mmc_store_path == Path("/mmc") and cfg.rust_bin_path == Path("/rust")


def test_placeholders_and_blanks_are_unset(tmp_path):
    cfg = load(_json(tmp_path, {"ZARR_STORAGE_PATH": "/path/to/zarr/storage", "IMAGE_STORE_PATH": ""}))
    assert cfg.zarr_store is None and cfg.image_store is None


def test_annotations_alias(monkeypatch):
    monkeypatch.setenv("ANNOTATIONS_STORE_PATH", "/ann")
    assert load().annotation_store == Path("/ann")


def test_require_names_the_key():
    with pytest.raises(ValueError, match="ZARR_STORAGE_PATH"):
        SpidaConfig().require("zarr_store")


def test_group_resolves_once_for_options_and_ctx(tmp_path):
    @click.group(cls=ConfigDefaultGroup)
    @click.option("--config", "config", default=None)
    @click.pass_context
    def cli(ctx, config):
        pass

    @cli.command()
    @click.option("--zarr_store", "zarr_store", default=None)
    @click.pass_context
    def show(ctx, zarr_store):
        from spida.io.options import run_config
        click.echo(f"{zarr_store} {ctx.obj['config'].zarr_store} {run_config(ctx, zarr_store=zarr_store).zarr_store}")

    cfg = str(_json(tmp_path, {"ZARR_STORAGE_PATH": "/file"}))
    res = CliRunner().invoke(cli, ["--config", cfg, "show"])
    assert res.exit_code == 0, res.output
    assert res.output.strip() == "None /file /file"  # no default_map: the option stays None
    res = CliRunner().invoke(cli, ["--config", cfg, "show", "--zarr_store", "/cli"])
    assert res.output.strip() == "/cli /file /cli"
