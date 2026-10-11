"""Shared CLI options that name a region (--exp/--reg) or point at its folders directly (--*_path)."""

from __future__ import annotations

import dataclasses

import click

from spida.config import SpidaConfig, _as_path
from spida.io.layout import Dataset, dataset_from

_DATASET_OPTIONS = [
    click.option("--exp", "exp", default=None, type=str, help="Experiment name."),
    click.option("--reg", "reg", default=None, type=str, help="Region name."),
]


_OVERRIDE_OPTIONS = [
    click.option("--zarr_path", "zarr_path", default=None, type=click.Path(), help="Override the default zarr store."),
    click.option("--raw_path", "raw_path", default=None, type=click.Path(), help="Override the default raw MERSCOPE region dir."),
    click.option("--seg_path", "seg_path", default=None, type=click.Path(), help="Override the default segmentation dir."),
    click.option("--figures_path", "figures_path", default=None, type=click.Path(), help="Override the default image dir."),
]


_TABLE_OPTIONS = [
    click.option("--anndata_path", "anndata_path", default=None, type=click.Path(), help="Override the default AnnData (.h5ad) files dir."),
    click.option("--annotations_path", "annotations_path", default=None, type=click.Path(), help="Override the default annotation intermediates dir."),
]


def dataset_options(f):
    """--exp / --reg: the experiment and region names."""
    for opt in reversed(_DATASET_OPTIONS):
        f = opt(f)
    return f


def dataset_overrides(f):
    """--zarr_path / --raw_path / --seg_path / --figures_path. Place last, just above @click.pass_context."""
    for opt in reversed(_OVERRIDE_OPTIONS):
        f = opt(f)
    return f


def table_options(f):
    """--anndata_path / --annotations_path, for commands that read or write AnnData files. Place last."""
    for opt in reversed(_TABLE_OPTIONS):
        f = opt(f)
    return f


def run_config(ctx, **roots) -> SpidaConfig:
    """The group's SpidaConfig with any per-command root options (e.g. --zarr_store) applied."""
    cfg = ctx.obj["config"] if ctx.obj and "config" in ctx.obj else SpidaConfig()
    known = {f.name for f in dataclasses.fields(SpidaConfig)}
    return dataclasses.replace(cfg, **{k: _as_path(v) for k, v in roots.items() if k in known and v is not None})


def dataset_from_options(cfg: SpidaConfig, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                         anndata_path=None, annotations_path=None) -> Dataset:
    try:
        return dataset_from(cfg, exp, reg, zarr=zarr_path, raw=raw_path, seg=seg_path, figures=figures_path,
                            anndata=anndata_path, annotations=annotations_path)
    except ValueError as e:
        raise click.UsageError(str(e))
