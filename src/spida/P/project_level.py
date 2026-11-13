"""
Combine many individual zarr objects into a single zarr object describing a single project. 
Omits the images to save space and compute time. 

Project JSON describes a project
    {
      "project_name" - required; output folder name
      "prefix" - required; segmentation
      "suffix" - table suffix (default "")
      "sample_id" - element-name template over exp, reg and obs (default "{exp}_{reg}")
      "required_obs" - obs every sample must set (default none)
      "obs_defaults" - obs for every sample
      "samples": [{"exp": ..., "reg": ..., "obs": {...}}, ...] - a sample's obs overrides obs_defaults
    }
"""

import json
import logging
import os
import re
import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from spida.io.layout import TreeLayout, table_key
from spida.io.store import is_writeable, obs_columns, read_table, require, sync_obs
from spida._constants import POINTS_KEY, SHAPES_KEY, TABLE_KEY

logger = logging.getLogger(__package__)

PROJECTS_DIR = "_projects"  # under each store root
TMP_SUFFIX = "__tmp"
_KEYS = {"project_name", "prefix", "suffix", "sample_id", "required_obs", "obs_defaults", "samples"}
_SAMPLE_KEYS = {"exp", "reg", "obs"}
_NAME = re.compile(r"[A-Za-z0-9_.-]+")  # allowed in element and folder names


@dataclass(frozen=True)
class Sample:
    exp: str
    reg: str
    id: str  # short name used in the project's element names and coordinate system
    obs: dict  # obs_defaults merged with the sample's obs


@dataclass(frozen=True)
class ProjectSpec:
    name: str
    prefix: str
    suffix: str
    samples: tuple
    raw: dict  # the JSON as read

    @property
    def table_key(self) -> str:
        """The source table in each sample's store, e.g. cellpose_cell_table_filt."""
        return table_key(self.prefix, self.suffix)

    def shapes(self, sample: Sample) -> str:
        return f"{self.prefix}_{sample.id}_polygons"

    def points(self, sample: Sample) -> str:
        return f"{self.prefix}_{sample.id}_transcripts"

    def table(self, sample: Sample) -> str:
        return f"{self.prefix}_{sample.id}_table{self.suffix}"


@dataclass(frozen=True)
class ProjectPaths:
    zarr: Path  # zarr_store/_projects/<name>
    h5ad: Path  # anndata_store/_projects/<name>/<prefix>_table<suffix>.h5ad
    figures: Path | None  # image_store/_projects/<name>


def load_project(
    project,
) -> ProjectSpec:
    """
    Read and check a project JSON before any work.
    Parameters:
    project (str | Path | dict): The project JSON file, or its contents.
    """
    raw = project if isinstance(project, dict) else json.loads(Path(project).read_text())
    errors = []
    if unknown := set(raw) - _KEYS:
        errors.append(f"unknown keys {sorted(unknown)} (allowed: {sorted(_KEYS)})")
    errors += [f"missing {k!r}" for k in ("project_name", "prefix", "samples") if not raw.get(k)]
    if errors:
        raise ValueError("project JSON: " + "; ".join(errors))

    name, template = raw["project_name"], raw.get("sample_id", "{exp}_{reg}")
    if not _NAME.fullmatch(name):
        errors.append(f"project_name {name!r} may only use letters, digits, '_', '.', '-'")
    defaults, required = raw.get("obs_defaults", {}), raw.get("required_obs", [])
    samples, seen_pairs, seen_ids = [], set(), {}
    for i, s in enumerate(raw["samples"]):
        where = f"sample {i} ({s.get('exp')} / {s.get('reg')})"
        if unknown := set(s) - _SAMPLE_KEYS:
            errors.append(f"{where}: unknown keys {sorted(unknown)}")
        if not s.get("exp") or not s.get("reg"):
            errors.append(f"{where}: needs exp and reg")
            continue
        if (s["exp"], s["reg"]) in seen_pairs:
            errors.append(f"{where}: listed twice")
        seen_pairs.add((s["exp"], s["reg"]))
        obs = {**defaults, **s.get("obs", {})}
        if missing := [c for c in required if c not in obs]:
            errors.append(f"{where}: required_obs not set: {missing}")
        try:
            sid = template.format(**{"exp": s["exp"], "reg": s["reg"], **obs})
        except KeyError as e:
            errors.append(f"{where}: sample_id {template!r} uses {e}, which this sample does not set")
            continue
        if not _NAME.fullmatch(sid):
            errors.append(f"{where}: sample_id {sid!r} may only use letters, digits, '_', '.', '-'")
        if sid in seen_ids:
            errors.append(f"{where}: sample_id {sid!r} also used by {seen_ids[sid]}")
        seen_ids[sid] = where
        samples.append(Sample(s["exp"], s["reg"], sid, obs))
    if errors:
        raise ValueError("project JSON:\n  " + "\n  ".join(errors))
    return ProjectSpec(name, raw["prefix"], raw.get("suffix", ""), tuple(samples), raw)


def project_paths(
    config,
    spec: ProjectSpec,
) -> ProjectPaths:
    """
    Where a project's outputs live.
    Parameters:
    config (SpidaConfig): Store roots (zarr_store and anndata_store required).
    spec (ProjectSpec): The project.
    """
    figures = Path(config.image_store) / PROJECTS_DIR / spec.name if config.image_store else None
    return ProjectPaths(
        zarr=Path(config.require("zarr_store")) / PROJECTS_DIR / spec.name,
        h5ad=Path(config.require("anndata_store")) / PROJECTS_DIR / spec.name / f"{spec.table_key}.h5ad",
        figures=figures,
    )


def _check_inputs(spec: ProjectSpec, layout: TreeLayout) -> dict:
    # every sample's store, shapes, points and table; raises listing all that are missing
    datasets, missing = {}, []
    for s in spec.samples:
        d = layout.dataset(s.exp, s.reg)
        try:
            k = d.keys(spec.prefix)
            require(d.zarr, [k[SHAPES_KEY], k[POINTS_KEY], spec.table_key])
        except (FileNotFoundError, ValueError) as e:
            missing.append(f"{s.id} ({s.exp} / {s.reg}): {e}")
        datasets[s.id] = d
    if missing:
        raise FileNotFoundError("project inputs:\n  " + "\n  ".join(missing))
    return datasets


def _warn_unsynced(spec: ProjectSpec, datasets: dict) -> None:
    # annotation columns on a sample's h5ad copy that its zarr table (the source used here) lacks
    import anndata as ad

    for s in spec.samples:
        path = datasets[s.id].table_h5ad(spec.table_key) if datasets[s.id].anndata else None
        if path is None or not path.exists():
            continue
        a = ad.read_h5ad(path, backed="r")
        cols = [c for c in a.obs.columns if c.startswith(("mmc_", "allcools_"))]
        a.file.close()
        if lost := sorted(set(cols) - set(obs_columns(datasets[s.id].zarr, spec.table_key))):
            logger.warning("%s: %d annotation columns are in %s but not in its zarr table, so not in the "
                           "project (run sync-annotations first): %s", s.id, len(lost), path.name, ", ".join(lost))


def _set_obs(table, sample: Sample) -> None:
    # project obs (strings as categoricals) + sample_id
    for col, value in {**sample.obs, "sample_id": sample.id}.items():
        table.obs[col] = pd.Categorical([value] * table.n_obs) if isinstance(value, str) else value


def _annotate(table, region, region_key: str, instance_key: str):
    # re-point a table at its shapes element(s)
    from spatialdata.models import TableModel

    table.uns.pop("spatialdata_attrs", None)
    return TableModel.parse(table, region=region, region_key=region_key, instance_key=instance_key)


def _to_sample_cs(element, sample_id: str) -> None:
    # the sample's "global" becomes its own coordinate system; others (e.g. pixel) are dropped
    from spatialdata.transformations import get_transformation, set_transformation

    transforms = get_transformation(element, get_all=True)
    if "global" not in transforms:
        raise ValueError(f"{sample_id}: element has no 'global' coordinate system ({list(transforms)})")
    set_transformation(element, {sample_id: transforms["global"]}, set_all=True)


def _stack_obsm(tables: list) -> dict:
    # obsm keys every sample has with the same width; the rest are left out with a warning
    import scipy.sparse as sp

    keys = [set(t.obsm) for t in tables]
    common = set.intersection(*keys)
    if partial := sorted(set.union(*keys) - common):
        logger.warning("obsm %s not in every sample; left out of the combined table", partial)
    out = {}
    for k in sorted(common):
        parts = [t.obsm[k] for t in tables]
        if all(isinstance(p, pd.DataFrame) for p in parts) and len({tuple(p.columns) for p in parts}) == 1:
            out[k] = pd.concat(parts)
        elif not any(isinstance(p, pd.DataFrame) for p in parts) and len({p.shape[1:] for p in parts}) == 1:
            out[k] = sp.vstack(parts).tocsr() if any(sp.issparse(p) for p in parts) else np.vstack(parts)
        else:
            logger.warning("obsm[%s] differs in shape across samples; left out of the combined table", k)
    return out


def _combine_tables(spec: ProjectSpec, tables: dict, region_key: str, instance_key: str):
    # concatenate the per-sample tables: X, layers, obs; obsm stacked; per-sample uns under uns["samples"]
    import anndata as ad

    tables = {sid: t.copy() for sid, t in tables.items()}
    ids = list(tables)
    var_sets = {sid: set(t.var_names) for sid, t in tables.items()}
    if len({frozenset(v) for v in var_sets.values()}) > 1:
        shared = set.intersection(*var_sets.values())
        logger.warning("samples have different genes; the combined table keeps the %d shared (%s)", len(shared),
                       ", ".join(f"{sid}: {len(v)}" for sid, v in var_sets.items()))
    obsm = _stack_obsm(list(tables.values()))
    uns = {"samples": {sid: {k: v for k, v in t.uns.items() if k != "spatialdata_attrs"}
                       for sid, t in tables.items()}}
    for t in tables.values():  # per-sample bases and graphs do not concatenate
        t.obsm.clear(); t.varm.clear(); t.obsp.clear(); t.uns.clear()

    adata = ad.concat([tables[sid] for sid in ids], join="inner", merge="same", index_unique=None)
    if adata.obs_names.has_duplicates:
        raise ValueError(f"{adata.obs_names.duplicated().sum()} cell names occur in more than one sample")
    for k, v in obsm.items():
        adata.obsm[k] = v.loc[adata.obs_names] if isinstance(v, pd.DataFrame) else v
    adata.uns.update(uns)
    adata.uns["spida_project"] = json.dumps(spec.raw)
    return _annotate(adata, [spec.shapes(s) for s in spec.samples], region_key, instance_key)


def read_combined(
    path,
    backed=None,
):
    """
    The combined h5ad, ready to pair with the project zarr's shapes (h5ad stores the region list as an array).
    Parameters:
    path (str | Path): The combined h5ad.
    backed (str, optional): Passed to anndata.read_h5ad.
    """
    import anndata as ad

    adata = ad.read_h5ad(path, backed=backed)
    attrs = adata.uns.get("spatialdata_attrs", {})
    if isinstance(attrs.get("region"), np.ndarray):
        attrs["region"] = attrs["region"].tolist()
    return adata


def combine_project(
    spec: ProjectSpec,
    config,
    overwrite: bool = False,
) -> ProjectPaths:
    """
    Build a project: the project zarr (per-sample shapes, points and tables; no images) and the combined h5ad.
    Inputs are each sample's zarr store; every input and output is checked before any data is read.
    Parameters:
    spec (ProjectSpec): The project (load_project).
    config (SpidaConfig): Store roots.
    overwrite (bool): Replace an existing project zarr / combined h5ad.
    """
    from spida.io.store import _sd

    sd = _sd()
    paths = project_paths(config, spec)
    datasets = _check_inputs(spec, TreeLayout(config))
    is_writeable(files=[paths.zarr, paths.h5ad], overwrite=overwrite)
    _warn_unsynced(spec, datasets)

    shapes, points, tables = {}, {}, {}
    region_key = instance_key = None
    for s in spec.samples:
        d = datasets[s.id]
        logger.info("Adding %s (%s / %s)", s.id, s.exp, s.reg)
        src = d.keys(spec.prefix)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            sdata = sd.read_zarr(d.zarr, selection=("shapes", "points"))
        shp, pts = sdata[src[SHAPES_KEY]], sdata[src[POINTS_KEY]]
        for element in (shp, pts):
            _to_sample_cs(element, s.id)
        table = read_table(d.zarr, spec.table_key)
        attrs = table.uns["spatialdata_attrs"]
        if region_key is None:
            region_key, instance_key = attrs["region_key"], attrs["instance_key"]
        elif (attrs["region_key"], attrs["instance_key"]) != (region_key, instance_key):
            raise ValueError(f"{s.id}: table links shapes by {attrs['region_key']}/{attrs['instance_key']}, "
                             f"other samples by {region_key}/{instance_key}")
        _set_obs(table, s)
        table.obs[region_key] = pd.Categorical([spec.shapes(s)] * table.n_obs)
        shapes[spec.shapes(s)], points[spec.points(s)] = shp, pts
        tables[s.id] = _annotate(table, spec.shapes(s), region_key, instance_key)

    combined = _combine_tables(spec, tables, region_key, instance_key)

    # project zarr: written beside, then swapped in
    tmp = paths.zarr.with_name(paths.zarr.name + TMP_SUFFIX)
    if tmp.exists():
        logger.warning("Removing %s left by an earlier failed write.", tmp)
        shutil.rmtree(tmp)
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        project = sd.SpatialData(shapes=shapes, points=points,
                                 tables={spec.table(s): tables[s.id] for s in spec.samples})
        project.attrs["spida"] = {"project": spec.raw}
        project.write(tmp)
    if paths.zarr.exists():
        shutil.rmtree(paths.zarr)
    tmp.rename(paths.zarr)
    logger.info("Wrote project store %s", paths.zarr)

    # combined h5ad: written beside, then swapped in
    paths.h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = paths.h5ad.with_name(f"{paths.h5ad.stem}{TMP_SUFFIX}.h5ad")
    combined.write_h5ad(tmp_h5ad)
    os.replace(tmp_h5ad, paths.h5ad)
    logger.info("Wrote combined table %s (%d cells x %d genes)", paths.h5ad, combined.n_obs, combined.n_vars)
    return paths


def sync_project(
    spec: ProjectSpec,
    config,
    columns,
    overwrite: bool = False,
) -> list[str]:
    """
    Copy obs ``columns`` from the combined h5ad into each sample's table in the project zarr (obs only).
    Cells are split by the region key (cells_region). Every table is checked before any is written.
    Parameters:
    spec (ProjectSpec): The project (load_project).
    config (SpidaConfig): Store roots.
    columns (Iterable[str]): obs columns of the combined h5ad to copy.
    overwrite (bool): Replace columns that already exist in the project tables.
    """
    import anndata as ad

    cols = list(columns or [])
    if not cols:
        raise ValueError("sync-project: name the columns to copy (--columns)")
    paths = project_paths(config, spec)
    require(paths.zarr, [spec.table(s) for s in spec.samples], files=[paths.h5ad])
    a = ad.read_h5ad(paths.h5ad, backed="r")
    obs, region_key = a.obs.copy(), a.uns["spatialdata_attrs"]["region_key"]
    a.file.close()
    if missing := [c for c in cols if c not in obs.columns]:
        raise KeyError(f"{paths.h5ad.name} has no columns {missing}")
    is_writeable(paths.zarr, columns={spec.table(s): cols for s in spec.samples}, overwrite=overwrite)

    for s in spec.samples:
        sub = obs[obs[region_key] == spec.shapes(s)]
        if sub.empty:
            logger.warning("%s: no cells in %s; its table gets missing values", s.id, paths.h5ad.name)
        sync_obs(paths.zarr, spec.table(s), sub, cols, source=f"{paths.h5ad.name} [{s.id}]")
    return cols
