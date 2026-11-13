"""Project level: combine-datasets (project JSON -> project zarr + combined h5ad) and sync-project."""
import json
import logging
import warnings

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")
sd = pytest.importorskip("spatialdata")

import anndata as ad
import dask.dataframe as dd
import geopandas as gpd
from click.testing import CliRunner
from shapely.geometry import Point
from spatialdata.models import PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Scale, get_transformation

from spida.config import SpidaConfig
from spida.io.layout import element_keys
from spida.P.project_level import combine_project, load_project, project_paths, read_combined, sync_project
from spida._constants import POINTS_KEY, SHAPES_KEY

N = 12
SAMPLES = [("exp1", "region_A", "R1", "D1"), ("exp1", "region_B", "R1", "D2"), ("exp2", "region_A", "R2", "D1")]


def _region_store(root, exp, reg, offset):
    rng = np.random.default_rng(offset)
    ns = f"{exp}_{reg}"
    k = element_keys("cellpose_cell", ns)
    ids = np.arange(N)
    shapes = gpd.GeoDataFrame(geometry=[Point(offset + i, i).buffer(0.4) for i in ids], index=ids)
    pts = PointsModel.parse(dd.from_pandas(pd.DataFrame(
        {"x": rng.uniform(0, N, 50), "y": rng.uniform(0, N, 50), "gene": pd.Categorical(["g0", "g1"] * 25)}),
        npartitions=1), transformations={"global": Scale([2.0, 2.0], axes=("x", "y"))})
    obs = pd.DataFrame({"EntityID": ids, "cells_region": pd.Categorical([k[SHAPES_KEY]] * N),
                        "experiment": exp, "region": reg, "donor": "old",
                        "base_leiden": pd.Categorical(["0", "1"] * (N // 2))},
                       index=[f"{i}.{exp}.{reg}.cellpose_cell" for i in ids])
    t = ad.AnnData(rng.poisson(3, (N, 4)).astype(np.float32), obs=obs,
                   var=pd.DataFrame({"kind": ["gene"] * 4}, index=[f"g{i}" for i in range(4)]))
    t.layers["raw"] = t.X.copy()
    t.obsm["spatial"] = rng.uniform(size=(N, 2))
    t.obsm["base_pca"] = rng.normal(size=(N, 3))
    t.obsm["blank"] = pd.DataFrame({"Blank-1": rng.poisson(1, N)}, index=t.obs_names)
    t.varm["base_pca"] = rng.normal(size=(4, 3))
    t.uns["cutoffs"] = {"nCount_RNA": offset}
    table = TableModel.parse(t, region=k[SHAPES_KEY], region_key="cells_region", instance_key="EntityID")
    s = sd.SpatialData(shapes={k[SHAPES_KEY]: ShapesModel.parse(shapes)}, points={k[POINTS_KEY]: pts},
                       tables={"cellpose_cell_table_filt": table})
    s.attrs["spida"] = {"namespace": ns, "exp": exp, "reg": reg}
    s.write(root / exp / reg)


@pytest.fixture
def cfg(tmp_path):
    for i, (exp, reg, _, _) in enumerate(SAMPLES):
        _region_store(tmp_path / "zarr", exp, reg, 100 * i)
    (tmp_path / "ad").mkdir()
    return SpidaConfig(zarr_store=tmp_path / "zarr", anndata_store=tmp_path / "ad")


def _raw(**extra):
    return {"project_name": "P1", "prefix": "cellpose_cell", "suffix": "_filt",
            "sample_id": "{brain_region}_{donor}", "required_obs": ["brain_region", "donor"],
            "obs_defaults": {"replicate": "salk", "major_structure": "X"},
            "samples": [{"exp": e, "reg": r, "obs": {"brain_region": b, "donor": d}} for e, r, b, d in SAMPLES],
            **extra}


@pytest.mark.parametrize("edit, match", [
    (lambda r: r.update(obs_default={}), "unknown keys"),
    (lambda r: r["samples"].append(dict(r["samples"][0])), "listed twice"),
    (lambda r: r["samples"][1]["obs"].pop("donor"), "required_obs not set"),
    (lambda r: r.update(sample_id="{brain_region}_{slab}"), "slab"),
    (lambda r: r.update(sample_id="{brain_region}"), "also used by"),
    (lambda r: r.update(project_name="a/b"), "project_name"),
])
def test_load_project_checks(edit, match):
    raw = _raw()
    edit(raw)
    with pytest.raises(ValueError, match=match):
        load_project(raw)


def test_load_project_merge(tmp_path):
    raw = _raw()
    raw["samples"][0]["obs"]["replicate"] = "ucsd"
    (tmp_path / "p.json").write_text(json.dumps(raw))
    spec = load_project(tmp_path / "p.json")
    assert [s.id for s in spec.samples] == ["R1_D1", "R1_D2", "R2_D1"]
    assert spec.samples[0].obs["replicate"] == "ucsd" and spec.samples[1].obs["replicate"] == "salk"
    assert load_project({k: v for k, v in raw.items() if k != "sample_id"}).samples[0].id == "exp1_region_A"


def test_combine(cfg):
    spec = load_project(_raw())
    paths = combine_project(spec, cfg)
    assert paths.zarr == cfg.zarr_store / "_projects" / "P1"
    assert paths.h5ad == cfg.anndata_store / "_projects" / "P1" / "cellpose_cell_table_filt.h5ad"

    p = sd.read_zarr(paths.zarr)
    assert not p.images
    assert set(p.shapes) == {f"cellpose_cell_{s.id}_polygons" for s in spec.samples}
    assert set(p.points) == {f"cellpose_cell_{s.id}_transcripts" for s in spec.samples}
    assert set(p.tables) == {f"cellpose_cell_{s.id}_table_filt" for s in spec.samples}
    assert p.attrs["spida"]["project"] == spec.raw
    t = p["cellpose_cell_R1_D2_table_filt"]
    assert t.uns["spatialdata_attrs"]["region"] == "cellpose_cell_R1_D2_polygons"
    assert set(t.obs["cells_region"]) == {"cellpose_cell_R1_D2_polygons"}
    assert set(t.obs["donor"]) == {"D2"} and set(t.obs["sample_id"]) == {"R1_D2"}
    assert set(t.obs["major_structure"]) == {"X"} and "base_pca" in t.varm  # per-sample table kept whole
    pts = p["cellpose_cell_R1_D2_transcripts"]
    assert set(get_transformation(pts, get_all=True)) == {"R1_D2"}
    assert get_transformation(pts, "R1_D2").scale.tolist() == [2.0, 2.0]

    c = read_combined(paths.h5ad)
    assert c.n_obs == N * len(SAMPLES) and c.n_vars == 4 and "kind" in c.var
    assert set(c.obsm) == {"spatial", "base_pca", "blank"} and not c.varm and not c.obsp
    assert set(c.layers) == {"raw"} and "base_leiden" in c.obs
    assert c.uns["samples"]["R2_D1"]["cutoffs"]["nCount_RNA"] == 200
    assert json.loads(c.uns["spida_project"]) == spec.raw
    assert sorted(c.uns["spatialdata_attrs"]["region"]) == sorted(p.shapes)
    assert np.array_equal(c[t.obs_names].obsm["base_pca"], t.obsm["base_pca"])

    # the combined table maps back onto the project's shapes
    joined = sd.SpatialData(shapes=dict(p.shapes), tables={"combined": c})
    v = sd.get_values("donor", sdata=joined, element_name="cellpose_cell_R2_D1_polygons", table_name="combined")
    assert len(v) == N and set(v["donor"]) == {"D1"}


def test_missing_input_raises_before_write(cfg):
    raw = _raw()
    raw["samples"].append({"exp": "exp9", "reg": "region_Z", "obs": {"brain_region": "R9", "donor": "D9"}})
    with pytest.raises(FileNotFoundError, match="R9_D9"):
        combine_project(load_project(raw), cfg)
    assert not (cfg.zarr_store / "_projects").exists() and not (cfg.anndata_store / "_projects").exists()


def test_refuses_existing_then_overwrite(cfg):
    spec = load_project(_raw())
    paths = combine_project(spec, cfg)
    t0 = paths.h5ad.stat().st_mtime_ns
    with pytest.raises(FileExistsError, match="P1"):
        combine_project(spec, cfg)
    assert paths.h5ad.stat().st_mtime_ns == t0
    combine_project(spec, cfg, overwrite=True)
    assert len(sd.read_zarr(paths.zarr).tables) == len(SAMPLES)
    assert not list(paths.zarr.parent.glob("*__tmp")) and not list(paths.h5ad.parent.glob("*__tmp*"))


def test_unsynced_warning(cfg, caplog):
    copy = ad.read_zarr(cfg.zarr_store / "exp1" / "region_B" / "tables" / "cellpose_cell_table_filt")
    copy.obs["mmc_REF_subclass"] = "Astro"
    out = cfg.anndata_store / "exp1" / "region_B" / "cellpose_cell_table_filt.h5ad"
    out.parent.mkdir(parents=True)
    copy.write_h5ad(out)
    with caplog.at_level(logging.WARNING):
        combine_project(load_project(_raw()), cfg)
    assert "R1_D2" in caplog.text and "mmc_REF_subclass" in caplog.text


def test_sync_project(cfg):
    spec = load_project(_raw())
    paths = combine_project(spec, cfg)
    c = read_combined(paths.h5ad)
    c.obs["leiden"] = pd.Categorical(np.arange(c.n_obs) % 3).astype(str)
    c = c[c.obs_names[2:]].copy()  # e.g. cells dropped at project level
    c.write_h5ad(paths.h5ad)

    with pytest.raises(ValueError, match="--columns"):
        sync_project(spec, cfg, [])
    with pytest.raises(KeyError, match="nope"):
        sync_project(spec, cfg, ["nope"])
    assert sync_project(spec, cfg, ["leiden"]) == ["leiden"]
    p = sd.read_zarr(paths.zarr)
    for s in spec.samples:
        t = p[spec.table(s)]
        expect = c.obs["leiden"].reindex(t.obs_names)
        assert t.obs["leiden"].astype(object).equals(expect.astype(object))
    assert p[spec.table(spec.samples[0])].obs["leiden"].isna().sum() == 2
    with pytest.raises(FileExistsError, match="leiden"):
        sync_project(spec, cfg, ["leiden"])
    sync_project(spec, cfg, ["leiden"], overwrite=True)


def test_cli(cfg, tmp_path):
    from spida.P.cli import cli

    (tmp_path / "p.json").write_text(json.dumps(_raw()))
    args = ["--zarr_store", str(cfg.zarr_store), "--anndata_store", str(cfg.anndata_store)]
    r = CliRunner().invoke(cli, ["combine-datasets", str(tmp_path / "p.json"), *args])
    assert r.exit_code == 0, r.output
    assert project_paths(cfg, load_project(_raw())).h5ad.exists()
    r = CliRunner().invoke(cli, ["combine-datasets", str(tmp_path / "p.json"), *args])
    assert r.exit_code != 0 and isinstance(r.exception, FileExistsError)
    r = CliRunner().invoke(cli, ["sync-project", str(tmp_path / "p.json"), "--columns", "donor", *args])
    assert r.exit_code != 0 and isinstance(r.exception, FileExistsError)  # donor is already on the tables
