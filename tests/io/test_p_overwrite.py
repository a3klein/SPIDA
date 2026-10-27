"""Phase-2 guards for P steps that share a table, call-region-tz and write-adata.

The heavy compute (setup, SOLO, RESOLVI, region calling) is replaced by small stand-ins that add the
same outputs; the guards, the store writes and the step-to-step interactions are real.
"""
import os
import json
import warnings

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")
sd = pytest.importorskip("spatialdata")

import anndata as ad
import dask.dataframe as dd
import geopandas as gpd
from shapely.geometry import Point
from spatialdata.models import PointsModel, ShapesModel, TableModel

from spida.io.layout import Dataset, element_keys
from spida._constants import POINTS_KEY, SHAPES_KEY, TABLE_KEY
from spida.P import main as m

K = element_keys("cellpose_cell", None)
TABLE, FILT = K[TABLE_KEY], K[TABLE_KEY] + "_filt"
N = 40


@pytest.fixture
def ds(tmp_path):
    rng = np.random.default_rng(0)
    shapes = gpd.GeoDataFrame({"EntityID": np.arange(N)}, geometry=[Point(i, i).buffer(1) for i in range(N)],
                              index=np.arange(N).astype(str))
    obs = pd.DataFrame({"cells_region": pd.Categorical([K[SHAPES_KEY]] * N), "EntityID": np.arange(N),
                        "pass_qc": np.arange(N) % 4 != 0}, index=np.arange(N).astype(str))
    table = TableModel.parse(ad.AnnData(rng.poisson(3, (N, 5)).astype(np.float32), obs=obs),
                             region=K[SHAPES_KEY], region_key="cells_region", instance_key="EntityID")
    pts = pd.DataFrame({"x": rng.uniform(0, 40, 300), "y": rng.uniform(0, 40, 300), "gene": rng.choice(list("AB"), 300)})
    s = sd.SpatialData(shapes={K[SHAPES_KEY]: ShapesModel.parse(shapes)}, tables={TABLE: table},
                       points={element_keys("default", None)[POINTS_KEY]: PointsModel.parse(dd.from_pandas(pts, npartitions=1))})
    s.attrs["spida"] = {"namespace": "", "exp": None, "reg": None}
    d = Dataset(zarr=tmp_path / "s.zarr", figures=tmp_path / "fig", anndata=tmp_path / "ad")
    s.write(d.zarr)
    return d


@pytest.fixture
def stand_ins(monkeypatch):
    import sys
    import types
    import spida.P.setup_adata as setup_mod
    import spida.utilities.ad_utils as ad_utils

    scvi_mod = types.ModuleType("spida.P.scvi_toolkit")  # stand-in: importing the real one loads scvi-tools / torch (~4 min)
    monkeypatch.setitem(sys.modules, scvi_mod.__name__, scvi_mod)

    def setup(adata, *a, **k):
        out = adata[adata.obs["pass_qc"]].copy()
        out.obs["base_leiden"] = pd.Categorical(["0"] * out.n_obs)
        return out

    def doublets(adata, threshold=0.5, **k):
        adata.obs["doublet_score"] = 0.1
        adata.obs["doublet_bool"] = False
        return adata

    def resolvi(adata, **k):
        for c in ("true_proportion", "diffusion_proportion", "background_proportion"):
            adata.obs[c] = 0.3
        return adata

    def embed(adata, key_added="", **k):
        adata.obs[f"{key_added}leiden"] = pd.Categorical(["0"] * adata.n_obs)

    monkeypatch.setattr(setup_mod, "run_setup", setup)
    scvi_mod.identify_doublets, scvi_mod.resolvi_cluster = doublets, resolvi
    monkeypatch.setattr(ad_utils, "_calc_embeddings", embed)


def _obs(d, key):
    return sd.read_zarr(d.zarr)[key].obs


def _setup(d, **k):
    m.setup_adata_region(d, "cellpose_cell", suffix="_filt", **k)


def _doublets(d, **k):
    m.remove_doublets_region(d, "cellpose_cell", suffix="_filt", **k)


def _resolvi(d, **k):
    m.resolvi_cluster_region(d, "cellpose_cell", suffix="_filt", model_kwargs={}, **k)


# ---- steps sharing a table ----

def test_setup_doublets_resolvi_do_not_block_each_other(ds, stand_ins):
    _setup(ds)
    _doublets(ds)
    _resolvi(ds)
    obs = _obs(ds, FILT)
    assert {"base_leiden", "doublet_bool", "true_proportion", "resolvi_leiden", "corr_leiden"} <= set(obs.columns)
    assert len(obs) == N - N // 4  # setup kept pass_qc cells


def test_resolvi_then_doublets_order_free(ds, stand_ins):
    _setup(ds)
    _resolvi(ds)
    _doublets(ds)
    assert {"doublet_bool", "true_proportion"} <= set(_obs(ds, FILT).columns)


@pytest.mark.parametrize("step,match", [(_setup, f"tables/{FILT}"),
                                        (_doublets, "doublet_score, doublet_bool"),
                                        (_resolvi, "true_proportion")])
def test_each_rerun_raises_on_its_own_outputs(ds, stand_ins, step, match):
    _setup(ds)
    _doublets(ds)
    _resolvi(ds)
    before = _obs(ds, FILT).copy()
    with pytest.raises(FileExistsError, match=match):
        step(ds)
    pd.testing.assert_frame_equal(_obs(ds, FILT), before)


def test_rerun_with_overwrite(ds, stand_ins):
    _setup(ds)
    _doublets(ds)
    _doublets(ds, overwrite=True)
    _resolvi(ds)
    _resolvi(ds, overwrite=True)
    obs = _obs(ds, FILT)
    assert {"doublet_bool", "true_proportion"} <= set(obs.columns)
    assert not any("__tmp" in p for p in sd.read_zarr(ds.zarr).elements_paths_on_disk())


def test_setup_overwrite_warns_about_dropped_columns(ds, stand_ins, caplog):
    _setup(ds)
    _doublets(ds)
    _resolvi(ds)
    with caplog.at_level("WARNING"):
        _setup(ds, overwrite=True)
    assert "doublet_score, doublet_bool" in caplog.text and "true_proportion" in caplog.text
    assert "doublet_bool" not in _obs(ds, FILT).columns  # rebuilt from the filtered table


def test_setup_no_more_silent_fallback(ds, stand_ins, caplog):
    _setup(ds)
    with caplog.at_level("WARNING"):
        with pytest.raises(FileExistsError):
            _setup(ds)
    assert "Skipping write operation" not in caplog.text


# ---- call-region-tz ----

def test_call_region_tz_guard_and_overwrite(ds, monkeypatch, tmp_path):
    import spida.P.call_regions as cr

    calls = []
    monkeypatch.setattr(cr, "call_regions", lambda **k: calls.append(k))
    gpkg = tmp_path / "wm.gpkg"
    m.call_region_tz(ds, geoms_name="wm_region", save_geoms_path=str(gpkg))
    assert calls[-1]["overwrite"] is False
    # existing shapes, or an existing gpkg, block a re-run before region calling starts
    s = sd.read_zarr(ds.zarr)
    s["wm_region"] = ShapesModel.parse(s[K[SHAPES_KEY]].copy())
    s.write_element("wm_region")
    with pytest.raises(FileExistsError, match="shapes/wm_region"):
        m.call_region_tz(ds, geoms_name="wm_region")
    gpkg.write_text("x")
    with pytest.raises(FileExistsError, match="wm.gpkg"):
        m.call_region_tz(ds, geoms_name="other", save_geoms_path=str(gpkg))
    n = len(calls)
    m.call_region_tz(ds, geoms_name="wm_region", overwrite=True)
    assert len(calls) == n + 1 and calls[-1]["overwrite"] is True


def test_call_regions_real_write(tmp_path):
    """Real region calling (no stub) on synthetic points: first write, overwrite swap, and the new
    shapes in the caller's in-memory sdata (the --plot figures read them there; gold PU/THM2 2026-10-17).
    Refusing a re-run is the up-front guard, tested in test_call_region_tz_guard_and_overwrite."""
    from spatialdata.transformations import Identity
    from spida.P.call_regions import call_region_worker, call_regions

    rng = np.random.default_rng(0)
    n = 40000
    x, y = rng.uniform(0, 1000, n), rng.uniform(0, 1000, n)
    gene = np.where((x < 400) & (rng.uniform(size=n) < 0.5), "BCAS1", rng.choice(["A", "B", "C"], n))
    pts = PointsModel.parse(dd.from_pandas(pd.DataFrame({"x": x, "y": y, "gene": gene}), npartitions=1),
                            feature_key="gene", transformations={"global": Identity(), "pixel": Identity()})
    z = tmp_path / "tz.zarr"
    sd.SpatialData(points={"tx": pts}).write(z)
    call_regions(zarr_path=z, geoms_name="wm", points_key="tx", transfer_genes="BCAS1", hex_size=50)
    assert sd.read_zarr(z)["wm"].union_all().centroid.x < 400  # found the BCAS1-rich side

    s = sd.read_zarr(z)
    call_region_worker(sdata=s, geoms_name="wm", points_key="tx", transfer_genes=["BCAS1"], hex_size=50,
                       gmm_ncomp=3, overwrite=True)
    r = sd.read_zarr(z)
    assert "wm" in s.shapes and "wm" in r.shapes  # in memory and on disk
    assert not any("__tmp" in p for p in r.elements_paths_on_disk())


# ---- write-adata ----

def test_write_adata_guard(ds, tmp_path, monkeypatch):
    from spida.config import SpidaConfig
    import spida.io.layout as layout

    monkeypatch.setattr(m, "dataset_from", lambda cfg, e, r: ds)
    out = ds.table_h5ad(TABLE)
    m.write_adata(SpidaConfig(), "E", reg_name="R", prefix_names=["cellpose_cell"])
    assert out.is_file()
    t0 = out.stat().st_mtime_ns
    with pytest.raises(FileExistsError, match=out.name):
        m.write_adata(SpidaConfig(), "E", reg_name="R", prefix_names=["cellpose_cell"])
    # all-regions mode skips instead of raising
    monkeypatch.setattr(layout.TreeLayout, "units", lambda self, exp, **k: ["R"])
    m.write_adata(SpidaConfig(), "E", prefix_names=["cellpose_cell"])
    assert out.stat().st_mtime_ns == t0
    os.utime(out, ns=(t0 - 10**10, t0 - 10**10))  # coarse clock: make "newer" unambiguous
    t0 = out.stat().st_mtime_ns
    m.write_adata(SpidaConfig(), "E", reg_name="R", prefix_names=["cellpose_cell"], overwrite=True)
    assert out.stat().st_mtime_ns > t0


# ---- inputs checked before any work (phase 3e) ----

def test_setup_without_filter_stops_before_reading(ds, monkeypatch):
    import anndata as ad
    s = sd.read_zarr(ds.zarr)
    t = s[TABLE].copy()
    t.obs = t.obs.drop(columns=["pass_qc"])
    from spida.io.store import write_element
    write_element(ds.zarr, TABLE, t, overwrite=True)
    monkeypatch.setattr(m, "read_table", lambda *a: pytest.fail("read the table"))
    with pytest.raises(FileNotFoundError, match="obs: pass_qc"):
        _setup(ds)


def test_doublets_without_setup_stops(ds, monkeypatch):
    monkeypatch.setattr(m, "read_table", lambda *a: pytest.fail("read the table"))
    with pytest.raises(FileNotFoundError, match=FILT):
        _doublets(ds)
