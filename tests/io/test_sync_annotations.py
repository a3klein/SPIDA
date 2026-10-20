"""sync-annotations: annotation columns from the h5ad copy into the zarr table's obs (phase 2, step 3)."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")
sd = pytest.importorskip("spatialdata")

import anndata as ad
import geopandas as gpd
from click.testing import CliRunner
from shapely.geometry import Point
from spatialdata.models import ShapesModel, TableModel

from spida.io.layout import Dataset, element_keys
from spida._constants import SHAPES_KEY, TABLE_KEY
from spida.I._io import sync_annotations

K = element_keys("cellpose_cell", None)
KEY = K[TABLE_KEY] + "_filt"
N = 30


@pytest.fixture
def ds(tmp_path):
    rng = np.random.default_rng(0)
    idx = np.arange(N).astype(str)
    shapes = gpd.GeoDataFrame({"EntityID": np.arange(N)}, geometry=[Point(i, i).buffer(1) for i in range(N)], index=idx)
    obs = pd.DataFrame({"cells_region": pd.Categorical([K[SHAPES_KEY]] * N), "EntityID": np.arange(N),
                        "base_leiden": pd.Categorical(["0", "1"] * (N // 2))}, index=idx)
    t = ad.AnnData(rng.poisson(3, (N, 4)).astype(np.float32), obs=obs)
    t.layers["raw"] = t.X.copy()
    t.obsm["spatial"] = rng.uniform(size=(N, 2))
    table = TableModel.parse(t, region=K[SHAPES_KEY], region_key="cells_region", instance_key="EntityID")
    s = sd.SpatialData(shapes={K[SHAPES_KEY]: ShapesModel.parse(shapes)}, tables={KEY: table})
    s.attrs["spida"] = {"namespace": "", "exp": None, "reg": None}
    d = Dataset(zarr=tmp_path / "s.zarr", anndata=tmp_path / "ad")
    s.write(d.zarr)
    copy = sd.read_zarr(d.zarr)[KEY].copy()  # the h5ad copy, as load_table makes it
    copy.obs["mmc_AIT_Subclass"] = pd.Categorical(np.where(np.arange(N) % 3, "Astro", "Oligo"))
    copy.obs["mmc_AIT_Subclass_transfer_score"] = np.linspace(0, 1, N)
    copy.obs["allcools_HMBA_Subclass"] = "D1"
    copy.obs["notes"] = "x"  # not an annotation column
    d.table_h5ad(KEY).parent.mkdir(parents=True)
    copy.write_h5ad(d.table_h5ad(KEY))
    return d


def _table(d):
    return sd.read_zarr(d.zarr)[KEY]


def test_first_sync_obs_only(ds):
    before = _table(ds)
    cols = sync_annotations(ds, "cellpose_cell")
    assert cols == ["mmc_AIT_Subclass", "mmc_AIT_Subclass_transfer_score", "allcools_HMBA_Subclass"]
    t = _table(ds)
    assert "notes" not in t.obs and list(t.obs["mmc_AIT_Subclass"]) == list(ad.read_h5ad(ds.table_h5ad(KEY)).obs["mmc_AIT_Subclass"])
    assert np.array_equal(t.X, before.X) and np.array_equal(t.layers["raw"], before.layers["raw"])
    assert np.array_equal(t.obsm["spatial"], before.obsm["spatial"])
    assert t.uns["spatialdata_attrs"] == before.uns["spatialdata_attrs"]
    pd.testing.assert_frame_equal(t.obs[before.obs.columns], before.obs)


def test_guard_and_overwrite(ds):
    sync_annotations(ds, "cellpose_cell", reference="AIT")
    with pytest.raises(FileExistsError, match="mmc_AIT_Subclass"):
        sync_annotations(ds, "cellpose_cell", reference="AIT")
    sync_annotations(ds, "cellpose_cell", reference="HMBA")  # another reference is not blocked
    c = ad.read_h5ad(ds.table_h5ad(KEY))
    c.obs["mmc_AIT_Subclass"] = pd.Categorical(["OPC"] * N)
    c.write_h5ad(ds.table_h5ad(KEY))
    sync_annotations(ds, "cellpose_cell", reference="AIT", overwrite=True)
    assert set(_table(ds).obs["mmc_AIT_Subclass"]) == {"OPC"}


def test_explicit_columns(ds):
    assert sync_annotations(ds, "cellpose_cell", columns=["notes"]) == ["notes"]
    assert "notes" in _table(ds).obs and "mmc_AIT_Subclass" not in _table(ds).obs
    with pytest.raises(KeyError, match="nope"):
        sync_annotations(ds, "cellpose_cell", columns=["nope"])


def test_reordered_and_subset_adata_file(ds, tmp_path, caplog):
    c = ad.read_h5ad(ds.table_h5ad(KEY))
    sub = c[c.obs_names[::-1][:20]].copy()  # reversed order, 10 cells missing
    sub.obs["mmc_X_flag"] = True
    f = tmp_path / "sub.h5ad"
    sub.write_h5ad(f)
    with caplog.at_level("WARNING"):
        sync_annotations(ds, "cellpose_cell", adata_file=f, reference="X")
    assert "10 cells" in caplog.text
    obs = _table(ds).obs
    assert obs["mmc_X_flag"].isna().sum() == 10 and obs.loc[sub.obs_names, "mmc_X_flag"].all()


def test_extra_cells_raise(ds, tmp_path):
    c = ad.read_h5ad(ds.table_h5ad(KEY))
    c.obs_names = [f"new{i}" for i in range(N)]
    f = tmp_path / "other.h5ad"
    c.write_h5ad(f)
    with pytest.raises(ValueError, match="not in the zarr table"):
        sync_annotations(ds, "cellpose_cell", adata_file=f)


def test_nothing_to_sync(ds):
    with pytest.raises(ValueError, match="no annotation columns"):
        sync_annotations(ds, "cellpose_cell", reference="absent")


def test_cli_and_experiment_skip(ds, tmp_path, monkeypatch):
    from spida.I import cli as i_cli
    import spida.io.layout as layout
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({}))
    args = ["--config", str(cfg), "sync-annotations", "cellpose_cell", "--zarr_path", str(ds.zarr),
            "--anndata_path", str(ds.anndata), "--reference", "AIT"]
    assert CliRunner().invoke(i_cli.cli, args).exit_code == 0
    assert isinstance(CliRunner().invoke(i_cli.cli, args).exception, FileExistsError)
    monkeypatch.setattr(layout.TreeLayout, "units", lambda self, exp, **k: ["R"])
    monkeypatch.setattr(layout.TreeLayout, "dataset", lambda self, exp, reg, **k: ds)
    res = CliRunner().invoke(i_cli.cli, ["--config", str(cfg), "sync-annotations-experiment", "cellpose_cell",
                                         "--exp", "E", "--reference", "AIT"])
    assert res.exit_code == 0, res.output  # already synced -> skipped, not an error


def test_backup_to_spatialdata_off_by_default():
    from spida.I import cli as i_cli
    opt = next(p for p in i_cli.cli.commands["allcools-integration-region"].params if p.name == "backup_to_spatialdata")
    assert opt.default is False
