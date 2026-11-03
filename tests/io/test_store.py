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
from spatialdata.models import Image2DModel, Labels2DModel, PointsModel, ShapesModel, TableModel

from spida.io.store import element_on_disk, is_writeable, update_obs, write_element

N = 30


def _sdata(rng):
    xy = rng.uniform(0, 100, (N, 2))
    shapes = gpd.GeoDataFrame({"EntityID": np.arange(N)}, geometry=[Point(x, y).buffer(2) for x, y in xy],
                              index=np.arange(N).astype(str))
    pts = pd.DataFrame({"x": rng.uniform(0, 100, 200), "y": rng.uniform(0, 100, 200),
                        "z": rng.uniform(0, 10, 200), "gene": rng.choice(list("ABCD"), 200)})
    obs = pd.DataFrame({"region": pd.Categorical(["cells"] * N), "EntityID": np.arange(N)}, index=np.arange(N).astype(str))
    table = TableModel.parse(ad.AnnData(rng.poisson(2, (N, 4)).astype(np.float32), obs=obs),
                             region="cells", region_key="region", instance_key="EntityID")
    return sd.SpatialData(
        images={"img": Image2DModel.parse(rng.integers(0, 255, (1, 128, 128)).astype(np.uint8),
                                          dims=("c", "y", "x"), scale_factors=[2], chunks=64)},
        labels={"lab": Labels2DModel.parse(rng.integers(0, 5, (128, 128)).astype(np.int32), dims=("y", "x"), chunks=64)},
        shapes={"cells": ShapesModel.parse(shapes)},
        points={"tx": PointsModel.parse(dd.from_pandas(pts, npartitions=2), coordinates={"x": "x", "y": "y", "z": "z"},
                                        feature_key="gene")},
        tables={"table": table},
    )


@pytest.fixture
def store(tmp_path):
    path = tmp_path / "s.zarr"
    _sdata(np.random.default_rng(0)).write(path)
    return path


def _mtimes(path, skip):
    # newest mtime per element, all but ``skip``
    out = {}
    for kind in ("images", "labels", "points", "shapes", "tables"):
        for el in (path / kind).iterdir():
            if el.is_dir() and el.name != skip:
                out[f"{kind}/{el.name}"] = max(f.stat().st_mtime for f in el.rglob("*") if f.is_file())
    return out


# ---- is_writeable ----

def test_is_writeable_new_outputs(store, tmp_path):
    assert is_writeable(store, ["new_table"], files=[tmp_path / "nope.parquet"])


def test_is_writeable_lists_every_existing_output(store, tmp_path):
    f = tmp_path / "polygons.parquet"
    f.write_text("x")
    with pytest.raises(FileExistsError) as e:
        is_writeable(store, ["table", "cells", "new"], files=[f])
    msg = str(e.value)
    assert "tables/table" in msg and "shapes/cells" in msg and str(f) in msg and "new" not in msg.replace("--overwrite", "")


def test_is_writeable_overwrite(store):
    assert is_writeable(store, ["table"], overwrite=True)


def test_is_writeable_empty_dir_is_free(tmp_path):
    (tmp_path / "seg").mkdir()
    assert is_writeable(files=[tmp_path / "seg"])
    (tmp_path / "seg" / "a.csv").write_text("x")
    with pytest.raises(FileExistsError):
        is_writeable(files=[tmp_path / "seg"])


def test_is_writeable_missing_store(tmp_path):
    assert is_writeable(tmp_path / "absent.zarr", ["table"])


# ---- write_element ----

def test_write_new_element(store):
    s = sd.read_zarr(store)
    write_element(store, "cells2", ShapesModel.parse(s["cells"].copy()))
    assert "cells2" in sd.read_zarr(store).shapes


def test_write_existing_without_overwrite_raises_and_keeps(store):
    before = _mtimes(store, skip=None)
    t = sd.read_zarr(store)["table"].copy()
    t.obs["new"] = 1
    with pytest.raises(FileExistsError):
        write_element(store, "table", t)
    assert _mtimes(store, skip=None) == before
    assert "new" not in sd.read_zarr(store)["table"].obs


def _new_table(s):
    t = s["table"].copy()
    t.obs["new"] = 1
    return t


def _new_points(s):  # lazily backed by the files being replaced
    p = s["tx"]
    p["gene2"] = p["gene"]
    return PointsModel.parse(p)


def _new_shapes(s):
    g = s["cells"].copy()
    g["new"] = 1
    return ShapesModel.parse(g)


def _new_image(s):
    return Image2DModel.parse(np.full((1, 128, 128), 7, np.uint8), dims=("c", "y", "x"), scale_factors=[2], chunks=64)


def _new_labels(s):
    return Labels2DModel.parse(np.ones((128, 128), np.int32), dims=("y", "x"), chunks=64)


def _check(r, key):
    if key == "table":
        return "new" in r["table"].obs and r["table"].n_obs == N
    if key == "tx":
        return "gene2" in r["tx"].columns and len(r["tx"]) == 200
    if key == "cells":
        return "new" in r["cells"].columns
    if key == "img":
        return list(r["img"].keys()) == ["scale0", "scale1"] and int(r["img"]["scale0"]["image"].values.min()) == 7
    return int(r["lab"].values.min()) == 1


@pytest.mark.parametrize("key,make", [("table", _new_table), ("tx", _new_points), ("cells", _new_shapes),
                                      ("img", _new_image), ("lab", _new_labels)])
def test_overwrite_swaps_in_place(store, key, make):
    others = _mtimes(store, skip=key)
    write_element(store, key, make(sd.read_zarr(store)), overwrite=True)
    r = sd.read_zarr(store)
    assert _check(r, key)
    assert element_on_disk(store, f"{key}__tmp") is None
    assert _mtimes(store, skip=key) == others  # nothing else touched


def test_stale_tmp_is_removed(store):
    s = sd.read_zarr(store)
    s["cells__tmp"] = ShapesModel.parse(s["cells"].copy())
    s.write_element("cells__tmp")
    write_element(store, "cells", _new_shapes(sd.read_zarr(store)), overwrite=True)
    r = sd.read_zarr(store)
    assert "cells__tmp" not in r.shapes and "new" in r["cells"].columns


# ---- reads ----

def test_read_table_matches_spatialdata_read(store):
    from spida.io.store import open_store, read_table
    a, b = open_store(store)["table"], read_table(store, "table")
    pd.testing.assert_frame_equal(a.obs, b.obs)
    assert np.array_equal(np.asarray(a.X), np.asarray(b.X))
    assert a.uns["spatialdata_attrs"] == b.uns["spatialdata_attrs"]


@pytest.mark.parametrize("how", ["write_element", "update_obs"])
def test_read_table_round_trip(store, how):
    """A table read with anndata, changed and written back, is a valid spatialdata table (P write-backs)."""
    from spida.io.store import read_table
    t = read_table(store, "table")
    t.obs["label"] = pd.Categorical(np.where(np.arange(N) % 2, "a", "b"))
    if how == "write_element":
        write_element(store, "table", t, overwrite=True)
    else:
        update_obs(store, "table", t.obs)
    r = sd.read_zarr(store)
    TableModel().validate(r["table"])
    assert list(r["table"].obs["label"]) == list(t.obs["label"])
    assert len(r["cells"]) == N and len(r["tx"]) == 200


def test_write_element_new_table_from_read_table(store):
    """setup_adata pattern: a subset of a read table written under a new key."""
    from spida.io.store import read_table
    t = read_table(store, "table")
    write_element(store, "table_filt", t[t.obs["EntityID"] % 2 == 0].copy())
    r = sd.read_zarr(store)
    TableModel().validate(r["table_filt"])
    assert r["table_filt"].n_obs == N // 2 and r["table_filt"].uns["spatialdata_attrs"]["region"] == "cells"


# ---- update_obs ----

def test_update_obs_only_touches_obs(store):
    t0 = sd.read_zarr(store)["table"]
    obs = t0.obs.copy()
    obs["label"] = pd.Categorical(np.where(np.arange(N) % 2, "a", "b"))
    update_obs(store, "table", obs)
    t = sd.read_zarr(store)["table"]
    assert list(t.obs["label"]) == list(obs["label"])
    assert np.array_equal(np.asarray(t.X), np.asarray(t0.X))
    assert t.uns["spatialdata_attrs"] == t0.uns["spatialdata_attrs"]
    r = sd.read_zarr(store)
    assert len(r["tx"]) == 200 and len(r["cells"]) == N and "img" in r.images and "lab" in r.labels


def test_update_obs_index_must_match(store):
    obs = sd.read_zarr(store)["table"].obs.iloc[::-1].copy()
    with pytest.raises(ValueError, match="index"):
        update_obs(store, "table", obs)


def test_update_obs_keeps_link_columns(store):
    obs = sd.read_zarr(store)["table"].obs.drop(columns=["EntityID"])
    with pytest.raises(ValueError, match="EntityID"):
        update_obs(store, "table", obs)


def test_update_obs_missing_table(store):
    with pytest.raises(KeyError):
        update_obs(store, "nope", pd.DataFrame())


# ---- column guard ----

def test_obs_columns_reads_attrs(store):
    from spida.io.store import obs_columns
    assert obs_columns(store, "table") == list(sd.read_zarr(store)["table"].obs.columns)
    assert obs_columns(store, "nope") == []


def test_is_writeable_columns(store):
    with pytest.raises(FileExistsError, match="tables/table obs: EntityID"):
        is_writeable(store, columns={"table": ["EntityID", "pass_qc"]})
    assert is_writeable(store, columns={"table": ["pass_qc"], "nope": ["EntityID"]})
    assert is_writeable(store, columns={"table": ["EntityID"]}, overwrite=True)


def test_write_table_swaps(store):
    from spida.io.store import write_table
    t = sd.read_zarr(store)["table"].copy()
    t.obs["new"] = 1
    write_table(store, "table", t, overwrite=True)
    r = sd.read_zarr(store)
    assert "new" in r["table"].obs and "table__tmp" not in r.tables


# ---- require (phase 3e) ----

def test_require_lists_every_missing_input(store, tmp_path):
    from spida.io.store import require
    f = tmp_path / "present.csv"
    f.write_text("x")
    assert require(store, ["table", "cells"], files=[f], columns={"table": ["EntityID"]})
    with pytest.raises(FileNotFoundError) as e:
        require(store, ["table", "nope"], files=[f, tmp_path / "absent.csv"],
                columns={"table": ["EntityID", "pass_qc"], "gone": ["x"]})
    msg = str(e.value)
    assert "nope" in msg and "absent.csv" in msg and "tables/table obs: pass_qc" in msg and "tables/gone" in msg
    assert "present.csv" not in msg and "EntityID" not in msg


def test_require_empty_dir_is_missing(tmp_path):
    from spida.io.store import require
    (tmp_path / "images").mkdir()
    with pytest.raises(FileNotFoundError, match="images"):
        require(files=[tmp_path / "images"])
