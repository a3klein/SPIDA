"""Tests for the Stage-3 spec-driven segmentation loader helpers
(``spida.S.io._seg_readers`` + the shared tail in ``ingest_exp``).

These exercise the reader logic and the EntityID-keying contract without a full
SpatialData zarr / merscope read: get_polygons, get_table, get_table_from_adata,
the multipolygon dedup, and canonical/legacy filename resolution.
"""

import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import anndata as ad
import pytest
from shapely.geometry import Polygon, MultiPolygon

import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity

from spida.S.io._seg_readers import (
    get_polygons, get_table, get_table_from_adata, _CELL_KEY,
)
from spida.S.io.ingest_exp import _resolve_file, _cast_multipolygons_to_polygons

DATA = Path(__file__).parent / "data" / "segmentation_utils"
TX = {"global": Identity()}


# --------------------------------------------------------------------------- #
# _resolve_file
# --------------------------------------------------------------------------- #
def test_resolve_prefers_canonical(tmp_path):
    (tmp_path / "canon.csv").write_text("x")
    (tmp_path / "legacy.csv").write_text("x")
    assert _resolve_file(tmp_path, "canon.csv", "legacy.csv") == tmp_path / "canon.csv"


def test_resolve_falls_back_to_legacy(tmp_path):
    (tmp_path / "legacy.csv").write_text("x")
    assert _resolve_file(tmp_path, "canon.csv", "legacy.csv") == tmp_path / "legacy.csv"


def test_resolve_returns_canonical_when_missing(tmp_path):
    # neither exists -> canonical path (so downstream errors name the expected file)
    assert _resolve_file(tmp_path, "canon.csv", "legacy.csv") == tmp_path / "canon.csv"
    assert _resolve_file(tmp_path, "canon.csv", None) == tmp_path / "canon.csv"


# --------------------------------------------------------------------------- #
# get_polygons
# --------------------------------------------------------------------------- #
def test_get_polygons_keyed_on_entityid():
    shapes = get_polygons(DATA / "sum_signals_boundaries.parquet", TX, cell_id_col=_CELL_KEY)
    assert shapes.index.name == _CELL_KEY
    assert all(isinstance(i, str) for i in shapes.index[:50])
    assert shapes.geometry.name == "geometry"           # renamed from "Geometry"
    assert "z" not in shapes.columns and "z_layer" in shapes.columns  # z avoided -> z_layer
    assert "global" in shapes.attrs["transform"]                       # is a parsed ShapesModel


def test_get_polygons_geometry_valid_multipolygon():
    shapes = get_polygons(DATA / "sum_signals_boundaries.parquet", TX, cell_id_col=_CELL_KEY)
    assert shapes.geometry.is_valid.all()
    assert set(shapes.geometry.geom_type.unique()) <= {"MultiPolygon"}


# --------------------------------------------------------------------------- #
# get_table (CSV counts + obs)
# --------------------------------------------------------------------------- #
def _write_cbg_and_obs(tmp_path):
    cbg = pd.DataFrame(
        {"GeneA": [3, 0, 5], "GeneB": [1, 2, 0], "Blank-1": [7, 7, 7]},
        index=pd.Index([0, 1, 2], name=_CELL_KEY),
    )
    cbg.to_csv(tmp_path / "cell_by_gene.csv")
    obs = pd.DataFrame(
        {"volume": [10.0, 20.0, 30.0], "DAPI_raw": [1.0, 2.0, 3.0],
         "center_x": [0.0, 1.0, 2.0], "center_y": [0.0, 1.0, 2.0]},
        index=pd.Index([0, 1, 2], name=_CELL_KEY),
    )
    obs.to_csv(tmp_path / "cell_metadata.csv")
    return tmp_path / "cell_by_gene.csv", tmp_path / "cell_metadata.csv"


def test_get_table_entityid_and_blanks(tmp_path):
    cbg_path, obs_path = _write_cbg_and_obs(tmp_path)
    table = get_table(cbg_path, obs_path, "reg", "exp", "exp_reg", "shapes_k",
                      cell_id_col=_CELL_KEY)
    attrs = table.uns["spatialdata_attrs"]
    assert attrs["instance_key"] == _CELL_KEY
    assert attrs["region_key"] == "cells_region"
    assert list(table.var_names) == ["GeneA", "GeneB"]       # Blank moved out of X
    assert "blank" in table.obsm and table.obsm["blank"].shape[1] == 1
    assert "spatial" in table.obsm
    assert "volume" in table.obs.columns and "DAPI_raw" in table.obs.columns
    assert all(isinstance(i, str) for i in table.obs_names)


# --------------------------------------------------------------------------- #
# get_table_from_adata (proseg zarr path)
# --------------------------------------------------------------------------- #
def _fake_counts_adata(n_cells=4):
    X = np.arange(n_cells * 3, dtype=float).reshape(n_cells, 3)
    var = pd.DataFrame({"gene": ["GeneA", "GeneB", "Blank-1"]},
                       index=["GeneA", "GeneB", "Blank-1"])
    obs = pd.DataFrame({"cell": np.arange(n_cells),
                        "centroid_x": np.arange(n_cells, dtype=float),
                        "centroid_y": np.arange(n_cells, dtype=float)},
                       index=[str(i) for i in range(n_cells)])
    import scipy.sparse as sp
    return ad.AnnData(sp.csr_matrix(X), var=var, obs=obs)


def test_get_table_from_adata_native_obs():
    adata = _fake_counts_adata()
    table = get_table_from_adata(adata, None, "reg", "exp", "exp_reg", "shapes_k",
                                 cell_id_col=_CELL_KEY)
    attrs = table.uns["spatialdata_attrs"]
    assert attrs["instance_key"] == _CELL_KEY
    assert list(table.var_names) == ["GeneA", "GeneB"]       # Blank-1 excluded
    assert "transcript_count" in table.obs.columns
    assert "blank" in table.obsm
    # transcript_count == row sum over non-blank genes
    expected = adata[:, ["GeneA", "GeneB"]].X.toarray().sum(axis=1)
    assert np.allclose(table.obs["transcript_count"].values, expected)


def test_get_table_from_adata_obs_override_alignment():
    adata = _fake_counts_adata(n_cells=4)
    # obs override in REVERSED cell order; loader passes it pre-aligned to cell order
    order = adata.obs["cell"].astype(int).values
    meta = pd.DataFrame({"volume": order.astype(float)},
                        index=pd.Index(order[::-1], name=_CELL_KEY)).reindex(order)
    table = get_table_from_adata(adata, meta, "reg", "exp", "exp_reg", "shapes_k",
                                 cell_id_col=_CELL_KEY)
    # label-based alignment: cell i keeps the volume meta held at label i.
    # meta was built reversed (label i -> volume order[-1-i]), so after reindex(order)
    # the per-row volume is order[::-1].
    assert np.allclose(table.obs["volume"].values, order[::-1].astype(float))


# --------------------------------------------------------------------------- #
# _cast_multipolygons_to_polygons dedup
# --------------------------------------------------------------------------- #
def test_cast_dedup_keeps_largest_per_entityid():
    big = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])       # area 100
    tiny = Polygon([(20, 20), (20, 21), (21, 21)])            # area 0.5
    gdf = gpd.GeoDataFrame(
        {_CELL_KEY: [1], "geometry": [MultiPolygon([big, tiny])]}, geometry="geometry"
    )
    gdf.index = gdf[_CELL_KEY].astype(str)
    shapes = ShapesModel.parse(gdf, transformations=TX)
    sdata = sd.SpatialData(shapes={"k": shapes})
    out = _cast_multipolygons_to_polygons(sdata, "k", subset_field=[_CELL_KEY])
    g = out["k"]
    assert len(g) == 1
    assert g.geometry.area.iloc[0] == pytest.approx(100.0)
