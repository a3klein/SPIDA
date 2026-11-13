"""Filtering reads schema names; proseg-native centroid_x / centroid_y only through a deprecated fallback."""
import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from spida._constants import CELL_ID, CELL_X, CELL_Y, CELL_VOLUME
from spida.io.layout import RegionMeta
from spida.P.filtering import Filter

META = RegionMeta("exp", "reg", "proseg_cell", "donor")


def _adata(x_col, y_col, n=4):
    obs = pd.DataFrame({"EntityID": [str(i) for i in range(n)], x_col: np.arange(n, dtype=float),
                        y_col: np.arange(n, dtype=float) * 2, "volume": np.full(n, 100.0),
                        "centroid_z": np.zeros(n), "cluster": np.ones(n, dtype=int)},
                       index=[f"{i}.exp.reg.proseg_cell" for i in range(n)])
    return ad.AnnData(np.ones((n, 3), dtype=np.float32), obs=obs)


def test_schema_names_no_warning(caplog):
    with caplog.at_level(logging.WARNING):
        m = Filter(_adata("center_x", "center_y"), META).load_metadata()
    assert {CELL_ID, CELL_X, CELL_Y, CELL_VOLUME} <= set(m.columns)
    assert m[CELL_Y].to_list() == [0.0, 2.0, 4.0, 6.0]
    assert "fallback" not in caplog.text


def test_legacy_proseg_centroids_read_with_warning(caplog):
    with caplog.at_level(logging.WARNING):
        m = Filter(_adata("centroid_x", "centroid_y"), META).load_metadata()
    assert m[CELL_X].to_list() == [0.0, 1.0, 2.0, 3.0]
    assert "centroid_x" in caplog.text and "fallback will be removed" in caplog.text


def test_neither_name_still_fails():
    with pytest.raises(Exception, match=CELL_X):
        Filter(_adata("pos_x", "pos_y"), META).load_metadata()
