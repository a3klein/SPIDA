"""Regression tests for re-running filter_cells_region on an already-filtered region.

`run_filtering` merges its computed feature frame back into `adata.obs`, so a
region that has been filtered once has BOTH the source and the renamed target
columns (EntityID and CELL_ID, center_x and CENTER_X, ...). `load_metadata` then
renamed EntityID -> CELL_ID a second time and polars raised

    DuplicateError: column 'CELL_ID' is duplicate

which killed job 1574 in the cpsam_v2 campaign. The practical effect was that a
region could be filtered exactly once and never re-filtered -- so changing a
cutoff and re-running was impossible for any region already published.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest

from spida._constants import DEFAULT_PRESET
from spida.P.filtering import Filter


def _adata(obs: pd.DataFrame) -> ad.AnnData:
    a = ad.AnnData(np.zeros((len(obs), 2), dtype="float32"))
    a.obs = obs
    return a


def _fresh_obs() -> pd.DataFrame:
    """obs as it looks straight out of load-segmentation-region."""
    return pd.DataFrame(
        {
            "EntityID": ["c1", "c2", "c3"],
            "center_x": [1.0, 2.0, 3.0],
            "center_y": [4.0, 5.0, 6.0],
            "volume": [200.0, 300.0, 400.0],
            "fov": [0, 0, 1],
        }
    )


def _refiltered_obs() -> pd.DataFrame:
    """obs after one completed run: renamed columns are merged back alongside."""
    df = _fresh_obs()
    df["CELL_ID"] = df["EntityID"]
    df["CENTER_X"] = df["center_x"]
    df["CENTER_Y"] = df["center_y"]
    df["pass_qc"] = True
    return df


def _filter(obs: pd.DataFrame) -> Filter:
    return Filter(_adata(obs), "exp", "reg", "seg", donor_name="donor",
                  PRESET=DEFAULT_PRESET)


def test_first_run_renames_normally():
    meta = _filter(_fresh_obs()).load_metadata()
    assert meta["CELL_ID"].to_list() == ["c1", "c2", "c3"]
    assert meta["CENTER_X"].to_list() == [1.0, 2.0, 3.0]


def test_rerun_on_already_filtered_obs_does_not_raise():
    """The actual 1574 failure: this used to raise polars DuplicateError."""
    meta = _filter(_refiltered_obs()).load_metadata()
    assert meta.height == 3


def test_rerun_keeps_the_existing_target_values():
    """Skipping the rename must not silently swap in different values."""
    obs = _refiltered_obs()
    meta = _filter(obs).load_metadata()
    assert meta["CELL_ID"].to_list() == ["c1", "c2", "c3"]
    assert meta["CENTER_X"].to_list() == [1.0, 2.0, 3.0]
    assert meta["volume"].to_list() == [200.0, 300.0, 400.0]


def test_identity_rename_is_not_treated_as_collision():
    """meta_map maps "volume" -> "volume"; skipping that would drop the column."""
    meta = _filter(_fresh_obs()).load_metadata()
    assert "volume" in meta.columns


def test_repeated_reruns_are_stable():
    """Filtering three times in a row must give the same metadata each time."""
    obs = _refiltered_obs()
    first = _filter(obs).load_metadata()
    for _ in range(2):
        again = _filter(obs).load_metadata()
        assert again.equals(first)


@pytest.mark.parametrize("missing", ["CENTER_X", "CELL_ID"])
def test_partial_prior_state(missing):
    """A half-written obs (only some targets present) must still resolve."""
    obs = _refiltered_obs().drop(columns=[missing])
    meta = _filter(obs).load_metadata()
    assert meta.height == 3
    assert meta["CELL_ID"].to_list() == ["c1", "c2", "c3"]
