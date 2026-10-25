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
from spida.io.layout import RegionMeta


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
    return Filter(_adata(obs), RegionMeta("exp", "reg", "seg", "donor"),
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


# ---- re-filter keeps the same columns (phase 2, 2026-10-15) ----
# reset_index() used to turn the composite index into an extra `Index` column on every
# second run (17 -> 18 -> 17 columns); the old index is now dropped, since it is rebuilt.

CUTOFFS = {"volume_min": 100, "volume_max": 4000, "vol_max_mult": 3, "max_vol_by_median": False,
           "n_count_min": 20, "n_count_max": 2000, "n_count_max_mult": 5, "n_count_by_median": False,
           "n_gene_min": 5, "n_gene_max": 1000, "n_blank_min": -1, "n_blank_max": 5,
           "ncpv_min_quantile": 0.02, "ncpv_max_quantile": 1, "DAPI_quantile": 0.05}


def _loaded(n=200):
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    ids = [str(1000 + i) for i in range(n)]
    obs = pd.DataFrame({"EntityID": ids, "center_x": rng.uniform(0, 100, n), "center_y": rng.uniform(0, 100, n),
                        "volume": rng.uniform(200, 2000, n), "fov": 0, "DAPI_high_pass": rng.uniform(0, 1, n)},
                       index=pd.Index(ids, name="index"))
    a = ad.AnnData(sp.csr_matrix(rng.poisson(3, (n, 50)).astype(np.float32)), obs=obs)
    a.obsm["blank"] = rng.poisson(0.2, (n, 3))
    return a


def test_refilter_is_stable(tmp_path):
    import json
    from spida.P.filtering import run_filtering

    cut = tmp_path / "cutoffs.json"
    cut.write_text(json.dumps(CUTOFFS))
    run = lambda a: run_filtering(a.copy(), RegionMeta("E", "R", "cellpose_cell", "D"), cutoffs_path=cut)
    first = run(_loaded())
    second = run(first)
    third = run(second)
    for again in (second, third):
        assert list(again.obs.columns) == list(first.obs.columns)
        assert list(again.obs_names) == list(first.obs_names)
        assert (again.obs["pass_qc"] == first.obs["pass_qc"]).all()
    assert "Index" not in first.obs.columns


def test_filter_guard_columns_are_produced(tmp_path):
    """Drift: every column filter-cells' guard checks is one run_filtering actually writes (phase 3d)."""
    import json
    from spida.P.filtering import run_filtering
    from spida.P.main import _FILTER_COLUMNS

    cut = tmp_path / "cutoffs.json"
    cut.write_text(json.dumps(CUTOFFS))
    out = run_filtering(_loaded(), RegionMeta("E", "R", "cellpose_cell", "D"), cutoffs_path=cut)
    assert set(_FILTER_COLUMNS) <= set(out.obs.columns)
