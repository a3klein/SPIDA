"""
Fast unit tests for cluster markers: call_markers_nsforest, summarize_nsforest_results(_from_uns),
call_markers dispatch, and the setup hook setup_adata._cluster_markers.
NSForest runs for real on a tiny table (few trees, one worker); failure paths are mocked.
"""
import logging
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

ns = pytest.importorskip("nsforest")

from spida.io.anndata import read_adata, write_adata  # noqa: E402
from spida.utilities import degs  # noqa: E402
from spida.utilities.degs import (  # noqa: E402
    call_markers,
    call_markers_nsforest,
    summarize_nsforest_results,
    summarize_nsforest_results_from_uns,
)

FAST = dict(n_trees=20, n_jobs=1)


def _adata(labels=("A", "B", "C", "D"), n_per=30, n_genes=24, seed=0):
    """Log-normalized; label i is marked by genes 3i..3i+2."""
    rng = np.random.default_rng(seed)
    lam = np.full((len(labels), n_genes), 0.3)
    for i in range(len(labels)):
        lam[i, 3 * i:3 * i + 3] = 8.0
    X = np.vstack([rng.poisson(lam[i], size=(n_per, n_genes)) for i in range(len(labels))]).astype(float)
    X = np.log1p(X / X.sum(1, keepdims=True) * 1e3)
    obs = pd.DataFrame({"ct": np.repeat(labels, n_per)}, index=[f"c{i}" for i in range(len(X))])
    return ad.AnnData(sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))


@pytest.fixture(scope="module")
def result():
    a = _adata()
    return a, call_markers_nsforest(a, "ct", uns_key="nsforest_markers", **FAST)


# ---- call_markers_nsforest ----

def test_finds_planted_markers(result):
    a, r = result
    assert set(r) == {"A", "B", "C", "D"}
    for i, ct in enumerate("ABCD"):
        nsf = r[ct].loc[r[ct]["is_nsforest_marker"], "names"]
        assert len(nsf) >= 1 and set(nsf) <= {f"g{3 * i}", f"g{3 * i + 1}", f"g{3 * i + 2}"}


def test_padding_and_columns(result):
    _, r = result
    for df in r.values():
        assert len(df) <= 5 and list(df["rank"]) == list(range(1, len(df) + 1))
        assert df["is_nsforest_marker"].iloc[: df["is_nsforest_marker"].sum()].all()  # NSForest genes first
        assert {"names", "cell_type", "cluster_size", "f_score", "precision", "recall", "onTarget"} <= set(df)


def test_uns_storage_and_input_untouched(result):
    a, r = result
    stacked = a.uns["nsforest_markers"]
    assert len(stacked) == sum(len(v) for v in r.values())
    assert a.uns["nsforest_markers_skipped"] == []
    assert not a.varm  # NSForest worked on a copy


def test_min_cells_drops_small_cluster():
    a = ad.concat([_adata(), _adata(labels=("E",), n_per=5, seed=1)], index_unique="-")
    assert "E" not in call_markers_nsforest(a, "ct", min_cells=10, **FAST)


def test_skipped_cluster_recorded(caplog):
    a = _adata()
    real = ns.nsforesting.NSForest

    def fail_b(adata, col, cluster_list, **kw):
        if cluster_list == ["B"]:
            raise ValueError("at least one array or dtype is required")
        return real(adata, col, cluster_list=cluster_list, **kw)

    with patch.object(ns.nsforesting, "NSForest", side_effect=fail_b), caplog.at_level(logging.WARNING):
        r = call_markers_nsforest(a, "ct", uns_key="m", **FAST)
    assert "B" not in r and a.uns["m_skipped"] == ["B"]
    assert "skipped B" in caplog.text


def test_all_clusters_fail_raises():
    with patch.object(ns.nsforesting, "NSForest", side_effect=ValueError("no genes")):
        with pytest.raises(RuntimeError, match="no results"):
            call_markers_nsforest(_adata(), "ct", **FAST)


def test_gene_selection_default_and_passthrough():
    seen = []
    real = ns.nsforesting.NSForest

    def spy(*a, **kw):
        seen.append(kw["gene_selection"])
        return real(*a, **kw)

    with patch.object(ns.nsforesting, "NSForest", side_effect=spy):
        call_markers_nsforest(_adata(), "ct", **FAST)
        call_markers_nsforest(_adata(), "ct", gene_selection="BinaryFirst_high", **FAST)
    assert set(seen[:4]) == {"BinaryFirst_moderate"} and set(seen[4:]) == {"BinaryFirst_high"}


def test_gene_selection_rejects_unknown():
    # NSForest documents "BinaryFirst_low" but silently falls back to "high" on it
    with pytest.raises(ValueError, match="gene_selection must be one of"):
        call_markers_nsforest(_adata(), "ct", gene_selection="BinaryFirst_low", **FAST)


def test_bad_inputs():
    with pytest.raises(ValueError, match="not found"):
        call_markers_nsforest(_adata(), "nope", **FAST)
    with pytest.raises(ValueError, match="at least 2"):
        call_markers_nsforest(_adata(labels=("A",)), "ct", **FAST)


# ---- summaries ----

def test_summary_and_h5ad_round_trip(result, tmp_path):
    a, r = result
    direct = summarize_nsforest_results(r, top_n=3)
    write_adata(tmp_path / "t.h5ad", a)
    from_uns = summarize_nsforest_results_from_uns(read_adata(tmp_path / "t.h5ad"), top_n=3)
    pd.testing.assert_frame_equal(direct.sort_values("cell_type").reset_index(drop=True),
                                  from_uns.sort_values("cell_type").reset_index(drop=True))
    assert all(len(x) <= 3 for x in direct["top_3_markers"])


def test_summary_nsforest_only(result):
    _, r = result
    s = summarize_nsforest_results(r, nsforest_only=True)
    for _, row in s.iterrows():
        assert row["top_5_markers"] == row["nsforest_markers"][:5]


def test_summary_from_uns_missing_key():
    with pytest.raises(KeyError):
        summarize_nsforest_results_from_uns(_adata(), uns_key="nope")


# ---- call_markers dispatch ----

@pytest.mark.parametrize("method", ["nsforest", "scran", "scanpy"])
def test_dispatch_default_uns_key(method):
    fn = {"nsforest": "call_markers_nsforest", "scran": "call_degs_scran", "scanpy": "call_degs_scanpy"}[method]
    with patch.object(degs, fn, return_value={}) as m:
        call_markers(_adata(), "ct", method=method, foo=1)
    m.assert_called_once()
    assert m.call_args.kwargs["uns_key"] == f"{method}_markers" and m.call_args.kwargs["foo"] == 1


def test_dispatch_bad_method():
    with pytest.raises(ValueError, match="method must be one of"):
        call_markers(_adata(), "ct", method="deseq")


# ---- setup hook ----

def test_cluster_markers_hook(caplog):
    from spida.P.setup_adata import _cluster_markers

    with patch.object(degs, "call_markers") as m:
        _cluster_markers(_adata(), "ct", None)
        m.assert_not_called()
        _cluster_markers(_adata(), "ct", "nsforest", n_jobs=3)
        assert m.call_args.kwargs == {"method": "nsforest", "n_jobs": 3, "celltype_col": "ct"}
        _cluster_markers(_adata(), "ct", "scran", n_jobs=3)
        assert "n_jobs" not in m.call_args.kwargs  # only nsforest takes n_jobs
    with patch.object(degs, "call_markers", side_effect=RuntimeError("boom")), caplog.at_level(logging.WARNING):
        _cluster_markers(_adata(), "ct", "nsforest")  # logged, not raised
    assert "nsforest markers on ct failed: boom" in caplog.text
