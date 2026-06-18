"""
Tests for spida.I.allcools private helpers.

_allcools_cef
  - ref-only mode (qry_col_level=None): returns expected gene index, calls CEF once
  - joint mode (qry_col_level set): returns union of ref and qry enriched genes
  - joint result is a superset of ref-only result
  - correct column name is passed to cluster_enriched_features for each dataset

_clust2clust_transfer  (proportion-threshold method)
  - assigns dominant ref label when proportion >= threshold and counts pass
  - returns "unknown" when proportion < threshold
  - returns "unknown" when ref_cell_count < threshold
  - returns "unknown" when qry_cell_count < threshold
  - saves full df_map to adata_comb.uns before filtering
  - output column is categorical

_allcools_clust2clust_transfer  (ALLCools confusion-matrix method)
  - assigns correct label to query cells for unambiguous clusters (1 ref type per group)
  - leaves query cells in small clusters as "unknown"
  - populates adata_comb.uns with integration results
"""

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp
from unittest.mock import patch


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ref_adata():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 20
    X = sp.csr_matrix(rng.poisson(5, size=(n_cells, n_genes)).astype(float))
    obs = pd.DataFrame(
        {"ref_type": np.repeat(["TypeA", "TypeB", "TypeC"], 20)},
        index=[f"ref_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def qry_adata():
    rng = np.random.default_rng(1)
    n_cells, n_genes = 40, 20
    X = sp.csr_matrix(rng.poisson(5, size=(n_cells, n_genes)).astype(float))
    obs = pd.DataFrame(
        {"qry_cluster": np.repeat(["ClustA", "ClustB"], 20)},
        index=[f"qry_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------

def _make_static_mock(enriched_indices):
    """Return a mock that marks the given gene indices as enriched."""
    def _mock(adata, cluster_col, top_n, alpha, stat_plot, method):
        flags = pd.Series(False, index=adata.var_names)
        for i in enriched_indices:
            if i < len(adata.var_names):
                flags.iloc[i] = True
        adata.var[f"{cluster_col}_enriched_features"] = flags
    return _mock


def _make_sequential_mock(ref_indices, qry_indices):
    """Return a mock that marks ref_indices on the first call and qry_indices on the second."""
    call_count = {"n": 0}

    def _mock(adata, cluster_col, top_n, alpha, stat_plot, method):
        indices = ref_indices if call_count["n"] == 0 else qry_indices
        flags = pd.Series(False, index=adata.var_names)
        for i in indices:
            if i < len(adata.var_names):
                flags.iloc[i] = True
        adata.var[f"{cluster_col}_enriched_features"] = flags
        call_count["n"] += 1

    return _mock


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestAllcoolsCef:

    @pytest.fixture(autouse=True)
    def _skip_if_no_allcools(self):
        pytest.importorskip("ALLCools")

    def test_ref_only_returns_correct_genes(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        mock = _make_static_mock([0, 1, 2, 5])
        with patch("ALLCools.clustering.cluster_enriched_features", mock):
            result = _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level=None,
                top_n=4,
            )

        assert set(result) == {f"gene_{i}" for i in [0, 1, 2, 5]}

    def test_ref_only_calls_cef_once(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        call_count = {"n": 0}

        def counting_mock(adata, cluster_col, top_n, alpha, stat_plot, method):
            call_count["n"] += 1
            flags = pd.Series(True, index=adata.var_names)
            adata.var[f"{cluster_col}_enriched_features"] = flags

        with patch("ALLCools.clustering.cluster_enriched_features", counting_mock):
            _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level=None,
                top_n=5,
            )

        assert call_count["n"] == 1

    def test_joint_returns_union(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        mock = _make_sequential_mock(ref_indices=[0, 1, 2], qry_indices=[3, 4])
        with patch("ALLCools.clustering.cluster_enriched_features", mock):
            result = _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level="qry_cluster",
                top_n=3,
            )

        assert set(result) == {f"gene_{i}" for i in [0, 1, 2, 3, 4]}

    def test_joint_is_superset_of_ref_only(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        with patch("ALLCools.clustering.cluster_enriched_features",
                   _make_sequential_mock([0, 1, 2], [3, 4])):
            joint = _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level="qry_cluster",
                top_n=3,
            )

        with patch("ALLCools.clustering.cluster_enriched_features",
                   _make_static_mock([0, 1, 2])):
            ref_only = _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level=None,
                top_n=3,
            )

        assert set(ref_only).issubset(set(joint))
        assert len(joint) >= len(ref_only)

    def test_joint_uses_correct_column_names(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        seen_cols = []

        def recording_mock(adata, cluster_col, top_n, alpha, stat_plot, method):
            seen_cols.append(cluster_col)
            flags = pd.Series(False, index=adata.var_names)
            flags.iloc[0] = True
            adata.var[f"{cluster_col}_enriched_features"] = flags

        with patch("ALLCools.clustering.cluster_enriched_features", recording_mock):
            _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level="qry_cluster",
                top_n=5,
            )

        assert seen_cols == ["ref_type", "qry_cluster"]

    def test_joint_overlapping_genes_not_duplicated(self, ref_adata, qry_adata):
        from spida.I.allcools import _allcools_cef

        # both ref and qry mark gene_0 as enriched; union should not duplicate it
        mock = _make_sequential_mock(ref_indices=[0, 1], qry_indices=[0, 2])
        with patch("ALLCools.clustering.cluster_enriched_features", mock):
            result = _allcools_cef(
                ref_adata.copy(), qry_adata.copy(),
                ref_col_level="ref_type",
                qry_col_level="qry_cluster",
                top_n=2,
            )

        assert len(result) == len(set(result))
        assert set(result) == {f"gene_{i}" for i in [0, 1, 2]}


# ---------------------------------------------------------------------------
# Fixtures for c2c transfer tests
# ---------------------------------------------------------------------------

def _make_adata_comb():
    """Combined AnnData with 4 integrated Leiden clusters.

    cluster_0: 8 ref TypeA + 2 ref TypeB  →  TypeA at 0.80,  5 qry cells  (passes all defaults)
    cluster_1: 4 ref TypeA + 6 ref TypeB  →  TypeB at 0.60,  5 qry cells  (proportion too low)
    cluster_2: 10 ref TypeA               →  TypeA at 1.00,  2 qry cells  (qry_cell_count too low)
    cluster_3: 3 ref TypeA                →  TypeA at 1.00,  5 qry cells  (ref_cell_count too low)
    """
    rng = np.random.default_rng(42)
    rows = []
    # cluster_0
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_0", "CellType": "TypeA"}] * 8
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_0", "CellType": "TypeB"}] * 2
    rows += [{"Modality": "query", "integrated_leiden": "cluster_0", "CellType": "unknown"}] * 5
    # cluster_1
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_1", "CellType": "TypeA"}] * 4
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_1", "CellType": "TypeB"}] * 6
    rows += [{"Modality": "query", "integrated_leiden": "cluster_1", "CellType": "unknown"}] * 5
    # cluster_2 — qry count too low
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_2", "CellType": "TypeA"}] * 10
    rows += [{"Modality": "query", "integrated_leiden": "cluster_2", "CellType": "unknown"}] * 2
    # cluster_3 — ref count too low
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_3", "CellType": "TypeA"}] * 3
    rows += [{"Modality": "query", "integrated_leiden": "cluster_3", "CellType": "unknown"}] * 5

    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
    n_genes = 5
    X = sp.csr_matrix(np.zeros((len(rows), n_genes)))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["integrated_leiden"] = adata.obs["integrated_leiden"].astype("category")
    return adata


# ---------------------------------------------------------------------------
# _clust2clust_transfer (proportion method)
# ---------------------------------------------------------------------------

class TestClust2ClustTransfer:

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.adata = _make_adata_comb()

    def test_assigns_above_threshold(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.8,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        # cluster_0 passes (0.80 TypeA, 10 ref, 5 qry)
        assert (qry.loc[qry["integrated_leiden"] == "cluster_0", "c2c_allcools_label_CellType"] == "TypeA").all()

    def test_unknown_below_proportion(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.8,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        # cluster_1: TypeB proportion = 0.6, below 0.8
        assert (qry.loc[qry["integrated_leiden"] == "cluster_1", "c2c_allcools_label_CellType"] == "unknown").all()

    def test_unknown_below_qry_cell_count(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.8,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        # cluster_2: only 2 query cells, below qry_cell_count=5
        assert (qry.loc[qry["integrated_leiden"] == "cluster_2", "c2c_allcools_label_CellType"] == "unknown").all()

    def test_unknown_below_ref_cell_count(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.8,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        # cluster_3: only 3 ref cells, below ref_cell_count=5
        assert (qry.loc[qry["integrated_leiden"] == "cluster_3", "c2c_allcools_label_CellType"] == "unknown").all()

    def test_saves_map_to_uns(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
        )
        assert "integrated_leiden_CellType_map" in result.uns
        df_map = result.uns["integrated_leiden_CellType_map"]
        assert set(df_map.index) == {"cluster_0", "cluster_1", "cluster_2", "cluster_3"}

    def test_result_is_categorical(self):
        from spida.I.allcools import _clust2clust_transfer
        result = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
        )
        assert result.obs["c2c_allcools_label_CellType"].dtype.name == "category"

    def test_lower_proportion_threshold_assigns_more(self):
        from spida.I.allcools import _clust2clust_transfer
        strict = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.8,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        lenient = _clust2clust_transfer(
            self.adata.copy(),
            ref_cell_type_column="CellType",
            joint_cluster_column="integrated_leiden",
            ref_proportion=0.5,
            ref_cell_count=5,
            qry_cell_count=5,
        )
        qry_strict = strict.obs.loc[strict.obs["Modality"] == "query", "c2c_allcools_label_CellType"]
        qry_lenient = lenient.obs.loc[lenient.obs["Modality"] == "query", "c2c_allcools_label_CellType"]
        n_known_strict = (qry_strict != "unknown").sum()
        n_known_lenient = (qry_lenient != "unknown").sum()
        assert n_known_lenient >= n_known_strict


# ---------------------------------------------------------------------------
# _allcools_clust2clust_transfer (ALLCools confusion-matrix method)
# ---------------------------------------------------------------------------

def _make_adata_comb_clean():
    """Combined AnnData with 2 cleanly separated Leiden clusters (> threshold cell counts).

    cluster_0: 60 ref TypeA, 30 qry cells (NaN ref-type, as in real data)
    cluster_1: 60 ref TypeB, 30 qry cells
    """
    rows = []
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_0", "CellType": "TypeA"}] * 60
    rows += [{"Modality": "query", "integrated_leiden": "cluster_0", "CellType": np.nan}] * 30
    rows += [{"Modality": "ref",   "integrated_leiden": "cluster_1", "CellType": "TypeB"}] * 60
    rows += [{"Modality": "query", "integrated_leiden": "cluster_1", "CellType": np.nan}] * 30

    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
    obs["integrated_leiden"] = obs["integrated_leiden"].astype("category")
    n_genes = 5
    X = sp.csr_matrix(np.zeros((len(rows), n_genes)))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _mock_cmc_one_to_one(confusion_matrix, **kwargs):
    """Each Leiden cluster maps to its own group (trivially unambiguous)."""
    # confusion_matrix.T was passed; rows=cell_types, cols=Leiden clusters
    # ref_group: Leiden clusters → group id
    # query_group: cell_types → group id
    leiden_clusters = confusion_matrix.columns  # columns of .T = Leiden clusters
    cell_types = confusion_matrix.index         # rows of .T = cell_types
    ref_group = pd.Series(range(len(leiden_clusters)), index=leiden_clusters)
    query_group = pd.Series(range(len(cell_types)), index=cell_types)
    return query_group, ref_group, confusion_matrix, None, 1.0


class TestAllcoolsClust2ClustTransfer:

    @pytest.fixture(autouse=True)
    def _skip_if_no_allcools(self):
        pytest.importorskip("ALLCools")

    def test_assigns_unambiguous_clusters(self):
        from spida.I.allcools import _allcools_clust2clust_transfer
        adata = _make_adata_comb_clean()
        with patch("spida.I.allcools.confusion_matrix_clustering", _mock_cmc_one_to_one):
            result = _allcools_clust2clust_transfer(
                adata,
                ref_cell_type_column="CellType",
                joint_cluster_column="integrated_leiden",
                qry_only_cluster_threshold=5,
                ref_only_cluster_threshold=5,
            )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        assert (qry.loc[qry["integrated_leiden"] == "cluster_0", "c2c_allcools_label_CellType"] == "TypeA").all()
        assert (qry.loc[qry["integrated_leiden"] == "cluster_1", "c2c_allcools_label_CellType"] == "TypeB").all()

    def test_small_clusters_excluded(self):
        """Clusters below qry_only_cluster_threshold are dropped before the confusion
        matrix and their query cells remain "unknown".

        Note: data is built with value_counts(dropna=True), so only ref cells (which
        have a real CellType) contribute to the count.  qry_only_cluster_threshold
        therefore filters on *ref cell count* per Leiden cluster.
        """
        from spida.I.allcools import _allcools_clust2clust_transfer
        rows = []
        # cluster_0: 60 ref TypeA → passes threshold
        rows += [{"Modality": "ref",   "integrated_leiden": "cluster_0", "CellType": "TypeA"}] * 60
        rows += [{"Modality": "query", "integrated_leiden": "cluster_0", "CellType": np.nan}] * 30
        # cluster_small: only 3 ref cells — below qry_only_cluster_threshold=10
        rows += [{"Modality": "ref",   "integrated_leiden": "cluster_small", "CellType": "TypeA"}] * 3
        rows += [{"Modality": "query", "integrated_leiden": "cluster_small", "CellType": np.nan}] * 30
        obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
        obs["integrated_leiden"] = obs["integrated_leiden"].astype("category")
        X = sp.csr_matrix(np.zeros((len(rows), 5)))
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
        adata = ad.AnnData(X=X, obs=obs, var=var)

        with patch("spida.I.allcools.confusion_matrix_clustering", _mock_cmc_one_to_one):
            result = _allcools_clust2clust_transfer(
                adata,
                ref_cell_type_column="CellType",
                joint_cluster_column="integrated_leiden",
                qry_only_cluster_threshold=10,
                ref_only_cluster_threshold=5,
            )
        qry = result.obs.loc[result.obs["Modality"] == "query"]
        # cluster_small is dropped (3 ref cells < 10), its query cells stay "unknown"
        assert (qry.loc[qry["integrated_leiden"] == "cluster_small", "c2c_allcools_label_CellType"] == "unknown").all()
        # cluster_0 should be labeled
        assert (qry.loc[qry["integrated_leiden"] == "cluster_0", "c2c_allcools_label_CellType"] == "TypeA").all()

    def test_populates_uns(self):

        from spida.I.allcools import _allcools_clust2clust_transfer
        adata = _make_adata_comb_clean()
        with patch("spida.I.allcools.confusion_matrix_clustering", _mock_cmc_one_to_one):
            result = _allcools_clust2clust_transfer(
                adata,
                ref_cell_type_column="CellType",
                joint_cluster_column="integrated_leiden",
                qry_only_cluster_threshold=5,
                ref_only_cluster_threshold=5,
            )
        assert "c2c_allcools_integration_results" in result.uns
        keys = result.uns["c2c_allcools_integration_results"].keys()
        assert "integration_group_cluster" in keys


# ---------------------------------------------------------------------------
# filt_adata_allcools
# ---------------------------------------------------------------------------

def _make_adata_with_scores():
    """AnnData with allcools_CellType and allcools_CellType_transfer_score columns.

    cell_0–3:  score 0.9  (TypeA)   → kept
    cell_4–5:  score 0.3  (TypeB)   → filtered (below 0.5)
    cell_6:    score 0.5  (TypeA)   → kept (at boundary, not strictly less than)
    cell_7:    NaN score  (TypeA)   → left unchanged (ref cell, no score)
    """
    n = 8
    obs = pd.DataFrame({
        "allcools_CellType": ["TypeA"] * 4 + ["TypeB"] * 2 + ["TypeA"] * 2,
        "allcools_CellType_transfer_score": [0.9, 0.9, 0.9, 0.9, 0.3, 0.3, 0.5, np.nan],
    }, index=[f"cell_{i}" for i in range(n)])
    X = sp.csr_matrix(np.zeros((n, 3)))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(3)])
    return ad.AnnData(X=X, obs=obs, var=var)


class TestFiltAdataAllcools:

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.adata = _make_adata_with_scores()

    def test_below_threshold_marked_unknown(self):
        from spida.I.allcools import filt_adata_allcools
        filt_adata_allcools(self.adata, "allcools_CellType", transfer_score_threshold=0.5)
        assert (self.adata.obs.loc[["cell_4", "cell_5"], "allcools_CellType_filt"] == "unknown").all()

    def test_above_threshold_kept(self):
        from spida.I.allcools import filt_adata_allcools
        filt_adata_allcools(self.adata, "allcools_CellType", transfer_score_threshold=0.5)
        assert (self.adata.obs.loc[["cell_0", "cell_1", "cell_2", "cell_3"], "allcools_CellType_filt"] == "TypeA").all()

    def test_boundary_score_kept(self):
        # score == threshold is NOT below, so should be kept
        from spida.I.allcools import filt_adata_allcools
        filt_adata_allcools(self.adata, "allcools_CellType", transfer_score_threshold=0.5)
        assert self.adata.obs.loc["cell_6", "allcools_CellType_filt"] == "TypeA"

    def test_nan_score_unchanged(self):
        # cell_7 has NaN score — no score to filter on, label should pass through
        from spida.I.allcools import filt_adata_allcools
        filt_adata_allcools(self.adata, "allcools_CellType", transfer_score_threshold=0.5)
        assert self.adata.obs.loc["cell_7", "allcools_CellType_filt"] == "TypeA"

    def test_custom_score_and_filt_col(self):
        from spida.I.allcools import filt_adata_allcools
        self.adata.obs["custom_score"] = [0.1] * 4 + [0.9] * 4
        filt_adata_allcools(
            self.adata, "allcools_CellType",
            transfer_score_threshold=0.5,
            score_col="custom_score",
            filt_col="my_filt",
        )
        assert (self.adata.obs.loc[["cell_0", "cell_1", "cell_2", "cell_3"], "my_filt"] == "unknown").all()
        assert (self.adata.obs.loc[["cell_4", "cell_5", "cell_6"], "my_filt"] != "unknown").all()

    def test_result_is_categorical(self):
        from spida.I.allcools import filt_adata_allcools
        filt_adata_allcools(self.adata, "allcools_CellType")
        assert self.adata.obs["allcools_CellType_filt"].dtype.name == "category"

    def test_higher_threshold_filters_more(self):
        from spida.I.allcools import filt_adata_allcools
        adata_strict = self.adata.copy()
        adata_lenient = self.adata.copy()
        filt_adata_allcools(adata_strict, "allcools_CellType", transfer_score_threshold=0.95)
        filt_adata_allcools(adata_lenient, "allcools_CellType", transfer_score_threshold=0.5)
        n_unknown_strict = (adata_strict.obs["allcools_CellType_filt"] == "unknown").sum()
        n_unknown_lenient = (adata_lenient.obs["allcools_CellType_filt"] == "unknown").sum()
        assert n_unknown_strict >= n_unknown_lenient
