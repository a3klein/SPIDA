"""
Smoke tests for spida.utilities.degs:  call_degs_scanpy, call_degs_scran,
call_degs_pydeseq2, and their summarize_*_from_uns counterparts.

Covers:
  - Basic execution with all non-optional parameters varied
  - uns_key storage (in-memory)
  - h5ad round-trip: write → read → uns key intact
  - zarr round-trip: write → read → uns key intact
  - Equality: summary rebuilt from uns matches summary from direct results
"""

import importlib
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def small_adata():
    """
    120 cells × 30 genes, 3 cell types (40 each), 4 samples (10 per ct×sample).
    TypeA overexpresses genes 0-4, TypeB overexpresses genes 5-9.
    Layers: 'counts' (raw int), 'lognorm' (log-normalised).
    """
    rng = np.random.default_rng(42)
    n_cells, n_genes = 120, 30

    cell_types = np.repeat(['TypeA', 'TypeB', 'TypeC'], 40)
    samples    = np.tile(['S1', 'S2', 'S3', 'S4'], 30)

    counts = rng.negative_binomial(8, 0.5, size=(n_cells, n_genes)).astype(float)
    counts[cell_types == 'TypeA', :5]  = rng.negative_binomial(60, 0.5, size=(40, 5))
    counts[cell_types == 'TypeB', 5:10] = rng.negative_binomial(60, 0.5, size=(40, 5))

    lib_size = counts.sum(axis=1, keepdims=True).clip(1)
    lognorm  = np.log1p(counts / lib_size * 1e4)

    obs = pd.DataFrame(
        {'cell_type': cell_types, 'sample': samples},
        index=[f'cell_{i}' for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])

    return ad.AnnData(
        X=sp.csr_matrix(lognorm),
        obs=obs,
        var=var,
        layers={
            'counts': sp.csr_matrix(counts.astype(int)),
            'lognorm': sp.csr_matrix(lognorm),
        },
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _lists_equal(a, b):
    return sorted(a) == sorted(b)


def _assert_summaries_equal(s1, s2, gene_col):
    s1 = s1.set_index('cell_type').sort_index()
    s2 = s2.set_index('cell_type').sort_index()
    assert list(s1.index) == list(s2.index), "cell_type index mismatch"
    numeric_cols = [c for c in s1.columns if c != gene_col and s1[c].dtype.kind in 'iuf']
    for col in numeric_cols:
        pd.testing.assert_series_equal(
            s1[col].reset_index(drop=True),
            s2[col].reset_index(drop=True),
            check_names=False,
            rtol=1e-5,
        )
    for ct in s1.index:
        assert _lists_equal(s1.loc[ct, gene_col], s2.loc[ct, gene_col]), \
            f"Gene list mismatch for {ct}"


# ---------------------------------------------------------------------------
# call_degs_scanpy
# ---------------------------------------------------------------------------

class TestCallDegsScampy:

    def test_basic(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        results = call_degs_scanpy(small_adata, celltype_col='cell_type')
        assert isinstance(results, dict)
        assert set(results.keys()) == {'TypeA', 'TypeB', 'TypeC'}
        for ct, df in results.items():
            assert 'names' in df.columns
            assert 'logfoldchanges' in df.columns
            assert 'significant' in df.columns

    def test_with_layer(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        results = call_degs_scanpy(
            small_adata, celltype_col='cell_type', layer='lognorm'
        )
        assert isinstance(results, dict)
        assert len(results) > 0

    def test_method_param(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        for method in ('wilcoxon', 't-test', 't-test_overestim_var'):
            results = call_degs_scanpy(
                small_adata, celltype_col='cell_type', method=method
            )
            assert isinstance(results, dict)

    def test_thresholds(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        results = call_degs_scanpy(
            small_adata, celltype_col='cell_type',
            logfc_threshold=0.5, pval_threshold=0.01,
        )
        assert isinstance(results, dict)

    def test_min_max_cells(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        results = call_degs_scanpy(
            small_adata, celltype_col='cell_type',
            min_cells=5, max_cells=100,
        )
        assert isinstance(results, dict)

    def test_uns_key_stored(self, small_adata):
        from spida.utilities.degs import call_degs_scanpy
        adata = small_adata.copy()
        call_degs_scanpy(adata, celltype_col='cell_type', uns_key='scanpy_markers')
        assert 'scanpy_markers' in adata.uns
        stacked = adata.uns['scanpy_markers']
        assert isinstance(stacked, pd.DataFrame)
        assert 'cell_type' in stacked.columns
        assert set(stacked['cell_type'].unique()) == {'TypeA', 'TypeB', 'TypeC'}

    def test_h5ad_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_scanpy
        adata = small_adata.copy()
        call_degs_scanpy(adata, celltype_col='cell_type', uns_key='scanpy_markers')
        h5_path = tmp_path / 'test.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)
        assert 'scanpy_markers' in adata2.uns
        assert isinstance(adata2.uns['scanpy_markers'], pd.DataFrame)
        assert 'cell_type' in adata2.uns['scanpy_markers'].columns

    def test_zarr_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_scanpy
        adata = small_adata.copy()
        call_degs_scanpy(adata, celltype_col='cell_type', uns_key='scanpy_markers')
        zarr_path = tmp_path / 'test.zarr'
        adata.write_zarr(zarr_path)
        adata2 = ad.read_zarr(zarr_path)
        assert 'scanpy_markers' in adata2.uns
        assert isinstance(adata2.uns['scanpy_markers'], pd.DataFrame)
        assert 'cell_type' in adata2.uns['scanpy_markers'].columns

    def test_summarize_from_uns_equals_direct(self, small_adata, tmp_path):
        from spida.utilities.degs import (
            call_degs_scanpy, summarize_scanpy_results,
            summarize_scanpy_results_from_uns,
        )
        adata = small_adata.copy()
        results = call_degs_scanpy(
            adata, celltype_col='cell_type', uns_key='scanpy_markers'
        )
        # Write and read back to test full serialization round-trip
        h5_path = tmp_path / 'scanpy_eq.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)

        s_direct = summarize_scanpy_results(results, top_n=5)
        s_uns    = summarize_scanpy_results_from_uns(adata2, uns_key='scanpy_markers', top_n=5)

        _assert_summaries_equal(s_direct, s_uns, gene_col='top_upregulated')


# ---------------------------------------------------------------------------
# call_degs_scran
# ---------------------------------------------------------------------------

class TestCallDegsScran:

    @pytest.fixture(autouse=True)
    def skip_if_no_scranpy(self):
        pytest.importorskip('scranpy')

    def test_basic(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        results = call_degs_scran(small_adata, celltype_col='cell_type', num_threads=1)
        assert isinstance(results, dict)
        assert set(results.keys()) == {'TypeA', 'TypeB', 'TypeC'}
        for ct, df in results.items():
            assert 'names' in df.columns
            assert 'auc_mean' in df.columns
            assert 'auc_min_rank' in df.columns

    def test_with_layer(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        results = call_degs_scran(
            small_adata, celltype_col='cell_type', layer='lognorm', num_threads=1
        )
        assert isinstance(results, dict)

    def test_with_block(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        results = call_degs_scran(
            small_adata, celltype_col='cell_type',
            block_col='sample', num_threads=1,
        )
        assert isinstance(results, dict)

    def test_min_cells(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        results = call_degs_scran(
            small_adata, celltype_col='cell_type',
            min_cells=5, num_threads=1,
        )
        assert isinstance(results, dict)

    def test_max_cells_per_cluster(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        results = call_degs_scran(
            small_adata, celltype_col='cell_type',
            max_cells_per_cluster=20, num_threads=1,
        )
        assert isinstance(results, dict)

    def test_uns_key_stored(self, small_adata):
        from spida.utilities.degs import call_degs_scran
        adata = small_adata.copy()
        call_degs_scran(adata, celltype_col='cell_type', uns_key='scran_markers', num_threads=1)
        assert 'scran_markers' in adata.uns
        stacked = adata.uns['scran_markers']
        assert isinstance(stacked, pd.DataFrame)
        assert 'cell_type' in stacked.columns
        assert set(stacked['cell_type'].unique()) == {'TypeA', 'TypeB', 'TypeC'}

    def test_h5ad_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_scran
        adata = small_adata.copy()
        call_degs_scran(adata, celltype_col='cell_type', uns_key='scran_markers', num_threads=1)
        h5_path = tmp_path / 'test.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)
        assert 'scran_markers' in adata2.uns
        assert isinstance(adata2.uns['scran_markers'], pd.DataFrame)
        assert 'cell_type' in adata2.uns['scran_markers'].columns

    def test_zarr_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_scran
        adata = small_adata.copy()
        call_degs_scran(adata, celltype_col='cell_type', uns_key='scran_markers', num_threads=1)
        zarr_path = tmp_path / 'test.zarr'
        adata.write_zarr(zarr_path)
        adata2 = ad.read_zarr(zarr_path)
        assert 'scran_markers' in adata2.uns
        assert isinstance(adata2.uns['scran_markers'], pd.DataFrame)
        assert 'cell_type' in adata2.uns['scran_markers'].columns

    def test_summarize_from_uns_equals_direct(self, small_adata, tmp_path):
        from spida.utilities.degs import (
            call_degs_scran, summarize_scran_results,
            summarize_scran_results_from_uns,
        )
        adata = small_adata.copy()
        results = call_degs_scran(
            adata, celltype_col='cell_type', uns_key='scran_markers', num_threads=1
        )
        h5_path = tmp_path / 'scran_eq.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)

        s_direct = summarize_scran_results(results, top_n=5, auc_threshold=0.6)
        s_uns    = summarize_scran_results_from_uns(
            adata2, uns_key='scran_markers', top_n=5, auc_threshold=0.6
        )

        _assert_summaries_equal(s_direct, s_uns, gene_col='top_5_markers')


# ---------------------------------------------------------------------------
# call_degs_pydeseq2
# ---------------------------------------------------------------------------

class TestCallDegsPydeseq2:

    @pytest.fixture(autouse=True)
    def skip_if_no_pydeseq2(self):
        pytest.importorskip('pydeseq2')

    def test_basic(self, small_adata):
        from spida.utilities.degs import call_degs_pydeseq2
        results = call_degs_pydeseq2(
            small_adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
        )
        assert isinstance(results, pd.DataFrame)
        assert 'log2FoldChange' in results.columns
        assert 'padj' in results.columns
        assert 'significant' in results.columns
        assert (results['comparison'] == 'TypeA_vs_TypeB').all()

    def test_no_lfc_shrink(self, small_adata):
        from spida.utilities.degs import call_degs_pydeseq2
        results = call_degs_pydeseq2(
            small_adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            lfc_shrink=False,
        )
        assert isinstance(results, pd.DataFrame)

    def test_thresholds(self, small_adata):
        from spida.utilities.degs import call_degs_pydeseq2
        results = call_degs_pydeseq2(
            small_adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            pval_threshold=0.01,
            logfc_threshold=1.0,
        )
        assert isinstance(results, pd.DataFrame)

    def test_uns_key_stored(self, small_adata):
        from spida.utilities.degs import call_degs_pydeseq2
        adata = small_adata.copy()
        call_degs_pydeseq2(
            adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            uns_key='deseq2_AB',
        )
        assert 'deseq2_AB' in adata.uns
        df = adata.uns['deseq2_AB']
        assert isinstance(df, pd.DataFrame)
        assert 'log2FoldChange' in df.columns

    def test_h5ad_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_pydeseq2
        adata = small_adata.copy()
        results = call_degs_pydeseq2(
            adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            uns_key='deseq2_AB',
        )
        h5_path = tmp_path / 'test.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)
        assert 'deseq2_AB' in adata2.uns
        df2 = adata2.uns['deseq2_AB']
        assert isinstance(df2, pd.DataFrame)
        assert 'log2FoldChange' in df2.columns

    def test_zarr_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_pydeseq2
        adata = small_adata.copy()
        call_degs_pydeseq2(
            adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            uns_key='deseq2_AB',
        )
        zarr_path = tmp_path / 'test.zarr'
        adata.write_zarr(zarr_path)
        adata2 = ad.read_zarr(zarr_path)
        assert 'deseq2_AB' in adata2.uns
        assert 'log2FoldChange' in adata2.uns['deseq2_AB'].columns

    def test_results_equal_after_h5ad_roundtrip(self, small_adata, tmp_path):
        from spida.utilities.degs import call_degs_pydeseq2
        adata = small_adata.copy()
        results = call_degs_pydeseq2(
            adata,
            celltype_col='cell_type',
            group1='TypeA', group2='TypeB',
            sample_col='sample',
            layer='counts',
            uns_key='deseq2_AB',
        )
        h5_path = tmp_path / 'deseq2_eq.h5ad'
        adata.write_h5ad(h5_path)
        adata2 = ad.read_h5ad(h5_path)

        r_from_uns = adata2.uns['deseq2_AB'].copy()
        r_from_uns['significant'] = r_from_uns['significant'].astype(bool)

        shared_genes = results['names'].values
        r1 = results.set_index('names').loc[shared_genes]
        r2 = r_from_uns.set_index('names').loc[shared_genes]

        pd.testing.assert_series_equal(
            r1['log2FoldChange'], r2['log2FoldChange'], check_names=False, rtol=1e-5
        )
        pd.testing.assert_series_equal(
            r1['significant'].astype(bool), r2['significant'].astype(bool),
            check_names=False,
        )
