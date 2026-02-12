"""
Comprehensive test suite for MMC (MapMyCells) refactored implementation.

Tests cover:
1. Configuration management (MMCConfig)
2. Preprocessing pipeline (MMCPreprocessor)
3. Annotation pipeline (MMCAnnotator)
4. Label transfer functionality
5. CLI wrappers
6. Backwards compatibility with legacy functions
"""

import os
import tempfile
import logging
from pathlib import Path
from typing import Optional

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def temp_dirs(): # Creating the test directories
    """Create temporary directories for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        dirs = {
            "mmc_store": tmpdir / "mmc_store",
            "anndata_store": tmpdir / "anndata_store",
            "annotations_store": tmpdir / "annotations_store",
            "zarr_storage": tmpdir / "zarr_storage",
            "gene_panel": tmpdir / "gene_panel",
        }
        for d in dirs.values():
            d.mkdir(parents=True, exist_ok=True)
        
        yield dirs


@pytest.fixture
def sample_ref_adata():
    """Create a sample reference AnnData object."""
    n_obs = 100
    n_vars = 50
    
    X = sparse.csr_matrix(np.random.poisson(5, (n_obs, n_vars)))
    var_names = [f"Gene_{i}" for i in range(n_vars)]
    obs_names = [f"Cell_{i}" for i in range(n_obs)]
    
    adata = ad.AnnData(
        X=X,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame({
            "cell_type": np.random.choice(["Type_A", "Type_B", "Type_C"], n_obs),
            "level1": np.random.choice(["L1_X", "L1_Y"], n_obs),
            "level2": np.random.choice(["L2_A", "L2_B", "L2_C"], n_obs),
        }, index=obs_names)
    )
    
    return adata


@pytest.fixture
def sample_query_adata():
    """Create a sample query AnnData object."""
    n_obs = 50
    n_vars = 40
    
    X = sparse.csr_matrix(np.random.poisson(5, (n_obs, n_vars)))
    var_names = [f"Gene_{i}" for i in range(40)]  # Overlaps with ref genes
    obs_names = [f"Query_Cell_{i}" for i in range(n_obs)]
    
    adata = ad.AnnData(
        X=X,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame({
            "experiment": "exp1",
            "segmentation": "seg1",
            "donor": "donor1",
        }, index=obs_names)
    )
    
    return adata


@pytest.fixture
def sample_codebook(temp_dirs):
    """Create a sample codebook CSV file."""
    codebook_path = temp_dirs["gene_panel"] / "codebook_test.csv"
    codebook_df = pd.DataFrame({
        "name": [f"Gene_{i}" for i in range(45)] + ["Blank-1", "Blank-2", "Blank-3"],
        "id": range(48),
    })
    codebook_df.to_csv(codebook_path, index=False)
    return codebook_path


@pytest.fixture
def mmc_config(temp_dirs, monkeypatch):
    """Create MMCConfig with environment variables set."""
    from spida.I.mmc import MMCConfig
    
    # Set environment variables
    monkeypatch.setenv("MMC_DIR", str(temp_dirs["mmc_store"]))
    monkeypatch.setenv("ANNDATA_STORE_PATH", str(temp_dirs["anndata_store"]))
    monkeypatch.setenv("ANNOTATIONS_STORE_PATH", str(temp_dirs["annotations_store"]))
    monkeypatch.setenv("ZARR_STORAGE_PATH", str(temp_dirs["zarr_storage"]))
    monkeypatch.setenv("GENE_PANEL_PATH", str(temp_dirs["gene_panel"]))
    monkeypatch.setenv("MMC_N_CPU", "2")
    monkeypatch.setenv("MMC_BOOTSTRAP_FACTOR", "0.9")
    monkeypatch.setenv("MMC_BOOTSTRAP_ITERATIONS", "50")
    
    config = MMCConfig.from_env()
    return config


# ============================================================================
# Tests for MMCConfig
# ============================================================================

class TestMMCConfig:
    """Tests for MMCConfig dataclass."""
    
    def test_config_from_env(self, mmc_config):
        """Test loading config from environment variables."""
        assert mmc_config.mmc_store_path.exists()
        assert mmc_config.anndata_store_path.exists()
        assert mmc_config.annotations_store_path.exists()
        assert mmc_config.n_cpu == 2
        assert mmc_config.bootstrap_factor == 0.9
        assert mmc_config.bootstrap_iterations == 50
    
    def test_config_validate(self, mmc_config):
        """Test config validation."""
        # Should not raise
        mmc_config.validate()
    
    def test_get_mmc_store_dir(self, mmc_config):
        """Test getting/creating MMC store directory."""
        store_dir = mmc_config.get_mmc_store_dir("brain_region_1", "codebook_1")
        expected_dir = mmc_config.mmc_store_path / "brain_region_1-codebook_1"
        assert store_dir == expected_dir
        assert store_dir.exists()
    
    def test_get_marker_dir(self, mmc_config):
        """Test getting marker directory."""
        marker_dir = mmc_config.get_marker_dir("brain_region_1", "codebook_1")
        expected_dir = mmc_config.mmc_store_path / "brain_region_1-codebook_1" / "marker_dir"
        assert marker_dir == expected_dir
        assert marker_dir.exists()
    
    def test_get_precomp_stats_path(self, mmc_config):
        """Test getting precomputed stats path."""
        stats_path = mmc_config.get_precomp_stats_path("region_1", "codebook_1", "ref_data")
        assert stats_path.suffix == ".h5"
        assert "precomputed_stats" in stats_path.name
    
    def test_get_query_path(self, mmc_config):
        """Test getting query AnnData path."""
        query_path = mmc_config.get_query_path("region_1", "codebook_1")
        assert query_path.suffix == ".h5ad"
        assert query_path.name == "query.h5ad"
    
    def test_get_codebook_path(self, mmc_config, sample_codebook):
        """Test finding codebook file."""
        codebook_path = mmc_config.get_codebook_path("test")
        assert codebook_path.exists()
        assert "test" in codebook_path.name


# ============================================================================
# Tests for MMCPreprocessor
# ============================================================================

class TestMMCPreprocessor:
    """Tests for MMCPreprocessor component."""
    
    def test_create_empty_query_adata(self, mmc_config, sample_ref_adata, temp_dirs):
        """Test creating empty query AnnData."""
        from spida.I.mmc import MMCPreprocessor
        
        preprocessor = MMCPreprocessor(mmc_config)
        
        # Save reference
        ref_path = temp_dirs["gene_panel"] / "ref.h5ad"
        sample_ref_adata.write_h5ad(ref_path)
        
        # Create query
        gene_list = list(sample_ref_adata.var_names[:30])
        query_path = temp_dirs["mmc_store"] / "query_test.h5ad"
        
        preprocessor.create_empty_query_adata(ref_path, gene_list, query_path)
        
        # Verify
        assert query_path.exists()
        query_adata = ad.read_h5ad(query_path)
        assert query_adata.n_obs == 1
        assert query_adata.n_vars <= len(gene_list)
    
    def test_preprocessor_initialization(self, mmc_config):
        """Test MMCPreprocessor initialization."""
        from spida.I.mmc import MMCPreprocessor
        
        preprocessor = MMCPreprocessor(mmc_config)
        assert preprocessor.config == mmc_config
        assert preprocessor.logger is not None


# ============================================================================
# Tests for MMCAnnotator
# ============================================================================

class TestMMCAnnotator:
    """Tests for MMCAnnotator component."""
    
    def test_annotator_initialization(self, mmc_config):
        """Test MMCAnnotator initialization."""
        from spida.I.mmc import MMCAnnotator
        
        annotator = MMCAnnotator(mmc_config)
        assert annotator.config == mmc_config
        assert annotator.logger is not None
    
    def test_transfer_labels(self, mmc_config, sample_query_adata, temp_dirs):
        """Test label transfer functionality."""
        from spida.I.mmc import MMCAnnotator
        
        annotator = MMCAnnotator(mmc_config)
        
        # Create dummy results CSV
        results_csv = temp_dirs["annotations_store"] / "results.csv"
        results_df = pd.DataFrame({
            "level1_name": ["L1_X", "L1_Y"] * (sample_query_adata.n_obs // 2 + 1),
            "level1_bootstrapping_probability": np.random.rand(sample_query_adata.n_obs),
            "level2_name": ["L2_A", "L2_B", "L2_C"] * (sample_query_adata.n_obs // 3 + 1),
            "level2_bootstrapping_probability": np.random.rand(sample_query_adata.n_obs),
        })
        results_df = results_df.iloc[:sample_query_adata.n_obs]
        results_df.index = sample_query_adata.obs_names
        results_df.to_csv(results_csv, comment="#")
        
        # Transfer labels
        adata = annotator.transfer_labels(sample_query_adata, results_csv)
        
        # Verify
        assert "mmc_level1" in adata.obs.columns
        assert "mmc_level1_transfer_score" in adata.obs.columns
        assert "mmc_level2" in adata.obs.columns
        assert "mmc_level2_transfer_score" in adata.obs.columns
        assert len(adata.obs["mmc_level1"]) > 0


# ============================================================================
# Tests for Public API Functions
# ============================================================================

class TestPublicAPI:
    """Tests for high-level public API functions."""
    
    def test_setup_annotation_pipeline_validation(self, mmc_config, sample_ref_adata, temp_dirs):
        """Test setup_annotation_pipeline with basic validation."""
        from spida.I.mmc import setup_annotation_pipeline
        
        # Save reference
        ref_path = temp_dirs["gene_panel"] / "ref.h5ad"
        sample_ref_adata.write_h5ad(ref_path)
        
        # This should succeed without error (up to actual cell_type_mapper availability)
        try:
            result = setup_annotation_pipeline(
                ref_path=str(ref_path),
                brain_region="test_region",
                codebook="test_codebook",
                hierarchy_list=["level1", "level2"],
                config=mmc_config,
                n_valid=5,
                n_per_utility=5,
            )
            # If we get here, the function structure is sound
            assert isinstance(result, dict)
        except ImportError:
            # Expected if cell_type_mapper is not installed
            pytest.skip("cell_type_mapper not installed")
    
    def test_annotate_region_parameters(self, mmc_config):
        """Test annotate_region parameter handling."""
        from spida.I.mmc import annotate_region
        
        # We're testing parameter passing, not actual execution
        # so we mock what we can
        try:
            # This will fail at zarr read but validates parameter structure
            annotate_region(
                exp_name="test_exp",
                reg_name="test_region",
                prefix_name="test_prefix",
                brain_region="test_brain_region",
                codebook="test_codebook",
                config=mmc_config,
            )
        except FileNotFoundError:
            # Expected - zarr path doesn't exist
            pass


# ============================================================================
# Tests for Legacy Backward Compatibility
# ============================================================================

class TestBackwardCompatibility:
    """Tests for legacy function backward compatibility."""
    
    def test_setup_mmc_legacy(self, mmc_config, sample_ref_adata, temp_dirs, monkeypatch):
        """Test that legacy setup_mmc still works."""
        from spida.I.mmc import setup_mmc
        
        # Set environment
        monkeypatch.setenv("MMC_DIR", str(mmc_config.mmc_store_path))
        monkeypatch.setenv("ANNDATA_STORE_PATH", str(mmc_config.anndata_store_path))
        monkeypatch.setenv("ANNOTATIONS_STORE_PATH", str(mmc_config.annotations_store_path))
        monkeypatch.setenv("ZARR_STORAGE_PATH", str(mmc_config.zarr_storage_path))
        monkeypatch.setenv("GENE_PANEL_PATH", str(mmc_config.gene_panel_path))
        
        # Save reference
        ref_path = temp_dirs["gene_panel"] / "ref.h5ad"
        sample_ref_adata.write_h5ad(ref_path)
        
        # Create codebook
        codebook_path = temp_dirs["gene_panel"] / "test_codebook.csv"
        pd.DataFrame({"name": list(sample_ref_adata.var_names)}).to_csv(codebook_path, index=False)
        
        try:
            result = setup_mmc(
                ref_path=str(ref_path),
                brain_region="test_region",
                codebook="test",
                hierarchy_list=["level1", "level2"],
                codebook_path=str(codebook_path),
                mmc_store_path=str(mmc_config.mmc_store_path),
            )
            assert result == 0
        except ImportError:
            pytest.skip("cell_type_mapper not installed")


# ============================================================================
# Tests for CLI Wrappers
# ============================================================================

class TestCLIWrappers:
    """Tests for CLI wrapper functions."""
    
    def test_mmc_setup_cli_parameter_handling(self, mmc_config, sample_ref_adata, temp_dirs, monkeypatch):
        """Test mmc_setup_cli parameter handling."""
        from spida.I.mmc_cli import mmc_setup_cli
        
        # Set environment
        monkeypatch.setenv("MMC_DIR", str(mmc_config.mmc_store_path))
        monkeypatch.setenv("ANNDATA_STORE_PATH", str(mmc_config.anndata_store_path))
        monkeypatch.setenv("ANNOTATIONS_STORE_PATH", str(mmc_config.annotations_store_path))
        monkeypatch.setenv("ZARR_STORAGE_PATH", str(mmc_config.zarr_storage_path))
        monkeypatch.setenv("GENE_PANEL_PATH", str(mmc_config.gene_panel_path))
        
        # Save reference
        ref_path = temp_dirs["gene_panel"] / "ref.h5ad"
        sample_ref_adata.write_h5ad(ref_path)
        
        # Create codebook
        codebook_path = temp_dirs["gene_panel"] / "test_codebook.csv"
        pd.DataFrame({"name": list(sample_ref_adata.var_names)}).to_csv(codebook_path, index=False)
        
        try:
            result = mmc_setup_cli(
                ref_path=str(ref_path),
                brain_region="test_region",
                codebook="test",
                hierarchy_list=["level1", "level2"],
                codebook_path=str(codebook_path),
                n_cpu=2,
                n_valid=5,
                n_per_utility=5,
            )
            assert isinstance(result, dict) or result is None
        except ImportError:
            pytest.skip("cell_type_mapper not installed")
    
    def test_mmc_annotation_region_cli_parameter_handling(self, mmc_config, monkeypatch):
        """Test mmc_annotation_region_cli parameter handling."""
        from spida.I.mmc_cli import mmc_annotation_region_cli
        
        # Set environment
        monkeypatch.setenv("MMC_DIR", str(mmc_config.mmc_store_path))
        monkeypatch.setenv("ANNDATA_STORE_PATH", str(mmc_config.anndata_store_path))
        monkeypatch.setenv("ANNOTATIONS_STORE_PATH", str(mmc_config.annotations_store_path))
        monkeypatch.setenv("ZARR_STORAGE_PATH", str(mmc_config.zarr_storage_path))
        monkeypatch.setenv("GENE_PANEL_PATH", str(mmc_config.gene_panel_path))
        
        try:
            result = mmc_annotation_region_cli(
                exp_name="exp1",
                reg_name="region1",
                prefix_name="prefix1",
                brain_region="brain_region1",
                codebook="codebook1",
                n_cpu=2,
                bootstrap_factor=0.8,
                bootstrap_iterations=100,
                rng_seed=13,
            )
        except (FileNotFoundError, KeyError):
            # Expected - zarr paths and data don't exist for this test
            pass


# ============================================================================
# Integration Tests
# ============================================================================

class TestIntegration:
    """Integration tests combining multiple components."""
    
    def test_config_preprocessor_flow(self, mmc_config, sample_ref_adata, temp_dirs):
        """Test configuration and preprocessing flow together."""
        from spida.I.mmc import MMCPreprocessor
        
        # Verify config is properly set up
        assert mmc_config.mmc_store_path.exists()
        
        # Initialize preprocessor
        preprocessor = MMCPreprocessor(mmc_config)
        assert preprocessor.config == mmc_config
        
        # Test helper methods
        store_dir = mmc_config.get_mmc_store_dir("region", "codebook")
        assert store_dir.exists()
        
        marker_dir = mmc_config.get_marker_dir("region", "codebook")
        assert marker_dir.exists()
    
    def test_config_annotator_flow(self, mmc_config, sample_query_adata, temp_dirs):
        """Test configuration and annotation flow together."""
        from spida.I.mmc import MMCAnnotator
        
        # Verify config
        assert mmc_config.annotations_store_path.exists()
        
        # Initialize annotator
        annotator = MMCAnnotator(mmc_config)
        assert annotator.config == mmc_config
        
        # Test label transfer with dummy data
        results_csv = temp_dirs["annotations_store"] / "results.csv"
        results_df = pd.DataFrame({
            "cell_type_name": ["Type_A", "Type_B"] * (sample_query_adata.n_obs // 2 + 1),
            "cell_type_bootstrapping_probability": np.random.rand(sample_query_adata.n_obs),
        })
        results_df = results_df.iloc[:sample_query_adata.n_obs]
        results_df.index = sample_query_adata.obs_names
        results_df.to_csv(results_csv, comment="#")
        
        adata = annotator.transfer_labels(sample_query_adata, results_csv)
        assert "mmc_cell_type" in adata.obs.columns


# ============================================================================
# Main Test Execution
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
