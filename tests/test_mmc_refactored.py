"""
MMC (MapMyCells): MMCConfig, MMCPreprocessor, MMCAnnotator, label transfer.
"""

import tempfile
import logging
from pathlib import Path

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

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
            "annotations_store": tmpdir / "annotations_store",
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
def mmc_config(temp_dirs):
    """MMCConfig built from a SpidaConfig."""
    from spida.config import SpidaConfig
    from spida.I.mmc import MMCConfig

    cfg = SpidaConfig(mmc_store_path=temp_dirs["mmc_store"], gene_panel_path=temp_dirs["gene_panel"])
    return MMCConfig.from_spida(cfg, n_cpu=2, bootstrap_factor=0.9, bootstrap_iterations=50)


# ============================================================================
# Tests for MMCConfig
# ============================================================================

class TestMMCConfig:
    """Tests for MMCConfig dataclass."""

    def test_config_from_spida(self, mmc_config, temp_dirs):
        assert mmc_config.mmc_store_path == temp_dirs["mmc_store"]
        assert mmc_config.gene_panel_path == temp_dirs["gene_panel"]
        assert mmc_config.n_cpu == 2
        assert mmc_config.bootstrap_factor == 0.9
        assert mmc_config.bootstrap_iterations == 50

    def test_config_from_spida_requires_mmc_store(self):
        from spida.config import SpidaConfig
        from spida.I.mmc import MMCConfig

        with pytest.raises(ValueError):
            MMCConfig.from_spida(SpidaConfig())

    def test_config_validate(self, mmc_config):
        """Test config validation."""
        # Should not raise
        mmc_config.validate()

    def test_missing_gene_panel_is_not_found(self, temp_dirs):
        """No gene panel dir behaves like a missing gene file (mapping lookup falls back to {})."""
        from spida.I.mmc import MMCConfig, MMCPreprocessor

        config = MMCConfig(mmc_store_path=temp_dirs["mmc_store"])
        config.validate()
        assert MMCPreprocessor(config).load_gene_name_mapping("x") == {}


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

    @pytest.mark.parametrize("annot_name", ["mmc", "mmc_ref-A"])
    def test_transfer_labels(self, mmc_config, sample_query_adata, temp_dirs, annot_name):
        """Labels land in {annot_name}_{level} / {annot_name}_{level}_transfer_score."""
        from spida.I.mmc import MMCAnnotator

        n = sample_query_adata.n_obs
        results_csv = temp_dirs["annotations_store"] / "results.csv"
        pd.DataFrame({
            "level1_name": np.resize(["L1_X", "L1_Y"], n),
            "level1_bootstrapping_probability": np.random.rand(n),
            "level2_name": np.resize(["L2_A", "L2_B", "L2_C"], n),
            "level2_bootstrapping_probability": np.random.rand(n),
        }, index=sample_query_adata.obs_names).to_csv(results_csv)

        adata = MMCAnnotator(mmc_config).transfer_labels(sample_query_adata, results_csv, annot_name=annot_name)

        for level in ("level1", "level2"):
            assert adata.obs[f"{annot_name}_{level}"].notna().all()
            assert adata.obs[f"{annot_name}_{level}_transfer_score"].between(0, 1).all()
        assert adata.obs[f"{annot_name}_level1"].iloc[1] == "L1_Y"
