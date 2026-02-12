"""
MapMyCells Refactored Implementation

This module demonstrates the proposed refactoring of the MapMyCells integration.
It includes the new architecture with separated concerns and improved usability.

To integrate this:
1. Add these classes to mmc.py
2. Update main.py to use new functions
3. Keep old functions for backwards compatibility (deprecated)
"""

import os
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple, Dict
import glob

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

logger = logging.getLogger(__name__)


# ============================================================================
# Configuration Management
# ============================================================================

@dataclass
class MMCConfig:
    """Centralized configuration for MapMyCells operations."""
    
    # Core storage paths
    mmc_store_path: Path
    anndata_store_path: Path
    annotations_store_path: Path
    zarr_storage_path: Path
    gene_panel_path: Path
    
    # Optional overrides
    n_cpu: int = 8
    bootstrap_factor: float = 0.8
    bootstrap_iterations: int = 100
    rng_seed: int = 13
    ref_norm: str = "log2CPM"
    query_norm: str = "log2CPM"
    
    def __post_init__(self):
        """Convert string paths to Path objects."""
        if isinstance(self.mmc_store_path, str):
            self.mmc_store_path = Path(self.mmc_store_path)
        if isinstance(self.anndata_store_path, str):
            self.anndata_store_path = Path(self.anndata_store_path)
        if isinstance(self.annotations_store_path, str):
            self.annotations_store_path = Path(self.annotations_store_path)
        if isinstance(self.zarr_storage_path, str):
            self.zarr_storage_path = Path(self.zarr_storage_path)
        if isinstance(self.gene_panel_path, str):
            self.gene_panel_path = Path(self.gene_panel_path)
    
    def get_mmc_store_dir(self, brain_region: str, codebook: str) -> Path:
        """Get or create the MMC store directory for a given region/codebook."""
        store_dir = self.mmc_store_path / f"{brain_region}-{codebook}"
        store_dir.mkdir(parents=True, exist_ok=True)
        return store_dir
    
    def get_marker_dir(self, brain_region: str, codebook: str) -> Path:
        """Get marker directory."""
        marker_dir = self.get_mmc_store_dir(brain_region, codebook) / "marker_dir"
        marker_dir.mkdir(parents=True, exist_ok=True)
        return marker_dir
    
    def get_precomp_stats_path(
        self, brain_region: str, codebook: str, ref_name: str
    ) -> Path:
        """Get precomputed stats file path."""
        store_dir = self.get_mmc_store_dir(brain_region, codebook)
        return store_dir / f"{ref_name}_precomputed_stats.h5"
    
    def get_query_path(self, brain_region: str, codebook: str) -> Path:
        """Get query AnnData path."""
        store_dir = self.get_mmc_store_dir(brain_region, codebook)
        return store_dir / "query.h5ad"
    
    def get_codebook_path(self, codebook: str) -> Path:
        """Find codebook file by name."""
        matches = glob.glob(str(self.gene_panel_path / f"*{codebook}*.csv"))
        if not matches:
            raise FileNotFoundError(
                f"No codebook found matching '{codebook}' in {self.gene_panel_path}"
            )
        return Path(matches[0])
    
    @classmethod
    def from_env(cls) -> "MMCConfig":
        """Load configuration from environment variables."""
        required_vars = {
            "MMC_DIR": "mmc_store_path",
            "ANNDATA_STORE_PATH": "anndata_store_path",
            "ANNOTATIONS_STORE_PATH": "annotations_store_path",
            "ZARR_STORAGE_PATH": "zarr_storage_path",
            "GENE_PANEL_PATH": "gene_panel_path",
        }
        
        config_dict = {}
        for env_var, param_name in required_vars.items():
            value = os.getenv(env_var)
            if not value:
                raise ValueError(
                    f"Required environment variable '{env_var}' not set. "
                    f"Set it or provide MMCConfig explicitly."
                )
            config_dict[param_name] = Path(value)
        
        # Optional environment variables with defaults
        config_dict["n_cpu"] = int(os.getenv("MMC_N_CPU", "8"))
        config_dict["bootstrap_factor"] = float(os.getenv("MMC_BOOTSTRAP_FACTOR", "0.8"))
        config_dict["bootstrap_iterations"] = int(os.getenv("MMC_BOOTSTRAP_ITERATIONS", "100"))
        config_dict["rng_seed"] = int(os.getenv("MMC_RNG_SEED", "13"))
        config_dict["ref_norm"] = os.getenv("MMC_REF_NORM", "log2CPM")
        config_dict["query_norm"] = os.getenv("MMC_QUERY_NORM", "log2CPM")
        
        return cls(**config_dict)
    
    def validate(self) -> None:
        """Validate that all paths exist."""
        paths = {
            "mmc_store_path": self.mmc_store_path,
            "anndata_store_path": self.anndata_store_path,
            "annotations_store_path": self.annotations_store_path,
            "zarr_storage_path": self.zarr_storage_path,
            "gene_panel_path": self.gene_panel_path,
        }
        
        for path_name, path_obj in paths.items():
            if not path_obj.exists():
                raise FileNotFoundError(f"{path_name} does not exist: {path_obj}")


# ============================================================================
# Preprocessing Component
# ============================================================================

class MMCPreprocessor:
    """Handles all reference preprocessing for MapMyCells annotation."""
    
    def __init__(self, config: MMCConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def create_empty_query_adata(
        self,
        ref_path: str | Path,
        gene_list: List[str],
        output_path: str | Path,
    ) -> None:
        """
        Create an empty AnnData with reference genes intersected with query genes.
        
        This is used by cell_type_mapper to identify which genes to use.
        """
        self.logger.info(f"Creating empty query AnnData at {output_path}")
        
        ref_adata = ad.read_h5ad(ref_path, backed="r")
        var_index = ref_adata.var.index
        
        # Intersect with available genes
        query_index = var_index.intersection(pd.Index(gene_list))
        self.logger.info(f"Gene intersection: {len(query_index)} / {len(gene_list)} genes found")
        
        # Create minimal AnnData
        query = ad.AnnData(
            X=sparse.csr_matrix(np.zeros((1, len(query_index)))),
            var=pd.DataFrame(index=query_index),
        )
        
        query.write_h5ad(output_path)
        self.logger.info(f"Created query AnnData with {len(query_index)} genes")
    
    def precompute_reference_stats(
        self,
        ref_path: str | Path,
        hierarchy_list: List[str],
        brain_region: str,
        codebook: str,
    ) -> Path:
        """Precompute statistics for the reference AnnData."""
        from cell_type_mapper.cli.precompute_stats_scrattch import (  # type: ignore
            PrecomputationScrattchRunner,
        )
        
        ref_name = Path(ref_path).stem
        output_path = self.config.get_precomp_stats_path(brain_region, codebook, ref_name)
        
        self.logger.info(f"Precomputing reference stats -> {output_path}")
        
        precomputation_config = {
            "hierarchy": hierarchy_list,
            "h5ad_path": str(ref_path),
            "output_path": str(output_path),
            "n_processors": self.config.n_cpu,
            "normalization": self.config.ref_norm,
            "clobber": True,
        }
        
        runner = PrecomputationScrattchRunner(args=[], input_data=precomputation_config)
        runner.run()
        
        self.logger.info("Precomputation complete")
        return output_path
    
    def compute_reference_markers(
        self,
        precomp_path: str | Path,
        brain_region: str,
        codebook: str,
        query_path: str | Path = None,
        n_valid: int = 10,
    ) -> Path:
        """Compute reference markers for the hierarchy."""
        from cell_type_mapper.cli.reference_markers import (  # type: ignore
            ReferenceMarkerRunner,
        )
        
        marker_dir = self.config.get_marker_dir(brain_region, codebook)
        
        self.logger.info(f"Computing reference markers -> {marker_dir}")
        
        reference_config = {
            "precomputed_path_list": [str(precomp_path)],
            "n_valid": n_valid,
            "output_dir": str(marker_dir),
            "clobber": True,
            "query_path": str(query_path) if query_path else None,
        }
        
        runner = ReferenceMarkerRunner(args=[], input_data=reference_config)
        runner.run()
        
        markers_path = marker_dir / "reference_markers.h5"
        self.logger.info(f"Reference markers saved to {markers_path}")
        return markers_path
    
    def compute_query_markers(
        self,
        ref_markers_path: str | Path,
        query_path: str | Path,
        brain_region: str,
        codebook: str,
        n_per_utility: int = 10,
    ) -> Path:
        """Downsample markers to those present in query."""
        from cell_type_mapper.cli.query_markers import (  # type: ignore
            QueryMarkerRunner,
        )
        
        marker_dir = self.config.get_marker_dir(brain_region, codebook)
        output_path = marker_dir / "calc_markers.json"
        
        self.logger.info(f"Computing query markers -> {output_path}")
        
        query_config = {
            "output_path": str(output_path),
            "reference_marker_path_list": [str(ref_markers_path)],
            "n_per_utility": n_per_utility,
            "n_processors": self.config.n_cpu,
            "query_path": str(query_path),
        }
        
        runner = QueryMarkerRunner(args=[], input_data=query_config)
        runner.run()
        
        self.logger.info(f"Query markers saved to {output_path}")
        return output_path
    
    def setup_annotation_pipeline(
        self,
        ref_path: str | Path,
        brain_region: str,
        codebook: str,
        hierarchy_list: List[str],
        codebook_path: str | Path = None,
        **preprocessing_params,
    ) -> Dict[str, Path]:
        """
        Complete setup: creates query AnnData, precomputes stats, and computes markers.
        
        Returns dict with paths to all intermediate files.
        """
        self.logger.info("Starting annotation pipeline setup")
        self.logger.info(f"  Brain Region: {brain_region}")
        self.logger.info(f"  Codebook: {codebook}")
        self.logger.info(f"  Hierarchy: {hierarchy_list}")
        
        # Get codebook path if not provided
        if not codebook_path:
            codebook_path = self.config.get_codebook_path(codebook)
        
        # Read gene list from codebook
        gene_list = pd.read_csv(codebook_path)["name"]
        gene_list = gene_list[~gene_list.str.startswith("Blank-")]
        
        # Create query AnnData
        query_path = self.config.get_query_path(brain_region, codebook)
        self.create_empty_query_adata(ref_path, gene_list, query_path)
        
        # Precompute reference stats
        precomp_stats_path = self.precompute_reference_stats(
            ref_path, hierarchy_list, brain_region, codebook
        )
        
        # Compute reference markers
        n_valid = preprocessing_params.get("n_valid", 10)
        ref_markers_path = self.compute_reference_markers(
            precomp_stats_path, brain_region, codebook, query_path, n_valid
        )
        
        # Compute query markers
        n_per_utility = preprocessing_params.get("n_per_utility", 10)
        query_markers_path = self.compute_query_markers(
            ref_markers_path, query_path, brain_region, codebook, n_per_utility
        )
        
        result = {
            "precomp_stats_path": precomp_stats_path,
            "ref_markers_path": ref_markers_path,
            "query_markers_path": query_markers_path,
            "query_adata_path": query_path,
        }
        
        self.logger.info("Annotation pipeline setup complete")
        return result


# ============================================================================
# Annotation Component
# ============================================================================

class MMCAnnotator:
    """Handles query annotation with MapMyCells."""
    
    def __init__(self, config: MMCConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def annotate(
        self,
        query_adata_path: str | Path,
        brain_region: str,
        codebook: str,
        output_dir: str | Path,
        **annotation_params,
    ) -> Tuple[Path, Path]:
        """Run MapMyCells annotation on query data."""
        from cell_type_mapper.cli.from_specified_markers import (  # type: ignore
            FromSpecifiedMarkersRunner,
        )
        
        # Get paths
        marker_dir = self.config.get_marker_dir(brain_region, codebook)
        markers_path = marker_dir / "calc_markers.json"
        ref_name = Path(query_adata_path).stem
        precomp_path = self.config.get_precomp_stats_path(brain_region, codebook, ref_name)
        
        # Output paths
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        extended_results_path = output_dir / "extended_results.json"
        csv_results_path = output_dir / "csv_results.csv"
        
        self.logger.info(f"Running annotation for {query_adata_path}")
        self.logger.info(f"  Output: {output_dir}")
        
        # Annotation parameters
        bootstrap_factor = annotation_params.get("bootstrap_factor", self.config.bootstrap_factor)
        bootstrap_iterations = annotation_params.get("bootstrap_iterations", self.config.bootstrap_iterations)
        rng_seed = annotation_params.get("rng_seed", self.config.rng_seed)
        n_cpu = annotation_params.get("n_cpu", self.config.n_cpu)
        
        config = {
            "precomputed_stats": {"path": str(precomp_path)},
            "query_markers": {"serialized_lookup": str(markers_path)},
            "type_assignment": {
                "bootstrap_factor": bootstrap_factor,
                "bootstrap_iteration": bootstrap_iterations,
                "rng_seed": rng_seed,
                "n_processors": n_cpu,
                "normalization": self.config.query_norm,
            },
            "query_path": str(query_adata_path),
            "extended_result_path": str(extended_results_path),
            "csv_result_path": str(csv_results_path),
        }
        
        runner = FromSpecifiedMarkersRunner(args=[], input_data=config)
        runner.run()
        
        self.logger.info(f"Annotation complete. Results at {output_dir}")
        return extended_results_path, csv_results_path
    
    def transfer_labels(
        self,
        adata: ad.AnnData,
        csv_results_path: str | Path,
    ) -> ad.AnnData:
        """Transfer annotation labels to AnnData object."""
        self.logger.info(f"Transferring labels from {csv_results_path}")
        
        results_df = pd.read_csv(csv_results_path, comment="#", index_col=0)
        
        # Extract hierarchy levels from column names
        hierarchy_list = []
        for col in results_df.columns:
            if col.endswith("_bootstrapping_probability"):
                hierarchy_list.append(col.split("_bootstrapping_probability")[0])
        
        self.logger.info(f"Found hierarchy levels: {hierarchy_list}")
        
        # Add columns to obs
        for h in hierarchy_list:
            score_col = f"{h}_bootstrapping_probability"
            label_col = f"{h}_name"
            
            adata.obs[f"mmc_{h}_transfer_score"] = results_df[score_col]
            adata.obs[f"mmc_{h}"] = results_df[label_col]
        
        # Log what was added
        added_cols = [f"mmc_{h}" for h in hierarchy_list] + [f"mmc_{h}_transfer_score" for h in hierarchy_list]
        self.logger.info(f"Added columns to adata.obs: {added_cols}")
        
        return adata


# ============================================================================
# High-Level Public API
# ============================================================================

def setup_annotation_pipeline(
    ref_path: str | Path,
    brain_region: str,
    codebook: str,
    hierarchy_list: List[str],
    codebook_path: str | Path = None,
    config: MMCConfig = None,
    **preprocessing_params,
) -> Dict[str, Path]:
    """
    Complete setup for MapMyCells annotation.
    
    This creates all necessary reference data and markers for annotation.
    Should be run once per reference/codebook combination.
    
    Parameters:
    -----------
    ref_path : str or Path
        Path to reference single-cell RNA-seq AnnData
    brain_region : str
        Brain region identifier
    codebook : str
        Codebook identifier for matching
    hierarchy_list : List[str]
        Hierarchy levels for annotation
    codebook_path : str or Path, optional
        Explicit path to codebook CSV
    config : MMCConfig, optional
        Configuration object. If None, loads from environment
    **preprocessing_params
        Additional params (n_cpu, n_valid, n_per_utility, etc.)
    
    Returns:
    --------
    dict with paths to precomputed files
    """
    if config is None:
        config = MMCConfig.from_env()
    
    config.validate()
    
    preprocessor = MMCPreprocessor(config)
    return preprocessor.setup_annotation_pipeline(
        ref_path, brain_region, codebook, hierarchy_list, codebook_path, **preprocessing_params
    )
