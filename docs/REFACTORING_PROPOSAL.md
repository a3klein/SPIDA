# MapMyCells (mmc.py) Refactoring Proposal

## Overview
This document outlines a comprehensive refactoring of the MapMyCells integration module to improve code organization, maintainability, and extensibility.

---

## Current Issues Summary

### Structural Problems
1. **Mixed responsibilities** - Single functions handle data I/O, configuration, validation, and execution
2. **Hardcoded file paths** - Directory structures embedded in multiple functions make adaptation difficult
3. **Scattered environment variable handling** - No centralized configuration management
4. **Limited error handling** - Silent failures when expected files/directories don't exist
5. **No metadata tracking** - Cannot audit when/how annotations were run

### Code Quality Issues
1. **Kwargs mishandling** - Fixed in immediate improvements
2. **Parameter naming inconsistency** - Fixed (hierarchy_list standardization)
3. **Redundant wrapper functions** - Fixed (mmc_setup alias)

---

## Proposed Architecture

### Phase 1: Configuration Management

#### Create `MMCConfig` dataclass
```python
from dataclasses import dataclass
from pathlib import Path
import os
from typing import Optional

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
    
    # Path construction methods
    def get_mmc_store_dir(self, brain_region: str, codebook: str) -> Path:
        """Get or create the MMC store directory for a given region/codebook."""
        store_dir = self.mmc_store_path / f"{brain_region}-{codebook}"
        store_dir.mkdir(parents=True, exist_ok=True)
        return store_dir
    
    def get_marker_dir(self, brain_region: str, codebook: str) -> Path:
        """Get marker directory."""
        return self.get_mmc_store_dir(brain_region, codebook) / "marker_dir"
    
    def get_precomp_stats_path(self, brain_region: str, codebook: str, ref_name: str) -> Path:
        """Get precomputed stats file path."""
        store_dir = self.get_mmc_store_dir(brain_region, codebook)
        return store_dir / f"{ref_name}_precomputed_stats.h5"
    
    def get_query_path(self, brain_region: str, codebook: str) -> Path:
        """Get query AnnData path."""
        store_dir = self.get_mmc_store_dir(brain_region, codebook)
        return store_dir / "query.h5ad"
    
    def get_codebook_path(self, codebook: str) -> Path:
        """Find codebook file by name."""
        import glob
        matches = glob.glob(str(self.gene_panel_path / f"*{codebook}*.csv"))
        if not matches:
            raise FileNotFoundError(f"No codebook found matching '{codebook}'")
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
        for path_name, path_obj in [
            ("mmc_store_path", self.mmc_store_path),
            ("anndata_store_path", self.anndata_store_path),
            ("annotations_store_path", self.annotations_store_path),
            ("zarr_storage_path", self.zarr_storage_path),
            ("gene_panel_path", self.gene_panel_path),
        ]:
            if not path_obj.exists():
                raise FileNotFoundError(f"{path_name} does not exist: {path_obj}")
```

### Phase 2: Preprocessing Component

#### Create `MMCPreprocessor` class
```python
import logging
from typing import List
from pathlib import Path
import pandas as pd
import anndata as ad
from scipy import sparse
import numpy as np

logger = logging.getLogger(__name__)

class MMCPreprocessor:
    """Handles all reference preprocessing for MapMyCells annotation."""
    
    def __init__(self, config: MMCConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def create_empty_query_adata(
        self,
        ref_path: str | Path,
        gene_list: List[str],
        output_path: str | Path
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
        """
        Precompute statistics for the reference AnnData.
        
        Returns path to precomputed stats file.
        """
        from cell_type_mapper.cli.precompute_stats_scrattch import (
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
        """
        Compute reference markers for the hierarchy.
        
        Returns path to reference markers file.
        """
        from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner
        
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
        """
        Downsample markers to those present in query.
        
        Returns path to query markers JSON file.
        """
        from cell_type_mapper.cli.query_markers import QueryMarkerRunner
        
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
    ) -> dict:
        """
        Complete setup: creates query AnnData, precomputes stats, and computes markers.
        
        Returns dict with paths to all intermediate files:
        {
            "precomp_stats_path": Path,
            "ref_markers_path": Path,
            "query_markers_path": Path,
            "query_adata_path": Path,
        }
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
```

### Phase 3: Annotation Component

#### Create `MMCAnnotator` class
```python
import pandas as pd
import anndata as ad
from pathlib import Path
import logging
from typing import Tuple

logger = logging.getLogger(__name__)

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
        """
        Run MapMyCells annotation on query data.
        
        Returns tuple of (extended_results_path, csv_results_path)
        """
        from cell_type_mapper.cli.from_specified_markers import (
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
        """
        Transfer annotation labels to AnnData object.
        
        Adds columns to adata.obs for each hierarchy level.
        """
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
```

### Phase 4: High-Level API

#### Update public functions
```python
def setup_annotation_pipeline(
    ref_path: str | Path,
    brain_region: str,
    codebook: str,
    hierarchy_list: List[str],
    codebook_path: str | Path = None,
    config: MMCConfig = None,
    **preprocessing_params,
) -> dict:
    """
    Complete setup for MapMyCells annotation.
    
    This creates all necessary reference data and markers for annotation.
    Should be run once per reference/codebook combination.
    
    Parameters:
    -----------
    ref_path : str or Path
        Path to reference single-cell RNA-seq AnnData
    brain_region : str
        Brain region identifier (e.g., "Whole brain", "Cortex")
    codebook : str
        Codebook identifier for matching
    hierarchy_list : List[str]
        Hierarchy levels for annotation (e.g., ["class", "subclass", "type"])
    codebook_path : str or Path, optional
        Explicit path to codebook CSV. If None, looks in gene_panel_path
    config : MMCConfig, optional
        Configuration object. If None, loads from environment variables
    **preprocessing_params
        Additional params (n_cpu, n_valid, n_per_utility, etc.)
    
    Returns:
    --------
    dict with paths to:
        - precomp_stats_path
        - ref_markers_path
        - query_markers_path
        - query_adata_path
    
    Example:
    --------
    >>> config = MMCConfig.from_env()
    >>> paths = setup_annotation_pipeline(
    ...     ref_path="/path/to/reference.h5ad",
    ...     brain_region="Whole brain",
    ...     codebook="mouse_brain_atlas",
    ...     hierarchy_list=["class", "subclass", "type"],
    ...     config=config,
    ... )
    """
    if config is None:
        config = MMCConfig.from_env()
    
    config.validate()
    
    preprocessor = MMCPreprocessor(config)
    return preprocessor.setup_annotation_pipeline(
        ref_path, brain_region, codebook, hierarchy_list, codebook_path, **preprocessing_params
    )


def annotate_query(
    query_adata: ad.AnnData,
    brain_region: str,
    codebook: str,
    output_dir: str | Path = None,
    config: MMCConfig = None,
    write_results: bool = True,
    **annotation_params,
) -> ad.AnnData:
    """
    Annotate query spatial data with cell types.
    
    Requires setup_annotation_pipeline to have been run first.
    
    Parameters:
    -----------
    query_adata : AnnData
        Query AnnData object (from spatial data)
    brain_region : str
        Brain region (must match setup call)
    codebook : str
        Codebook (must match setup call)
    output_dir : str or Path, optional
        Directory to save results. If None, creates in annotations_store_path
    config : MMCConfig, optional
        Configuration object. If None, loads from environment
    write_results : bool, default True
        Whether to write results to disk
    **annotation_params
        Additional params (bootstrap_factor, bootstrap_iterations, rng_seed, n_cpu, etc.)
    
    Returns:
    --------
    AnnData with added annotation columns in obs:
        - mmc_{hierarchy_level}: Predicted cell type
        - mmc_{hierarchy_level}_transfer_score: Bootstrap confidence score
    
    Example:
    --------
    >>> import anndata as ad
    >>> config = MMCConfig.from_env()
    >>> adata = ad.read_h5ad("/path/to/query.h5ad")
    >>> adata = annotate_query(
    ...     adata,
    ...     brain_region="Whole brain",
    ...     codebook="mouse_brain_atlas",
    ...     config=config,
    ... )
    >>> print(adata.obs.columns)
    """
    if config is None:
        config = MMCConfig.from_env()
    
    config.validate()
    
    if output_dir is None:
        output_dir = config.annotations_store_path / "default_output"
    
    # Run annotation
    annotator = MMCAnnotator(config)
    
    # Create temporary query AnnData at expected location
    temp_query_path = config.get_query_path(brain_region, codebook)
    # Note: This assumes setup was run. Otherwise this will fail with clear error
    
    extended_path, csv_path = annotator.annotate(
        str(temp_query_path), brain_region, codebook, output_dir, **annotation_params
    )
    
    # Transfer labels
    query_adata = annotator.transfer_labels(query_adata, csv_path)
    
    if write_results:
        # Save annotated AnnData
        output_dir = Path(output_dir)
        output_adata_path = output_dir / "query_annotated.h5ad"
        query_adata.write_h5ad(output_adata_path)
    
    return query_adata


def annotate_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params,
) -> ad.AnnData:
    """
    Annotate a single region from spatial data.
    
    Loads data from SpatialData zarr store, annotates, and returns.
    
    Example:
    --------
    >>> adata = annotate_region(
    ...     exp_name="experiment1",
    ...     reg_name="region_001",
    ...     prefix_name="default",
    ...     brain_region="Whole brain",
    ...     codebook="mouse_brain_atlas",
    ... )
    """
    if config is None:
        config = MMCConfig.from_env()
    
    # Load from zarr
    zarr_path = config.zarr_storage_path / exp_name / reg_name
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")
    
    # Set metadata
    adata.uns["experiment"] = exp_name
    adata.uns["segmentation"] = reg_name
    adata.uns["donor"] = reg_name.replace("region_", "")
    
    # Annotate
    output_dir = config.annotations_store_path / exp_name / reg_name / "mmc"
    return annotate_query(
        adata, brain_region, codebook, output_dir, config, **annotation_params
    )


def annotate_experiment(
    exp_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params,
) -> list[ad.AnnData]:
    """
    Annotate all regions in an experiment.
    
    Example:
    --------
    >>> adatas = annotate_experiment(
    ...     exp_name="experiment1",
    ...     prefix_name="default",
    ...     brain_region="Whole brain",
    ...     codebook="mouse_brain_atlas",
    ... )
    """
    import glob
    
    if config is None:
        config = MMCConfig.from_env()
    
    # Get all regions
    exp_path = config.zarr_storage_path / exp_name
    regions = glob.glob(str(exp_path / "region_*"))
    
    results = []
    for reg_path in regions:
        reg_name = Path(reg_path).name
        adata = annotate_region(
            exp_name, reg_name, prefix_name, brain_region, codebook, config, **annotation_params
        )
        results.append(adata)
    
    return results
```

---

## Benefits of This Refactoring

### 1. **Separation of Concerns**
   - Configuration isolated in `MMCConfig`
   - Preprocessing logic in `MMCPreprocessor`
   - Annotation logic in `MMCAnnotator`
   - High-level orchestration in public functions

### 2. **Improved Testability**
   - Each class can be tested independently
   - Mock configuration easier
   - Dependency injection via config parameter

### 3. **Better Error Handling**
   - Config validation upfront
   - Clear error messages for missing files/paths
   - No silent failures

### 4. **Extensibility**
   - Easy to add new preprocessing steps
   - Can subclass `MMCPreprocessor`/`MMCAnnotator` for variations
   - Path construction centralized for easy customization

### 5. **Better User Experience**
   - Simpler public API
   - Clearer parameter expectations
   - Comprehensive logging throughout

### 6. **Maintainability**
   - Reduced function size (easier to understand)
   - Single responsibility per class
   - Standardized configuration approach

---

## Migration Path

### Step 1: Add new classes alongside existing code
- Keep current functions working
- Old functions delegate to new classes

### Step 2: Update CLI to use new API
- `main.py` calls new public functions
- Maintains backwards compatibility through wrapping

### Step 3: Deprecate old internal functions
- Add deprecation warnings to `_setup_mmc`, `_mmc_runner`
- Document migration in docstrings

### Step 4: Remove old code
- Only after all dependent code migrated

---

## Implementation Priority

1. **High Priority**: `MMCConfig` - enables all other improvements
2. **High Priority**: `MMCPreprocessor` - handles setup
3. **Medium Priority**: `MMCAnnotator` - handles annotation
4. **Medium Priority**: New public API functions
5. **Low Priority**: Deprecate/remove old functions

---

## Example Usage After Refactoring

```python
# Setup (run once per reference/codebook combo)
config = MMCConfig.from_env()
setup_annotation_pipeline(
    ref_path="/path/to/reference.h5ad",
    brain_region="Whole brain",
    codebook="mouse_brain_atlas",
    hierarchy_list=["class", "subclass", "type"],
    config=config,
)

# Annotate single region
adata = annotate_region(
    exp_name="experiment1",
    reg_name="region_001",
    prefix_name="default",
    brain_region="Whole brain",
    codebook="mouse_brain_atlas",
    config=config,
)

# Annotate entire experiment
adatas = annotate_experiment(
    exp_name="experiment1",
    prefix_name="default",
    brain_region="Whole brain",
    codebook="mouse_brain_atlas",
    config=config,
)

# Use results
print(adata.obs.columns)  # mmc_class, mmc_class_transfer_score, mmc_subclass, etc.
```

