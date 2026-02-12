# MapMyCells Refactoring - Implementation Summary

## Overview
This document summarizes the immediate improvements made to `mmc.py` and provides a detailed elaboration on the proposed refactoring.

---

## Part 1: Immediate Improvements (✅ Completed)

### Changes Made to mmc.py

#### 1. Fixed kwargs Propagation
**Problem**: Kwargs were incorrectly passed as a dict parameter instead of unpacked.
```python
# Before (WRONG)
_setup_mmc(..., kwargs=kwargs)
_mmc_runner(..., kwargs=kwargs)

# After (CORRECT)
_setup_mmc(..., **kwargs)
_mmc_runner(..., **kwargs)
```

#### 2. Standardized Parameter Names
**Problem**: The typo `heirarchy_list` was used instead of `hierarchy_list` throughout.
```python
# Before
def _setup_mmc(..., heirarchy_list: list, ...):
    precomputation_config = {"hierarchy": heirarchy_list, ...}

# After
def _setup_mmc(..., hierarchy_list: list, ...):
    precomputation_config = {"hierarchy": hierarchy_list, ...}
```

**All occurrences updated**:
- `_setup_mmc()` function signature
- `setup_mmc()` function signature
- `_mmc_runner()` function signature
- `run_mmc()` function signature
- `mmc_annotation_region()` function signature
- `mmc_annotation_experiment()` function signature
- All logging statements
- All function calls

#### 3. Removed Redundant Wrapper
**Problem**: `mmc_setup()` was an unnecessary wrapper that just called `setup_mmc()`.
```python
# Before (REDUNDANT)
def mmc_setup(...):
    """Setup function for MapMyCells integration."""
    setup_mmc(...)
    logging.info("DONE")
    return 0

# After (SIMPLE ALIAS)
mmc_setup = setup_mmc
```

#### 4. Improved Parameter Consistency
- Converted all `BRAIN_REGION` to `brain_region` (snake_case)
- Converted all `CODEBOOK` to `codebook` (snake_case)
- Improved docstrings to be more detailed
- Fixed `ref_norm` hardcoding to use actual parameter value

**Impact**: More Pythonic naming conventions, easier to read and maintain.

---

## Part 2: Proposed Refactoring - Detailed Elaboration

### Architecture Overview

The refactoring proposes separating the monolithic script into **four distinct layers**:

```
┌─────────────────────────────────────────────────┐
│        High-Level Public API                     │
│  setup_annotation_pipeline()                     │
│  annotate_query()                                │
│  annotate_region()                               │
│  annotate_experiment()                           │
└─────────────────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────┐
│  MMCPreprocessor              MMCAnnotator       │
│  - Reference setup            - Query annotation│
│  - Marker computation         - Label transfer   │
│  - Stats precomputation                          │
└─────────────────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────┐
│              MMCConfig                           │
│  Centralized configuration & path management     │
│  from_env() | validate()                         │
└─────────────────────────────────────────────────┘
```

### Layer 1: Configuration Management (`MMCConfig`)

**Purpose**: Eliminate scattered environment variable handling and hardcoded paths.

**Key Features**:

```python
@dataclass
class MMCConfig:
    # All paths in one place
    mmc_store_path: Path
    anndata_store_path: Path
    annotations_store_path: Path
    zarr_storage_path: Path
    gene_panel_path: Path
    
    # Algorithm parameters
    n_cpu: int = 8
    bootstrap_factor: float = 0.8
    bootstrap_iterations: int = 100
    rng_seed: int = 13
    ref_norm: str = "log2CPM"
    query_norm: str = "log2CPM"
```

**Methods**:

| Method | Purpose |
|--------|---------|
| `get_mmc_store_dir()` | Gets/creates store for a region-codebook pair |
| `get_marker_dir()` | Gets marker directory for hierarchical annotation |
| `get_precomp_stats_path()` | Returns path to precomputed reference stats |
| `get_query_path()` | Returns path to query AnnData |
| `get_codebook_path()` | Finds codebook file by name |
| `from_env()` | Loads all configuration from environment variables |
| `validate()` | Checks that all required paths exist |

**Benefits**:
- ✅ No more scattered `os.getenv()` calls
- ✅ Path construction logic centralized (easy to adapt for different storage layouts)
- ✅ Early validation catches configuration errors
- ✅ Easy to override for testing (pass explicit config instead of relying on env vars)

**Example Usage**:
```python
# Load from environment
config = MMCConfig.from_env()

# Or create explicitly
config = MMCConfig(
    mmc_store_path="/data/mmc",
    anndata_store_path="/data/anndatas",
    annotations_store_path="/data/annotations",
    zarr_storage_path="/data/zarr",
    gene_panel_path="/data/gene_panels",
)

# Use in code
marker_dir = config.get_marker_dir("Whole brain", "mouse_atlas")
query_path = config.get_query_path("Whole brain", "mouse_atlas")
```

### Layer 2: Preprocessing Component (`MMCPreprocessor`)

**Purpose**: Encapsulate all reference preprocessing logic in a single, testable class.

**Separation from Original**:

```python
# Old mmc.py (mixed responsibilities)
def _setup_mmc(...):
    # 1. Creates query AnnData
    # 2. Configures paths
    # 3. Validates environment
    # 4. Precomputes stats
    # 5. Computes markers
    # ... all in one function

# New architecture
class MMCPreprocessor:
    def create_empty_query_adata()       # Responsibility 1
    def precompute_reference_stats()     # Responsibility 2
    def compute_reference_markers()      # Responsibility 3
    def compute_query_markers()          # Responsibility 4
    def setup_annotation_pipeline()      # Orchestrates 1-4
```

**Key Methods**:

```python
class MMCPreprocessor:
    def __init__(self, config: MMCConfig):
        """Initialize with centralized configuration."""
        self.config = config
    
    def create_empty_query_adata(
        self,
        ref_path: Path,
        gene_list: List[str],
        output_path: Path
    ) -> None:
        """
        Create minimal AnnData with reference genes ∩ query genes.
        
        This is required by cell_type_mapper to know which genes to use.
        It's essentially a template for what genes the query has.
        """
        # Implementation details...
        # - Reads reference var_names
        # - Intersects with available genes
        # - Creates sparse (1 x n_genes) AnnData
        # - Logs intersection statistics
    
    def precompute_reference_stats(
        self,
        ref_path: Path,
        hierarchy_list: List[str],
        brain_region: str,
        codebook: str
    ) -> Path:
        """
        Precompute statistics for efficient annotation.
        
        This step:
        - Reads reference AnnData
        - Computes per-class statistics
        - Stores in HDF5 format for fast loading
        - Uses config.n_cpu for parallelization
        
        Result is reusable for multiple queries.
        """
        # Calls cell_type_mapper's PrecomputationScrattchRunner
        # Uses config paths and parameters
        # Returns path to .h5 file
    
    def compute_reference_markers(
        self,
        precomp_path: Path,
        brain_region: str,
        codebook: str,
        query_path: Path = None,
        n_valid: int = 10
    ) -> Path:
        """
        Find differentially expressed genes for each cell type.
        
        This step:
        - Uses precomputed stats
        - Identifies markers per class
        - Stores hierarchically (class → subclass → type)
        - Optional: restricts to genes in query
        
        Result guides annotation step.
        """
        # Calls cell_type_mapper's ReferenceMarkerRunner
        # Saves to marker_dir/reference_markers.h5
    
    def compute_query_markers(
        self,
        ref_markers_path: Path,
        query_path: Path,
        brain_region: str,
        codebook: str,
        n_per_utility: int = 10
    ) -> Path:
        """
        Downsample markers to most informative genes.
        
        This step:
        - Takes reference markers
        - Restricts to genes in query codebook
        - Selects top-N genes per hierarchy level
        - Optimizes for speed/accuracy trade-off
        
        Result is used during annotation.
        """
        # Calls cell_type_mapper's QueryMarkerRunner
        # Saves to marker_dir/calc_markers.json
    
    def setup_annotation_pipeline(
        self,
        ref_path: Path,
        brain_region: str,
        codebook: str,
        hierarchy_list: List[str],
        codebook_path: Path = None,
        **preprocessing_params
    ) -> Dict[str, Path]:
        """
        Orchestrate complete setup.
        
        Calls methods 1-4 in correct order, with error handling.
        Returns dict with all intermediate file paths.
        
        This is THE function to call for setup!
        """
        # 1. Load codebook genes
        # 2. Create query template
        # 3. Precompute reference stats
        # 4. Compute reference markers
        # 5. Compute query markers
        # 6. Return all paths
```

**Why This Matters**:
- ✅ Each step is testable independently
- ✅ Clear progress logging at each stage
- ✅ Can skip early steps if they already ran (check for file existence)
- ✅ Easy to add new preprocessing steps later
- ✅ Follows single responsibility principle

### Layer 3: Annotation Component (`MMCAnnotator`)

**Purpose**: Encapsulate query annotation logic separately from preprocessing.

```python
class MMCAnnotator:
    def __init__(self, config: MMCConfig):
        """Initialize with centralized configuration."""
        self.config = config
    
    def annotate(
        self,
        query_adata_path: Path,
        brain_region: str,
        codebook: str,
        output_dir: Path,
        **annotation_params
    ) -> Tuple[Path, Path]:
        """
        Run MapMyCells annotation.
        
        Uses precomputed markers from setup phase to annotate query.
        
        Parameters from annotation_params:
        - bootstrap_factor: fraction of cells to use per iteration
        - bootstrap_iterations: number of bootstrap samples
        - rng_seed: random seed for reproducibility
        - n_cpu: number of processors
        
        Returns paths to:
        - extended_results.json (cell-level annotations + scores)
        - csv_results.csv (aggregated results)
        """
        # Locates precomputed markers and stats from config
        # Runs cell_type_mapper's FromSpecifiedMarkersRunner
        # Returns result file paths
    
    def transfer_labels(
        self,
        adata: ad.AnnData,
        csv_results_path: Path
    ) -> ad.AnnData:
        """
        Transfer annotation results to AnnData.obs.
        
        For each hierarchy level, adds columns:
        - mmc_{level}: predicted cell type name
        - mmc_{level}_transfer_score: bootstrap confidence (0-1)
        
        Example output:
        adata.obs columns include:
        - mmc_class: "Excitatory neuron"
        - mmc_class_transfer_score: 0.95
        - mmc_subclass: "Layer 2/3 pyramidal"
        - mmc_subclass_transfer_score: 0.87
        - mmc_type: "L2/3 intratelencephalic"
        - mmc_type_transfer_score: 0.72
        """
        # Parse CSV results
        # Extract hierarchy levels from column names
        # Add to adata.obs
        # Return modified adata
```

**Why Separate**:
- ✅ Preprocessing runs infrequently (setup once per reference)
- ✅ Annotation runs frequently (once per region/experiment)
- ✅ Different parameter sets for each phase
- ✅ Can parallelize annotations across multiple queries using same reference

### Layer 4: High-Level Public API

**Purpose**: Provide simple, user-friendly functions for common operations.

```python
# Setup (run once per reference-codebook combination)
def setup_annotation_pipeline(
    ref_path: Path,
    brain_region: str,
    codebook: str,
    hierarchy_list: List[str],
    codebook_path: Path = None,
    config: MMCConfig = None,
    **preprocessing_params
) -> Dict[str, Path]:
    """
    One-line setup for MapMyCells.
    
    Returns paths to all generated files for reproducibility.
    """
    # Implementation
    config = config or MMCConfig.from_env()
    config.validate()
    preprocessor = MMCPreprocessor(config)
    return preprocessor.setup_annotation_pipeline(
        ref_path, brain_region, codebook, hierarchy_list, codebook_path, **preprocessing_params
    )

# Annotate a single AnnData
def annotate_query(
    query_adata: ad.AnnData,
    brain_region: str,
    codebook: str,
    output_dir: Path = None,
    config: MMCConfig = None,
    write_results: bool = True,
    **annotation_params
) -> ad.AnnData:
    """
    Annotate spatial data.
    
    Takes AnnData, adds annotation columns, returns modified AnnData.
    Optionally writes results to disk.
    """
    # Implementation
    config = config or MMCConfig.from_env()
    config.validate()
    
    annotator = MMCAnnotator(config)
    temp_query_path = config.get_query_path(brain_region, codebook)
    
    extended_path, csv_path = annotator.annotate(
        temp_query_path, brain_region, codebook, output_dir, **annotation_params
    )
    query_adata = annotator.transfer_labels(query_adata, csv_path)
    
    if write_results:
        query_adata.write_h5ad(Path(output_dir) / "query_annotated.h5ad")
    
    return query_adata

# Annotate a single region from SpatialData
def annotate_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params
) -> ad.AnnData:
    """
    Annotate one region from experiment.
    
    - Loads from zarr store
    - Runs annotation
    - Returns annotated AnnData
    """
    # Implementation
    config = config or MMCConfig.from_env()
    zarr_path = config.zarr_storage_path / exp_name / reg_name
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")
    
    output_dir = config.annotations_store_path / exp_name / reg_name / "mmc"
    return annotate_query(
        adata, brain_region, codebook, output_dir, config, **annotation_params
    )

# Annotate entire experiment (all regions)
def annotate_experiment(
    exp_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params
) -> List[ad.AnnData]:
    """
    Annotate all regions in experiment.
    
    Loops through regions, annotates each, returns list of AnnDatas.
    """
    # Implementation
    config = config or MMCConfig.from_env()
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

## Part 3: Implementation Timeline

### Phase 1 (Immediate): Already Done ✅
- [x] Fix kwargs propagation
- [x] Standardize parameter names
- [x] Remove redundant wrapper

### Phase 2 (Short-term): Add New Classes
1. Add `MMCConfig` class to mmc.py
2. Add `MMCPreprocessor` class to mmc.py
3. Add `MMCAnnotator` class to mmc.py
4. Add new public API functions (wrap old ones)

**No breaking changes** - old functions still work!

### Phase 3 (Medium-term): Update CLI
- Update `main.py` to use new public API
- Deprecate old internal functions with warnings

### Phase 4 (Long-term): Cleanup
- Remove deprecated functions
- Update tests to use new API
- Update documentation

---

## Part 4: Usage Comparison

### Before Refactoring
```python
# Multiple low-level calls, environment dependency
setup_mmc(ref_path, BRAIN_REGION, CODEBOOK, heirarchy_list)
mmc_annotation_region(exp_name, reg_name, prefix_name, BRAIN_REGION, CODEBOOK)
```

### After Refactoring
```python
# Simple, high-level calls
config = MMCConfig.from_env()  # Explicit config
setup_annotation_pipeline(ref_path, brain_region, codebook, hierarchy_list, config=config)
annotate_region(exp_name, reg_name, prefix_name, brain_region, codebook, config=config)
```

---

## Part 5: Benefits Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Configuration** | Scattered env vars | Centralized `MMCConfig` |
| **Error Handling** | Silent failures | Early validation |
| **Testability** | Hard (env dependent) | Easy (config injectable) |
| **Path Construction** | Hardcoded strings | `Config` methods |
| **Logging** | Sporadic | Comprehensive (each step) |
| **Single Responsibility** | Mixed (100+ lines per function) | Separated (10-30 lines each) |
| **Documentation** | Minimal | Detailed docstrings |
| **Extensibility** | Hard (monolithic) | Easy (separate classes) |
| **Parameter Names** | Inconsistent (heirarchy_list) | Consistent (hierarchy_list) |

---

## Files Created

1. **`REFACTORING_PROPOSAL.md`** - Complete architectural proposal
2. **`mmc_refactored.py`** - Reference implementation of new classes

Both files are located in `/home/x-aklein2/projects/aklein/SPIDA/src/spida/I/`

---

## Next Steps

1. **Review** the refactoring proposal
2. **Copy** classes from `mmc_refactored.py` into `mmc.py`
3. **Update** `main.py` CLI to use new public functions
4. **Test** with existing workflows
5. **Deprecate** old functions gradually
6. **Remove** deprecated code after grace period

