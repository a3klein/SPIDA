# Before & After Code Examples

This document shows concrete examples of how the refactoring improves code clarity and maintainability.

---

## Example 1: Configuration Management

### Before (Scattered Environment Variables)
```python
# Old code scattered across functions
def setup_mmc(ref_path, brain_region, codebook, hierarchy_list, ...):
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")  # ❌ Hardcoded env var
        if not mmc_store_path:
            raise ValueError("Please provide a mmc_store_path...")
    
    if not codebook_path:
        gene_panel_path = os.getenv("GENE_PANEL_PATH")  # ❌ Another env var
        if not gene_panel_path:
            raise ValueError("Please provide a codebook path...")
        codebook_path = glob.glob(f"{gene_panel_path}/*{codebook}*.csv")[0]  # ❌ Fragile path construction

def run_mmc(query_adata, brain_region, codebook, ...):
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")  # ❌ Duplicated code
        if not mmc_store_path:
            raise ValueError("...")  # ❌ Duplicated error handling
    
    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")  # ❌ Different env var
        if not anndata_store_path:
            raise ValueError("...")
    
    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")  # ❌ Another one
        if not annotations_store_path:
            raise ValueError("...")

def mmc_annotation_experiment(exp_name, prefix_name, brain_region, codebook, ...):
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")  # ❌ Defaults to hardcoded path
    exp_path = Path(f"{zarr_store}/{exp_name}")
    # ...
```

**Problems**:
- ❌ Environment variable handling repeated in multiple functions
- ❌ No validation that all required vars are set
- ❌ Path construction logic is fragile and hardcoded
- ❌ Default values hardcoded in code
- ❌ Hard to test (must set env vars)
- ❌ Hard to use with different configurations

### After (Centralized Configuration)
```python
# New code - all in one place
@dataclass
class MMCConfig:
    # All paths declared upfront
    mmc_store_path: Path
    anndata_store_path: Path
    annotations_store_path: Path
    zarr_storage_path: Path
    gene_panel_path: Path
    
    # All parameters with sensible defaults
    n_cpu: int = 8
    bootstrap_factor: float = 0.8
    bootstrap_iterations: int = 100
    rng_seed: int = 13
    ref_norm: str = "log2CPM"
    query_norm: str = "log2CPM"
    
    @classmethod
    def from_env(cls) -> "MMCConfig":
        """Load from environment with comprehensive validation."""
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
        
        config_dict["n_cpu"] = int(os.getenv("MMC_N_CPU", "8"))
        config_dict["bootstrap_factor"] = float(os.getenv("MMC_BOOTSTRAP_FACTOR", "0.8"))
        # ... more optional vars
        
        return cls(**config_dict)
    
    def get_codebook_path(self, codebook: str) -> Path:
        """Find codebook file - centralized logic."""
        matches = glob.glob(str(self.gene_panel_path / f"*{codebook}*.csv"))
        if not matches:
            raise FileNotFoundError(f"No codebook found matching '{codebook}'")
        return Path(matches[0])
    
    def get_marker_dir(self, brain_region: str, codebook: str) -> Path:
        """Construct marker directory path."""
        marker_dir = self.mmc_store_path / f"{brain_region}-{codebook}" / "marker_dir"
        marker_dir.mkdir(parents=True, exist_ok=True)
        return marker_dir
    
    def validate(self) -> None:
        """Validate all paths exist."""
        for name, path in [
            ("mmc_store_path", self.mmc_store_path),
            ("gene_panel_path", self.gene_panel_path),
            # ...
        ]:
            if not path.exists():
                raise FileNotFoundError(f"{name} does not exist: {path}")


# Usage: simple, testable, and explicit
def setup_mmc(ref_path, brain_region, codebook, hierarchy_list, config=None):
    if config is None:
        config = MMCConfig.from_env()
    
    config.validate()  # ✅ All validation in one place
    
    codebook_path = config.get_codebook_path(codebook)  # ✅ Centralized logic
    marker_dir = config.get_marker_dir(brain_region, codebook)  # ✅ Centralized logic
    
    # Rest of function...

# Usage in tests: inject configuration
def test_mmc_setup():
    config = MMCConfig(
        mmc_store_path=Path("/tmp/test_mmc"),
        anndata_store_path=Path("/tmp/test_adata"),
        annotations_store_path=Path("/tmp/test_annot"),
        zarr_storage_path=Path("/tmp/test_zarr"),
        gene_panel_path=Path("/tmp/test_panels"),
    )
    
    # No need to set environment variables!
    setup_mmc(ref_path, "brain", "codebook", ["class"], config=config)
```

**Benefits**:
- ✅ All configuration in one place
- ✅ Comprehensive validation upfront
- ✅ Easy to test (pass explicit config)
- ✅ Can use different configurations without env vars
- ✅ No code duplication
- ✅ Clear defaults

---

## Example 2: Single Responsibility

### Before (Mixed Concerns)
```python
# ❌ One function does EVERYTHING
def _setup_mmc(
    ref_path, mmc_store_path, brain_region, codebook,
    codebook_path, hierarchy_list, ref_norm, **kwargs
):
    """Setup preprocessing - but what does it do exactly?"""
    
    out_dir = f"{mmc_store_path}/{brain_region}-{codebook}"
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    marker_dir = f"{out_dir}/marker_dir"
    query_path = f"{out_dir}/query.h5ad"
    
    gene_list = pd.read_csv(codebook_path)["name"]
    gene_list = gene_list[~gene_list.str.startswith("Blank-")]
    
    # Task 1: Create empty AnnData
    adata = ad.read_h5ad(ref_path, backed="r")
    var_index = adata.var.index
    query_index = var_index.intersection(pd.Index(gene_list))
    query = ad.AnnData(X=sparse.csr_matrix(np.zeros((1, len(query_index)))), ...)
    query.write_h5ad(query_path)
    
    # Task 2: Precompute stats
    n_cpu = kwargs.get("n_cpu", 8)
    precomp_output_path = f"{out_dir}/{Path(ref_path).stem}_precomputed_stats.h5"
    precomputation_config = {"hierarchy": hierarchy_list, ...}
    precomputation_runner = PrecomputationScrattchRunner(args=[], input_data=precomputation_config)
    precomputation_runner.run()
    
    # Task 3: Compute reference markers
    n_valid = kwargs.get("n_valid", 10)
    reference_config = {"precomputed_path_list": [precomp_output_path], ...}
    reference_runner = ReferenceMarkerRunner(args=[], input_data=reference_config)
    reference_runner.run()
    
    # Task 4: Compute query markers
    n_per_utility = kwargs.get("n_per_utility", 10)
    query_config = {"output_path": f"{marker_dir}/calc_markers.json", ...}
    query_runner = QueryMarkerRunner(args=[], input_data=query_config)
    query_runner.run()
```

**Problems**:
- ❌ 50+ lines doing 4 different things
- ❌ Hard to understand flow
- ❌ Hard to test individual steps
- ❌ Kwargs extraction scattered throughout
- ❌ Can't reuse individual steps
- ❌ Can't easily add error recovery for specific steps

### After (Separated Concerns)
```python
# ✅ Each class has ONE responsibility
class MMCPreprocessor:
    def __init__(self, config: MMCConfig):
        self.config = config
    
    # ✅ Step 1: Only create query AnnData
    def create_empty_query_adata(
        self, ref_path: Path, gene_list: List[str], output_path: Path
    ) -> None:
        """Create empty AnnData with gene intersection."""
        ref_adata = ad.read_h5ad(ref_path, backed="r")
        var_index = ref_adata.var.index
        query_index = var_index.intersection(pd.Index(gene_list))
        query = ad.AnnData(X=sparse.csr_matrix(np.zeros((1, len(query_index)))), ...)
        query.write_h5ad(output_path)
    
    # ✅ Step 2: Only precompute stats
    def precompute_reference_stats(
        self, ref_path: Path, hierarchy_list: List[str],
        brain_region: str, codebook: str
    ) -> Path:
        """Precompute reference statistics."""
        output_path = self.config.get_precomp_stats_path(brain_region, codebook, Path(ref_path).stem)
        precomputation_config = {"hierarchy": hierarchy_list, ...}
        runner = PrecomputationScrattchRunner(args=[], input_data=precomputation_config)
        runner.run()
        return output_path
    
    # ✅ Step 3: Only compute reference markers
    def compute_reference_markers(
        self, precomp_path: Path, brain_region: str, codebook: str,
        query_path: Path = None, n_valid: int = 10
    ) -> Path:
        """Compute reference markers."""
        marker_dir = self.config.get_marker_dir(brain_region, codebook)
        reference_config = {"precomputed_path_list": [str(precomp_path)], ...}
        runner = ReferenceMarkerRunner(args=[], input_data=reference_config)
        runner.run()
        return marker_dir / "reference_markers.h5"
    
    # ✅ Step 4: Only compute query markers
    def compute_query_markers(
        self, ref_markers_path: Path, query_path: Path,
        brain_region: str, codebook: str, n_per_utility: int = 10
    ) -> Path:
        """Downsample markers to query genes."""
        marker_dir = self.config.get_marker_dir(brain_region, codebook)
        query_config = {"output_path": str(marker_dir / "calc_markers.json"), ...}
        runner = QueryMarkerRunner(args=[], input_data=query_config)
        runner.run()
        return marker_dir / "calc_markers.json"
    
    # ✅ Step 5: Orchestrate all steps
    def setup_annotation_pipeline(
        self, ref_path: Path, brain_region: str, codebook: str,
        hierarchy_list: List[str], codebook_path: Path = None, **params
    ) -> Dict[str, Path]:
        """Complete setup - calls steps 1-4 in order."""
        gene_list = pd.read_csv(codebook_path)["name"]
        gene_list = gene_list[~gene_list.str.startswith("Blank-")]
        
        query_path = self.config.get_query_path(brain_region, codebook)
        self.create_empty_query_adata(ref_path, gene_list, query_path)
        
        precomp_path = self.precompute_reference_stats(ref_path, hierarchy_list, brain_region, codebook)
        ref_markers = self.compute_reference_markers(precomp_path, brain_region, codebook, query_path)
        query_markers = self.compute_query_markers(ref_markers, query_path, brain_region, codebook)
        
        return {
            "precomp_stats_path": precomp_path,
            "ref_markers_path": ref_markers,
            "query_markers_path": query_markers,
        }

# ✅ Usage: clear, testable, reusable
preprocessor = MMCPreprocessor(config)

# Can test each step independently
preprocessor.create_empty_query_adata(ref_path, genes, output)
preprocessor.precompute_reference_stats(ref_path, hierarchy, region, codebook)

# Or run complete pipeline
results = preprocessor.setup_annotation_pipeline(...)
```

**Benefits**:
- ✅ Each method has ONE clear purpose (10-20 lines)
- ✅ Can test each step independently
- ✅ Can reuse individual steps
- ✅ Flow is clear and easy to follow
- ✅ Easy to add error recovery for specific steps
- ✅ Easy to add new preprocessing steps

---

## Example 3: Parameter Handling

### Before (Kwargs Magic)
```python
# ❌ Fragile and unclear
def _setup_mmc(..., **kwargs):
    n_cpu = kwargs.get("n_cpu", 8)  # ❌ Magic number, where from?
    n_per_utility = kwargs.get("n_per_utility", 10)  # ❌ What if typo in kwarg?
    n_valid = kwargs.get("n_valid", 10)  # ❌ No IDE autocomplete
    
    precomputation_config = {
        "n_processors": n_cpu,
        # ...
    }

def _mmc_runner(..., **kwargs):
    n_cpu = kwargs.get("n_cpu", 1)  # ❌ Different default!
    bootstrap_factor = kwargs.get("bootstrap_factor", 0.8)
    bootstrap_iteration = kwargs.get("bootstrap_iteration", 100)  # ❌ Typo: should be "bootstrap_iterations"?
    rng = kwargs.get("rng_seed", 13)

# Usage: what does this actually pass?
setup_mmc(ref, region, codebook, hierarchy, bootstrap_factor=0.8, n_cpu=16)
```

### After (Explicit Configuration)
```python
# ✅ Clear, type-checked, with IDE support
@dataclass
class MMCConfig:
    # All parameters declared with types and defaults
    n_cpu: int = 8
    bootstrap_factor: float = 0.8
    bootstrap_iterations: int = 100  # ✅ Consistent naming
    rng_seed: int = 13
    ref_norm: str = "log2CPM"
    query_norm: str = "log2CPM"

# ✅ IDE autocomplete, type checking, clear usage
preprocessor = MMCPreprocessor(config)
preprocessor.precompute_reference_stats(
    ref_path, hierarchy, brain_region, codebook
)

# ✅ Annotation parameters explicit in method signature
annotator.annotate(
    query_path, brain_region, codebook, output_dir,
    bootstrap_factor=0.8,  # ✅ IDE suggests this
    bootstrap_iterations=100,  # ✅ Clear name
    rng_seed=13,  # ✅ Type checked
    n_cpu=16  # ✅ Explicit parameter
)
```

**Benefits**:
- ✅ IDE autocomplete and type checking
- ✅ Documentation visible in IDE
- ✅ No silent parameter typos
- ✅ Clear defaults
- ✅ Easy to find where parameters are used

---

## Example 4: Public API

### Before (Multiple Similar Functions)
```python
# ❌ Three functions that do almost the same thing
def mmc_annotation_region(exp_name, reg_name, prefix_name, BRAIN_REGION, CODEBOOK, ...):
    """Annotate one region."""
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")
    return run_mmc(adata, BRAIN_REGION, CODEBOOK, ...)

def mmc_annotation_experiment(exp_name, prefix_name, BRAIN_REGION, CODEBOOK, ...):
    """Annotate entire experiment."""
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions:
        reg_name = reg.split("/")[-1]
        mmc_annotation_region(exp_name, reg_name, prefix_name, BRAIN_REGION, CODEBOOK, ...)

# ❌ Confusing wrapper function
def mmc_setup(ref_path, heirarchy_list, BRAIN_REGION, CODEBOOK, ...):
    """Wrapper that just calls setup_mmc."""
    setup_mmc(ref_path, BRAIN_REGION, CODEBOOK, codebook_path, heirarchy_list, ...)
    logging.info("DONE")
    return 0
```

**Problems**:
- ❌ `mmc_setup` is just an alias (why have two names?)
- ❌ Region vs experiment functions similar but separate
- ❌ Parameter naming inconsistent (BRAIN_REGION vs brain_region)
- ❌ Return values inconsistent (tuple vs 0)
- ❌ Uses different underlying functions (_run_mmc vs run_mmc)

### After (Consistent, Clear API)
```python
# ✅ One-line setup
def setup_annotation_pipeline(
    ref_path: Path,
    brain_region: str,
    codebook: str,
    hierarchy_list: List[str],
    codebook_path: Path = None,
    config: MMCConfig = None,
    **preprocessing_params
) -> Dict[str, Path]:
    """Setup for MapMyCells annotation (once per reference)."""
    config = config or MMCConfig.from_env()
    config.validate()
    preprocessor = MMCPreprocessor(config)
    return preprocessor.setup_annotation_pipeline(
        ref_path, brain_region, codebook, hierarchy_list, codebook_path, **preprocessing_params
    )

# ✅ Annotate a single spatial region
def annotate_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params
) -> ad.AnnData:
    """Annotate one region from spatial data."""
    config = config or MMCConfig.from_env()
    zarr_path = config.zarr_storage_path / exp_name / reg_name
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table")
    
    output_dir = config.annotations_store_path / exp_name / reg_name / "mmc"
    return annotate_query(adata, brain_region, codebook, output_dir, config, **annotation_params)

# ✅ Annotate entire experiment (all regions)
def annotate_experiment(
    exp_name: str,
    prefix_name: str,
    brain_region: str,
    codebook: str,
    config: MMCConfig = None,
    **annotation_params
) -> List[ad.AnnData]:
    """Annotate all regions in an experiment."""
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

# ✅ Usage: clear, consistent
config = MMCConfig.from_env()

# Setup once
setup_annotation_pipeline("ref.h5ad", "brain", "atlas", ["class", "subclass", "type"], config=config)

# Annotate one region
adata = annotate_region("exp1", "region_001", "default", "brain", "atlas", config=config)

# Annotate all regions
adatas = annotate_experiment("exp1", "default", "brain", "atlas", config=config)
```

**Benefits**:
- ✅ Single responsibility per function
- ✅ Consistent parameter naming
- ✅ Consistent return types
- ✅ No confusing aliases
- ✅ Clear separation of concern (setup vs annotate)
- ✅ Easy to document
- ✅ Easy to test

---

## Summary of Improvements

| Issue | Before | After |
|-------|--------|-------|
| **Configuration** | Scattered env vars in multiple functions | Centralized `MMCConfig` class |
| **Path Construction** | Hardcoded f-strings throughout | `Config.get_*_path()` methods |
| **Parameter Passing** | Kwargs dict with magic keys | Dataclass with type hints |
| **Function Size** | 50+ lines doing multiple things | 10-20 lines doing one thing |
| **Testability** | Requires mocking env vars | Can pass explicit config |
| **Code Reuse** | Can't reuse individual steps | Each step callable independently |
| **API Clarity** | Multiple similar functions | Clear setup/annotate distinction |
| **IDE Support** | No autocomplete for kwargs | Full autocomplete and type checking |
| **Error Messages** | Silent failures | Early validation with clear errors |
| **Documentation** | Unclear what parameters do | Explicit in type hints and docstrings |

