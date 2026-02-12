# SPIDA I Module CLI Usage Guide

The **I (Integration/Annotation) module** provides command-line tools for integrating external annotation methods (MapMyCells, ALLCools, MOSCOT) with spatially-resolved genomic data. It includes utilities for data backup, reference setup, and label transfer.

## Table of Contents

1. [Installation & Setup](#installation--setup)
2. [Configuration](#configuration)
3. [Commands Overview](#commands-overview)
4. [Detailed Command Documentation](#detailed-command-documentation)
5. [Common Workflows](#common-workflows)
6. [Troubleshooting](#troubleshooting)

---

## Installation & Setup

### Prerequisites

- Python 3.10+
- SPIDA package installed: `pip install -e .`
- Pixi environment configured with mmc dependencies

### Activate the I Module CLI

```bash
# Using pixi
pixi run -e mmc python -m spida.I

# Or directly with Python
python -m spida.I --help
```

### Verify Installation

```bash
python -m spida.I --help
```

You should see output listing all available commands (backup-adata, allcools-integration, mmc-setup, etc.)

---

## Configuration

The I module uses a hierarchical configuration system:

**Priority Order:** CLI arguments > .env config file > Environment variables > Defaults

### Environment Variables

Set these in your `.env` file or shell:

```bash
# Required paths
MMC_DIR=/path/to/mmc/markers
ANNDATA_STORE_PATH=/path/to/anndata/store
ANNOTATIONS_STORE_PATH=/path/to/annotations/store
ZARR_STORAGE_PATH=/path/to/zarr/storage
GENE_PANEL_PATH=/path/to/gene/panels

# Optional: MMC parameters (all have sensible defaults)
MMC_N_CPU=8                          # Number of CPUs for parallelization
MMC_BOOTSTRAP_FACTOR=0.8             # Bootstrap sampling factor
MMC_BOOTSTRAP_ITERATIONS=100         # Number of bootstrap iterations
MMC_RNG_SEED=13                      # Random seed for reproducibility
MMC_REF_NORM=log2CPM                 # Reference normalization method
MMC_QUERY_NORM=log2CPM               # Query normalization method
```

### Using a Config File

Pass a config file via CLI:

```bash
python -m spida.I --config my_config.env mmc-setup ...
```

---

## Commands Overview

| Command | Purpose | Module |
|---------|---------|--------|
| `backup-adata` | Backup AnnData from spatialdata to disk | Integration |
| `backup-adata-experiment` | Batch backup for entire experiment | Integration |
| `allcools-integration-region` | ALLCools integration for one region | ALLCools |
| `allcools-integration-experiment` | ALLCools integration for experiment | ALLCools |
| `plot-allcools-integration` | Generate plots from ALLCools results | ALLCools |
| `mmc-setup` | Setup MapMyCells reference data | MapMyCells |
| `mmc-annotation-region` | MapMyCells annotation for one region | MapMyCells |
| `mmc-annotation-experiment` | MapMyCells annotation for experiment | MapMyCells |
| `moscot-integration` | MOSCOT trajectory analysis | MOSCOT (placeholder) |

---

## Detailed Command Documentation

### Backup Commands

#### `backup-adata` - Backup Single Region

Create standalone AnnData files from spatialdata for external annotation tools.

**Usage:**
```bash
python -m spida.I backup-adata EXP_NAME REG_NAME PREFIX_NAME \
  [--adata_path PATH]
```

**Arguments:**
- `EXP_NAME`: Experiment name
- `REG_NAME`: Region name
- `PREFIX_NAME`: Prefix for keys in spatialdata object

**Options:**
- `--adata_path`: Path to AnnData store (uses ANNDATA_STORE_PATH if not provided)

**Example:**
```bash
python -m spida.I backup-adata brain_exp region_001 default_prefix
```

#### `backup-adata-experiment` - Batch Backup

Backup AnnData objects for all regions in an experiment.

**Usage:**
```bash
python -m spida.I backup-adata-experiment EXP_NAME PREFIX_NAME \
  [--adata_path PATH]
```

**Arguments:**
- `EXP_NAME`: Experiment name
- `PREFIX_NAME`: Prefix for keys in spatialdata

**Example:**
```bash
python -m spida.I backup-adata-experiment brain_exp default_prefix
```

---

### MapMyCells Commands

MapMyCells is a single-cell annotation tool that uses reference marker genes to classify query cells.

#### `mmc-setup` - Setup Reference Data

**ONE-TIME OPERATION**: Run once per reference identifier combination.

Create reference data, marker genes, and preprocessing files needed for annotation. This step is computationally intensive but only needs to be done once.

**Usage:**
```bash
python -m spida.I mmc-setup REF_PATH IDENTIFIER HIERARCHY_LIST \
  [--gene_names_path PATH] \
  [--gene_name_mapping_path PATH] \
  [--mmc_store_path PATH] \
  [--ref_norm NORM] \
  [--n_cpu N] \
  [--n_valid N] \
  [--n_per_utility N]
```

**Arguments:**
- `REF_PATH`: Path to reference single-cell RNA-seq AnnData file
- `IDENTIFIER`: Unique identifier for this reference (e.g., `motor_cortex_v1`)
- `HIERARCHY_LIST`: Space-separated hierarchy levels (e.g., `level1 level2 level3`)

**Options:**
- `--gene_names_path`: Path to gene_names.txt file (auto-detected from GENE_PANEL_PATH if not provided)
- `--gene_name_mapping_path`: Path to gene_name_mapping.json file (auto-detected if not provided)
- `--mmc_store_path`: Path to MMC markers store (uses MMC_DIR if not provided)
- `--ref_norm`: Normalization method (default: `log2CPM`)
- `--n_cpu`: Number of CPUs (default: 8)
- `--n_valid`: Number of valid markers per hierarchy level (default: 10)
- `--n_per_utility`: Number of markers per utility after downsampling (default: 10)

**Gene Name Mapping Workflow:**

The setup process now includes automatic gene name mapping:

1. Loads `gene_name_mapping.json` mapping reference gene names to MERFISH gene names
2. Applies mapping to reference gene names
3. Loads gene list from `gene_names.txt` (which uses MERFISH naming)
4. Creates query template with intersection of mapped genes and gene_names.txt

**Example:**
```bash
python -m spida.I mmc-setup \
  /data/reference_rna.h5ad \
  motor_cortex_v1 \
  class subclass cell_type \
  --n_cpu 16 \
  --n_valid 15 \
  --n_per_utility 12
```

**File Requirements:**

Your GENE_PANEL_PATH should contain:
- `motor_cortex_v1_gene_names.txt` - Text file with gene names (one per line)
- `motor_cortex_v1_gene_name_mapping.json` - JSON mapping reference genes to MERFISH genes

**Output:**
- Precomputed reference statistics
- Reference marker genes
- Query marker templates
- Query AnnData template

---

#### `mmc-annotation-region` - Annotate Single Region

Annotate cells in a single region using MapMyCells reference data.

**Usage:**
```bash
python -m spida.I mmc-annotation-region \
  EXP_NAME REG_NAME PREFIX_NAME IDENTIFIER \
  [--suffix SUFFIX] \
  [--mmc_store_path PATH] \
  [--anndata_store PATH] \
  [--annotation_store PATH] \
  [--zarr_store PATH] \
  [--n_cpu N] \
  [--bootstrap_factor FACTOR] \
  [--bootstrap_iterations N] \
  [--rng_seed SEED]
```

**Arguments:**
- `EXP_NAME`: Experiment name
- `REG_NAME`: Region name
- `PREFIX_NAME`: Prefix for keys in spatialdata
- `IDENTIFIER`: Reference identifier (must match setup)

**Options:**
- `--suffix`: Suffix for keys in spatialdata (default: `_filt`)
- `--mmc_store_path`: Path to MMC markers store
- `--anndata_store`: Path to AnnData store
- `--annotation_store`: Path to annotations store
- `--zarr_store`: Path to zarr storage
- `--n_cpu`: Number of CPUs (default: 1)
- `--bootstrap_factor`: Bootstrap sampling factor (default: 0.8)
- `--bootstrap_iterations`: Number of bootstrap iterations (default: 100)
- `--rng_seed`: Random seed (default: 13)

**Example:**
```bash
python -m spida.I mmc-annotation-region \
  brain_exp region_001 merged motor_cortex_v1 \
  --n_cpu 4 \
  --bootstrap_iterations 150
```

**Output:**
- Annotated AnnData with columns: `mmc_<level>`, `mmc_<level>_transfer_score`
- Extended results JSON file
- CSV results file

---

#### `mmc-annotation-experiment` - Batch Annotate Experiment

Annotate all regions in an experiment with the same reference.

**Usage:**
```bash
python -m spida.I mmc-annotation-experiment \
  EXP_NAME PREFIX_NAME IDENTIFIER \
  [--suffix SUFFIX] \
  [--mmc_store_path PATH] \
  [--anndata_store_path PATH] \
  [--annotations_store_path PATH] \
  [--n_cpu N] \
  [--bootstrap_factor FACTOR] \
  [--bootstrap_iterations N] \
  [--rng_seed SEED]
```

**Arguments:**
- `EXP_NAME`: Experiment name
- `PREFIX_NAME`: Prefix for keys in spatialdata
- `IDENTIFIER`: Reference identifier (must match setup)

**Options:** Same as `mmc-annotation-region`

**Example:**
```bash
python -m spida.I mmc-annotation-experiment \
  brain_exp merged motor_cortex_v1 \
  --n_cpu 8 \
  --bootstrap_iterations 200
```

**Note:** This processes all regions matching `region_*` pattern in the experiment.

---

### ALLCools Commands

ALLCools integrates spatial transcriptomics with single-cell reference data using cell-to-cell similarity.

#### `allcools-integration-region` - ALLCools Single Region

**Usage:**
```bash
python -m spida.I allcools-integration-region \
  EXP_NAME REG_NAME PREFIX_NAME REF_PATH \
  [--suffix SUFFIX] \
  [--anndata_store_path PATH] \
  [--annotations_store_path PATH] \
  [additional ALLCools parameters...]
```

**Arguments:**
- `EXP_NAME`: Experiment name
- `REG_NAME`: Region name
- `PREFIX_NAME`: Prefix for keys in spatialdata
- `REF_PATH`: Path to reference single-cell AnnData

**Options:**
- `--suffix`: Suffix for table keys (default: `_filt`)
- `--anndata_store_path`: Path to AnnData store
- `--annotations_store_path`: Path to annotations store

**Example:**
```bash
python -m spida.I allcools-integration-region \
  brain_exp region_001 filtered /data/reference_rna.h5ad \
  --rna_cell_type_column supercluster_name
```

---

#### `allcools-integration-experiment` - ALLCools Batch

**Usage:**
```bash
python -m spida.I allcools-integration-experiment \
  EXP_NAME PREFIX_NAME REF_PATH \
  [--suffix SUFFIX] \
  [--anndata_store_path PATH] \
  [--annotations_store_path PATH]
```

---

#### `plot-allcools-integration` - Generate Plots

**Usage:**
```bash
python -m spida.I plot-allcools-integration \
  EXP_NAME SEG_NAME DONOR \
  [--anndata_store_path PATH] \
  [--annotations_store_path PATH] \
  [--output_path PATH] \
  [--plot_joint_embeddings] \
  [--plot_c2c_transfer] \
  [--ref_cell_type_column COLUMN]
```

---

### MOSCOT Commands

**Status:** Placeholder - integration in progress

```bash
python -m spida.I moscot-integration
# Returns: NotImplementedError
```

---

## Common Workflows

### Workflow 1: MapMyCells Complete Annotation

Complete workflow for setting up and annotating with MapMyCells.

**Step 1: Prepare Gene Files**

Create the required gene files in your GENE_PANEL_PATH:

```bash
# gene_names.txt - one gene per line
cat > motor_cortex_v1_gene_names.txt << EOF
Slc17a7
Gad2
Vip
Satb2
Rorb
...
EOF

# gene_name_mapping.json - map reference names to MERFISH names
cat > motor_cortex_v1_gene_name_mapping.json << EOF
{
  "ENSMUSG00000020917": "Slc17a7",
  "ENSMUSG00000031127": "Gad2",
  "ENSMUSG00000080184": "Vip",
  "Gene_name_1": "Satb2",
  "Gene_name_2": "Rorb"
}
EOF
```

**Step 2: Setup (run once)**
```bash
# 1a. Prepare your reference RNA-seq data
# File: reference_rna.h5ad

# 1b. Run setup
python -m spida.I mmc-setup \
  /path/to/reference_rna.h5ad \
  motor_cortex_v1 \
  class subclass cell_type \
  --n_cpu 16
```

**Step 3: Annotate Regions**
```bash
# Annotate single region
python -m spida.I mmc-annotation-region \
  my_exp region_001 merged motor_cortex_v1 \
  --n_cpu 8

# Or annotate entire experiment
python -m spida.I mmc-annotation-experiment \
  my_exp merged motor_cortex_v1 \
  --n_cpu 8
```

**Step 4: Inspect Results**
```bash
# Results are in: ANNOTATIONS_STORE_PATH/my_exp/region_001/mmc/
# Files:
#   - csv_results.csv          (cell x annotation_level)
#   - extended_results.json    (detailed results)
#   - adata_*.h5ad             (AnnData with mmc_* columns)
```

---

### Workflow 2: Multiple References / Identifiers

Use different references for different brain regions.

```bash
# Setup reference 1 (motor cortex)
python -m spida.I mmc-setup ref_motor_cortex.h5ad \
  motor_cortex_v1 \
  class subclass cell_type

# Setup reference 2 (visual cortex)
python -m spida.I mmc-setup ref_visual_cortex.h5ad \
  visual_cortex_v1 \
  class subclass cell_type

# Annotate motor regions with motor reference
python -m spida.I mmc-annotation-region \
  exp region_motor01 merged motor_cortex_v1

# Annotate visual regions with visual reference
python -m spida.I mmc-annotation-region \
  exp region_visual01 merged visual_cortex_v1
```

---

### Workflow 3: Tuning Bootstrap Parameters

Fine-tune annotation confidence via bootstrap parameters.

```bash
# Conservative (high confidence, fewer predictions)
python -m spida.I mmc-annotation-region \
  exp region_001 merged identifier_v1 \
  --bootstrap_factor 0.5 \
  --bootstrap_iterations 200

# Aggressive (more predictions, lower threshold)
python -m spida.I mmc-annotation-region \
  exp region_001 merged identifier_v1 \
  --bootstrap_factor 0.9 \
  --bootstrap_iterations 50
```

---

### Workflow 4: Large-Scale Batch Processing

Annotate large experiments efficiently.

```bash
# Annotate all regions with parallelization
python -m spida.I mmc-annotation-experiment \
  large_exp merged motor_cortex_v1 \
  --n_cpu 16 \
  --bootstrap_iterations 100

# This will process all region_* directories in parallel
```

---

## Troubleshooting

### Issue: "Required environment variable not set"

**Problem:** Missing environment variables

**Solution:**
```bash
# Option 1: Set in .env file
cat > .env << EOF
MMC_DIR=/path/to/mmc
ANNDATA_STORE_PATH=/path/to/anndata
ANNOTATIONS_STORE_PATH=/path/to/annotations
ZARR_STORAGE_PATH=/path/to/zarr
GENE_PANEL_PATH=/path/to/genes
EOF

# Option 2: Export in shell
export MMC_DIR=/path/to/mmc
export ANNDATA_STORE_PATH=/path/to/anndata
# ... etc

# Option 3: Pass via CLI
python -m spida.I --config my_config.env mmc-setup ...
```

---

### Issue: "No gene_names found matching"

**Problem:** Gene names file not found

**Solution:**
```bash
# Check gene panel directory
ls /path/to/gene/panels/

# Provide explicit path
python -m spida.I mmc-setup ref.h5ad identifier \
  level1 level2 \
  --gene_names_path /path/to/specific/gene_names.txt \
  --gene_name_mapping_path /path/to/specific/mapping.json
```

---

### Issue: Memory errors during setup

**Problem:** Out of memory during preprocessing

**Solution:**
```bash
# Use fewer CPUs or reduce bootstrap parameters
python -m spida.I mmc-setup ref.h5ad level1 level2 region codebook \
  --n_cpu 4 \
  --n_valid 5 \
  --n_per_utility 5
```

---

### Issue: Slow annotation

**Problem:** Annotation is taking too long

**Solution:**
```bash
# Increase CPU count
python -m spida.I mmc-annotation-region \
  exp region_001 merged region codebook \
  --n_cpu 16

# Reduce bootstrap iterations
python -m spida.I mmc-annotation-region \
  exp region_001 merged region codebook \
  --bootstrap_iterations 50
```

---

### Issue: Configuration confusion

**Problem:** Not sure which config is being used

**Solution:**
```bash
# Enable debug logging
python -m spida.I --config debug.env mmc-setup ... -v

# This will show:
# - Which env file is loaded
# - Which environment variables are set
# - Which defaults are applied
```

---

## Advanced Usage

### Custom Reference Hierarchies

MapMyCells supports multi-level hierarchies. Define as needed:

```bash
# Two-level hierarchy
python -m spida.I mmc-setup ref.h5ad class cell_type region codebook

# Four-level hierarchy
python -m spida.I mmc-setup ref.h5ad \
  kingdom phylum class species region codebook

# Single level
python -m spida.I mmc-setup ref.h5ad cell_type region codebook
```

---

### Reproducible Results

Use fixed random seeds for reproducibility:

```bash
# Set seed for reproducibility
python -m spida.I mmc-annotation-region \
  exp region_001 merged region codebook \
  --rng_seed 42 \
  --bootstrap_iterations 100
```

---

### Parallel Processing

Leverage multi-core systems:

```bash
# During setup (reference processing)
python -m spida.I mmc-setup ref.h5ad level1 level2 region codebook \
  --n_cpu 32

# During annotation (query processing)
python -m spida.I mmc-annotation-region \
  exp region_001 merged region codebook \
  --n_cpu 32
```

---

## Exit Codes

- `0`: Success
- `1`: General error (missing path, validation failure)
- `2`: Missing required argument
- Non-zero: Specific exception from underlying library

---

## See Also

- [SPIDA Main Documentation](../../README.md)
- [S Module (Segmentation)](../S/S_CLI_USAGE.md)
- [P Module (Preprocessing)](../P/P_CLI_USAGE.md)
- [Configuration Guide](../config.md)
- [MapMyCells Documentation](https://github.com/broadinstitute/celltype_mapper)

---

## Support

For issues or questions:

1. Check [Troubleshooting](#troubleshooting) section
2. Review `SPIDA/logs/` for error messages
3. Check environment variables with `python -m spida.I --help`
4. Open an issue on GitHub with:
   - Full command executed
   - Error message/stack trace
   - Environment info (`pixi env show`)
