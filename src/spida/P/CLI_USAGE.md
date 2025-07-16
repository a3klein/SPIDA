# SPIDA P Module CLI Usage Guide

The SPIDA P (Pipeline) module provides a command-line interface for filtering cells and setting up AnnData objects in the BICAN project. This CLI tool replaces the previous Fire-based interface with a more structured argparse-based approach.

## Installation and Setup

Before using the CLI, ensure you have the SPIDA package installed and the necessary environment variables set:

- `DEF_CUTOFFS_PATH`: Default path to filtering cutoffs JSON file
- `ZARR_STORAGE_PATH`: Path to ZARR storage directory
- `IMAGE_STORE_PATH`: Path to image storage directory

## General Usage

```bash
python -m spida.P.cli [command] [arguments] [options]
```

### Get Version Information

```bash
python -m spida.P.cli --version
```

### Get Help

```bash
python -m spida.P.cli --help
python -m spida.P.cli [command] --help
```

## Available Commands

### 1. Filter Cells for a Specific Region

Filter cells for a specific region in an experiment.

```bash
python -m spida.P.cli filter-cells-region EXP_NAME REG_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--cutoffs-path PATH`: Path to the cutoffs JSON file (default: uses DEF_CUTOFFS_PATH env var)
- `--plot`: Generate plots of the filtering results
- `--image-path PATH`: Path to save the plot (default: auto-generated)

**Example:**
```bash
python -m spida.P.cli filter-cells-region experiment1 region_001 default --plot
```

### 2. Filter Cells for All Regions

Filter cells for all regions in an experiment.

```bash
python -m spida.P.cli filter-cells-all EXP_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--cutoffs-path PATH`: Path to the cutoffs JSON file
- `--plot`: Generate plots for all regions
- `--image-path PATH`: Base path for saving plots

**Example:**
```bash
python -m spida.P.cli filter-cells-all experiment1 default --plot
```

### 3. Write AnnData Objects to Disk

Write AnnData objects to disk for a specific experiment and region.

```bash
python -m spida.P.cli write-adata EXP_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment

**Options:**
- `--reg-name REG_NAME`: Name of the region (optional, if not provided writes all regions)
- `--prefix-names PREFIX1 PREFIX2 ...`: List of prefix names (required)
- `--output-path PATH`: Output directory path (default: `/ceph/cephatlas/aklein/bican/data/anndatas/`)

**Examples:**
```bash
# Write specific region
python -m spida.P.cli write-adata experiment1 --reg-name region_001 --prefix-names default proseg

# Write all regions
python -m spida.P.cli write-adata experiment1 --prefix-names default proseg cellpose
```

### 4. Setup AnnData for a Specific Region

Setup the AnnData object for downstream analysis for a specific region. This involves normalizing data, calculating PCA, UMAP, t-SNE, and Leiden clusters.

```bash
python -m spida.P.cli setup-adata-region EXP_NAME REG_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--plot`: Generate plots of the setup results
- `--image-path PATH`: Path to save the plot

**Example:**
```bash
python -m spida.P.cli setup-adata-region experiment1 region_001 default --plot
```

### 5. Setup AnnData for All Regions

Setup the AnnData objects for all regions in an experiment.

```bash
python -m spida.P.cli setup-adata-all EXP_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--plot`: Generate plots for all regions
- `--image-path PATH`: Base path for saving plots

**Example:**
```bash
python -m spida.P.cli setup-adata-all experiment1 default --plot
```

### 6. Plot Filtering Results

Generate plots for filtering results of a specific region.

```bash
python -m spida.P.cli plot-filtering-region EXP_NAME REG_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--image-path PATH`: Path to save the plot

**Example:**
```bash
python -m spida.P.cli plot-filtering-region experiment1 region_001 default --image-path /path/to/output.pdf
```

### 7. Plot Setup Results

Generate plots for setup results of a specific region.

```bash
python -m spida.P.cli plot-setup-region EXP_NAME REG_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--image-path PATH`: Path to save the plot

**Example:**
```bash
python -m spida.P.cli plot-setup-region experiment1 region_001 default --image-path /path/to/setup_plots.pdf
```

## Typical Workflow

Here's a typical workflow for processing an experiment:

```bash
# 1. Filter cells for all regions
python -m spida.P.cli filter-cells-all experiment1 default --plot

# 2. Setup AnnData for all regions
python -m spida.P.cli setup-adata-all experiment1 default --plot

# 3. Write AnnData objects to disk
python -m spida.P.cli write-adata experiment1 --prefix-names default
```

For processing a single region:

```bash
# 1. Filter cells for a specific region
python -m spida.P.cli filter-cells-region experiment1 region_001 default --plot

# 2. Setup AnnData for the region
python -m spida.P.cli setup-adata-region experiment1 region_001 default --plot

# 3. Write AnnData to disk
python -m spida.P.cli write-adata experiment1 --reg-name region_001 --prefix-names default
```

## Environment Variables

The CLI uses the following environment variables:

- `DEF_CUTOFFS_PATH`: Default path to filtering cutoffs JSON file
  - Default: `/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json`
- `ZARR_STORAGE_PATH`: Path to ZARR storage directory
  - Default: `/data/aklein/bican_zarr`
- `IMAGE_STORE_PATH`: Path to image storage directory
  - Default: `/ceph/cephatlas/aklein/bican/images`

## Logging

The CLI provides detailed logging information. Logs include:
- Command arguments and their types
- Processing status messages
- Error information if commands fail

## Migration from Fire CLI

If you were previously using the Fire-based CLI from `main.py`, here are the equivalent commands:

| Old Fire Command | New Argparse Command |
|------------------|---------------------|
| `python main.py filter_cells_region exp reg prefix` | `python -m spida.P.cli filter-cells-region exp reg prefix` |
| `python main.py filter_cells_all exp prefix` | `python -m spida.P.cli filter-cells-all exp prefix` |
| `python main.py write_adata exp --reg_name=reg --prefix_names=[p1,p2]` | `python -m spida.P.cli write-adata exp --reg-name reg --prefix-names p1 p2` |
| `python main.py setup_adata_region exp reg prefix` | `python -m spida.P.cli setup-adata-region exp reg prefix` |
| `python main.py setup_adata_all exp prefix` | `python -m spida.P.cli setup-adata-all exp prefix` |

Note the key differences:
- Commands use hyphens instead of underscores
- List arguments are space-separated instead of bracketed
- Boolean flags use `--flag` instead of `--flag=True`
