# SPIDA I Module CLI Usage Guide

The SPIDA I (Integration) module provides a command-line interface for integration and annotation tasks in the BICAN project. This CLI tool handles AnnData object backup, ALLCools integration, MapMyCells annotation, and MOSCOT integration workflows.

## Installation and Setup

Before using the CLI, ensure you have the SPIDA package installed and the necessary environment variables set for different integration environments:

- For ALLCools integration: Set up the "preprocessing" environment
- For MapMyCells annotation: Ensure proper environment configuration

## General Usage

```bash
python -m spida.I.cli [command] [arguments] [options]
```

### Get Version Information

```bash
python -m spida.I.cli --version
```

### Get Help

```bash
python -m spida.I.cli --help
python -m spida.I.cli [command] --help
```

## Available Commands

The I module is organized into several functional categories:

### 1. Utilities

#### Backup AnnData for a Specific Region

Backup AnnData objects for a specific region due to environment incompatibilities.

```bash
python -m spida.I.cli backup-adata-region EXP_NAME REG_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--adata-path PATH`: Path to the store of AnnData objects (default: None)

**Example:**
```bash
python -m spida.I.cli backup-adata-region experiment1 region_001 default
```

#### Backup AnnData for an Entire Experiment

Backup AnnData objects for all regions in an experiment.

```bash
python -m spida.I.cli backup-adata-experiment EXP_NAME PREFIX_NAME [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object

**Options:**
- `--adata-path PATH`: Path to the store of AnnData objects (default: None)

**Example:**
```bash
python -m spida.I.cli backup-adata-experiment experiment1 default
```

### 2. ALLCools Integration

ALLCools integration requires the "preprocessing" environment.

#### ALLCools Integration for a Specific Region

Run ALLCools integration on a given experiment and region.

```bash
python -m spida.I.cli allcools-integration-region EXP_NAME REG_NAME PREFIX_NAME REF_PATH [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object
- `REF_PATH`: Path to the reference RNA AnnData object

**Options:**
- `--anndata-store-path PATH`: Path to the store of AnnData objects (default: None)
- `--annotations-store-path PATH`: Path to the store of annotation specific files (default: None)
- `-k`, `--kwargs KEY=VALUE`: Additional keyword arguments for ALLCools integration

**Example:**
```bash
python -m spida.I.cli allcools-integration-region experiment1 region_001 default /path/to/reference.h5ad \
    --anndata-store-path /data/anndatas \
    -k resolution=0.5 n_neighbors=15
```

#### ALLCools Integration for an Entire Experiment

Run ALLCools integration for all regions in an experiment.

```bash
python -m spida.I.cli allcools-integration-experiment EXP_NAME PREFIX_NAME REF_PATH [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object
- `REF_PATH`: Path to the reference RNA AnnData object

**Options:**
- `--anndata-store-path PATH`: Path to the store of AnnData objects (default: None)
- `--annotations-store-path PATH`: Path to the store of annotation specific files (default: None)
- `-k`, `--kwargs KEY=VALUE`: Additional keyword arguments for ALLCools integration

**Example:**
```bash
python -m spida.I.cli allcools-integration-experiment experiment1 default /path/to/reference.h5ad \
    --anndata-store-path /data/anndatas \
    --annotations-store-path /data/annotations
```

### 3. MapMyCells Annotation

#### MapMyCells Setup

Setup function for MapMyCells integration. This prepares the reference data and markers for annotation.

```bash
python -m spida.I.cli mmc-setup REF_PATH HIERARCHY_LIST BRAIN_REGION CODEBOOK [OPTIONS]
```

**Arguments:**
- `REF_PATH`: Path to the reference data
- `HIERARCHY_LIST`: List of hierarchy levels for the annotation (comma-separated)
- `BRAIN_REGION`: Brain region for the annotation
- `CODEBOOK`: Codebook for the annotation

**Options:**
- `--codebook-path PATH`: Path to the codebook file (default: None)
- `--mmc-store-path PATH`: Path to the store of MapMyCells markers (default: None)
- `--ref-norm METHOD`: Normalization method for the reference AnnData.X (default: log2CPM)
- `-k`, `--kwargs KEY=VALUE`: Additional keyword arguments

**Example:**
```bash
python -m spida.I.cli mmc-setup /path/to/reference.h5ad "class,subclass,type" "Whole brain" "mouse_brain_atlas" \
    --ref-norm log2CPM \
    --mmc-store-path /data/mmc_markers
```

#### MapMyCells Annotation for a Specific Region

Run MapMyCells annotation on a given experiment and region.

```bash
python -m spida.I.cli mmc-annotation-region EXP_NAME REG_NAME PREFIX_NAME BRAIN_REGION CODEBOOK [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `REG_NAME`: Name of the region
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object
- `BRAIN_REGION`: Brain region for the annotation
- `CODEBOOK`: Codebook for the annotation

**Options:**
- `--mmc-store-path PATH`: Path to the store of MapMyCells markers (default: None)
- `--anndata-store-path PATH`: Path to the store of AnnData objects (default: None)
- `--annotations-store-path PATH`: Path to the store of annotation specific files (default: None)
- `-k`, `--kwargs KEY=VALUE`: Additional keyword arguments

**Example:**
```bash
python -m spida.I.cli mmc-annotation-region experiment1 region_001 default "Whole brain" "mouse_brain_atlas" \
    --mmc-store-path /data/mmc_markers \
    --anndata-store-path /data/anndatas \
    --annotations-store-path /data/annotations
```

#### MapMyCells Annotation for an Entire Experiment

Run MapMyCells annotation for all regions in an experiment.

```bash
python -m spida.I.cli mmc-annotation-experiment EXP_NAME PREFIX_NAME BRAIN_REGION CODEBOOK [OPTIONS]
```

**Arguments:**
- `EXP_NAME`: Name of the experiment
- `PREFIX_NAME`: Prefix for the keys in the spatialdata object
- `BRAIN_REGION`: Brain region for the annotation
- `CODEBOOK`: Codebook for the annotation

**Options:**
- `--mmc-store-path PATH`: Path to the store of MapMyCells markers (default: None)
- `--anndata-store-path PATH`: Path to the store of AnnData objects (default: None)
- `--annotations-store-path PATH`: Path to the store of annotation specific files (default: None)
- `-k`, `--kwargs KEY=VALUE`: Additional keyword arguments

**Example:**
```bash
python -m spida.I.cli mmc-annotation-experiment experiment1 default "Whole brain" "mouse_brain_atlas" \
    --mmc-store-path /data/mmc_markers \
    --anndata-store-path /data/anndatas \
    --annotations-store-path /data/annotations
```

### 4. MOSCOT Integration

#### MOSCOT Integration (Placeholder)

This command is a placeholder for MOSCOT integration functionality that is not yet implemented.

```bash
python -m spida.I.cli moscot-integration
```

**Note:** This functionality is currently under development.

## Typical Workflows

### Complete Integration and Annotation Workflow

Here's a typical workflow for processing an experiment with both ALLCools and MapMyCells:

```bash
# 1. Backup AnnData objects (if needed for environment compatibility)
python -m spida.I.cli backup-adata-experiment experiment1 default

# 2. Setup MapMyCells markers
python -m spida.I.cli mmc-setup /path/to/reference.h5ad "class,subclass,type" "Whole brain" "mouse_brain_atlas" \
    --mmc-store-path /data/mmc_markers

# 3. Run ALLCools integration for the experiment
python -m spida.I.cli allcools-integration-experiment experiment1 default /path/to/reference.h5ad \
    --anndata-store-path /data/anndatas \
    --annotations-store-path /data/annotations

# 4. Run MapMyCells annotation for the experiment
python -m spida.I.cli mmc-annotation-experiment experiment1 default "Whole brain" "mouse_brain_atlas" \
    --mmc-store-path /data/mmc_markers \
    --anndata-store-path /data/anndatas \
    --annotations-store-path /data/annotations
```

### Single Region Processing

For processing a single region:

```bash
# 1. Backup AnnData for the region
python -m spida.I.cli backup-adata-region experiment1 region_001 default

# 2. Run ALLCools integration for the region
python -m spida.I.cli allcools-integration-region experiment1 region_001 default /path/to/reference.h5ad

# 3. Run MapMyCells annotation for the region
python -m spida.I.cli mmc-annotation-region experiment1 region_001 default "Whole brain" "mouse_brain_atlas"
```

## Environment Requirements

### ALLCools Integration
- Requires the "preprocessing" environment to be activated
- Ensure ALLCools and related dependencies are installed

### MapMyCells Annotation
- Requires proper MapMyCells environment setup
- Ensure access to reference data and codebooks

## Command Line Arguments

### Common Argument Types

- **Paths**: Use absolute paths for reliability
- **Lists**: For hierarchy lists, use comma-separated values (e.g., "class,subclass,type")
- **Key-Value Pairs**: Use `-k key=value` format for additional parameters

### Path Arguments
The CLI accepts various path arguments that can be:
- Absolute paths: `/full/path/to/file`
- Relative paths: `./relative/path`
- Environment variables will be used for defaults when available

## Logging

The CLI provides detailed logging information including:
- Command arguments and their types
- Processing status for each step
- Error messages with context
- Completion status

## Error Handling

Common issues and solutions:

1. **Environment Incompatibilities**: Use backup commands to save AnnData objects before switching environments
2. **Missing Reference Data**: Ensure reference paths are correct and files exist
3. **Permission Issues**: Check write permissions for output directories
4. **Memory Issues**: Consider processing regions individually for large experiments

## Integration with Other Modules

The I module works in conjunction with:
- **P Module**: Provides filtered and setup AnnData objects for integration
- **Segmentation Modules**: Provides cell boundaries and metadata
- **Visualization Modules**: Can plot results from integration and annotation

This CLI provides a comprehensive interface for the integration and annotation pipeline in the BICAN project, allowing for flexible processing of both individual regions and entire experiments.
