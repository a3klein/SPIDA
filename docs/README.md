#  SPIDA (Spatial Preprocessing Imputation Distribution Analysis)

A package to host MERSCOPE experiment analysis with modular preprocessing, analysis, and integration pipelines.

## Quick Start

### 1. Install Dependencies

This package uses [pixi](https://pixi.sh) for dependency management. Clone the repository and install your desired environment:

```bash
git clone <repository-url>
cd SPIDA
pixi install -e preprocessing
```

Available environments:
- **preprocessing**: Main environment for most functionality
- **preprocessing-gpu**: preprocessing but with GPU support
- **cellpose**: Cellpose segmentation
- **deepcell**: DeepCell/Mesmer segmentation
- **mmc**: MapMyCells environment for annotations
- **dist**: distribution analysis for downstream tasks

### 2. Configure Your Project

Create a `.env` file in the repository root or a `.json` file to store project configurations. See [Quick Start Guide](./docs/QUICK_START.md) for details.

### 3. Run Commands

```bash
pixi run -e preprocessing python -m spida.{module}.cli {command} {args}
```
OR 
```bash
pixi run -e preprocessing spida-{module} {command} {args}
```

## Documentation

- **[Quick Start Guide](./docs/QUICK_START.md)**: Get running in 5 minutes
- **[Configuration Guide](./docs/CONFIGURATION.md)**: Setting up a config file for your project
- **[S Module (Spatial Processing)](./docs/S_USAGE.md)**: IO, deconvolution, segmentation, ...
- **[P Module (Preprocessing)](./docs/P_USAGE.md)**: cell level QC + clustering, transcript level analysis
- **[I Module (Integration)](./docs/I_USAGE.md)**: Annotation, imputation, and integration
- **[D Module (Distribution)](./docs/D_USAGE.md)**: Distribution Analysis (coming soon...)
- **[Deconvolution Setup](./docs/DECONVOLUTION.md)**: Configure deconwolf for your imaging setup

## Module Overview

### S (Spatial Processing)

Handles image processing and segmentation:
- Setup spatialdata objects (`.zarr` stores)
- Run segmentation pipelines (cellpose, proseg, etc.)
- Deconvolve images using deconwolf

**Quick start**: `pixi run -e preprocessing spida-S --help`

### P (Preprocessing)

Handles cell and transcript analysis:
- Filter cells and generate QC metrics
- Extract anndata objects to separate stores
- Calculate dimensionality reduction and clustering
- Detect and resolve doublets

**Quick start**: `pixi run -e preprocessing spida-P --help`

### I (Integration/Imputation)

Handles data annotation and imputation:
- Annotate data (MapMyCells, ALLCools, etc.)
- Impute whole transcriptome RNA expression
- Impute methylation data onto spatial context

**Quick start**: `pixi run -e preprocessing spida-I --help`

## Project Structure

```
SPIDA/
├── src/spida/
│   ├── config.py              # Central configuration system
│   ├── settings.py            # Global settings management
│   ├── S/                      # Spatial processing module
│   ├── P/                      # Preprocessing module
│   └── I/                      # Integration module
|   └── D/                      # Distribution module
├── docs/                       # Documentation files
├── .env                        # Environment variables (create this)
├── pixi.toml                   # Pixi configuration
└── pyproject.toml              # Project metadata
```

## Advanced Topics

- [Configuration Reference](./docs/CONFIG_REFERENCE.md): All available configuration options
- [Deconvolution Setup](./docs/DECONVOLUTION.md): Configure deconwolf for your imaging setup
- [Architecture Overview](./docs/ARCHITECTURE.md): Understand the system design

## Support

For issues or questions, see the relevant module documentation or check existing issues.