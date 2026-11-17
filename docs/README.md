#  SPIDA (Spatial Preprocessing Imputation Distribution Analysis)

A package for MERSCOPE experiment analysis: ingestion, deconvolution, segmentation, QC,
clustering and annotation, from single regions to multi-sample projects.

## Quick Start

### 1. Install dependencies

SPIDA uses [pixi](https://pixi.sh). Clone the repository and install an environment:

```bash
git clone https://github.com/a3klein/SPIDA.git
cd SPIDA
pixi install -e preprocessing
```

Environments:
- **preprocessing**: main environment for most commands
- **preprocessing-gpu**: preprocessing with GPU support (deconvolution)
- **cellpose**: cellpose segmentation
- **mmc**: MapMyCells annotation
- **dist**: distribution analysis for downstream tasks

### 2. Configure your project

Create a JSON config file with your data paths and binary locations. See the
[Configuration Guide](./configuration.md).

### 3. Run commands

```bash
pixi run -e preprocessing spida-{module} --config config.json {command} --exp {EXPERIMENT} --reg {REGION} [args]
```

## Documentation

- **[Quick Start Guide](./quick_start.md)**: ingest a region and run the first QC
- **[Configuration Guide](./configuration.md)**: the config file and its keys
- **[Processing pipeline](./pipeline.md)**: the five steps from raw data to annotated tables, and SLURM job scripts
- **[Project-level datasets](./project_level.md)**: combining regions into one dataset
- **[Data layout](./data_layout.md)**: what is written where, and the `spida.io` layer
- **[Software setup](./software_setup.md)**: proseg, deconwolf and the optional VPT binary
- **[Deconvolution setup](./deconvolution.md)**: configure deconwolf for your imaging setup

## Module overview

### S (Spatial processing)

- Build SpatialData (`.zarr`) stores from MERSCOPE output
- Run segmentations (cellpose, proseg, user-provided) and load them into the store
- Deconvolve images with deconwolf

`pixi run -e preprocessing spida-S --help`

### P (Preprocessing)

- Tissue masks, transcript-level QC and hex clustering
- Cell filtering and QC metrics
- Normalisation, dimensionality reduction and clustering, per region and per project
- Project-level datasets (`combine-datasets`, `setup-dataset`, `sync-project`)

`pixi run -e preprocessing spida-P --help`

### I (Integration / annotation)

- Annotate tables (MapMyCells, ALLCools)
- Sync annotation columns back into the zarr store

`pixi run -e mmc spida-mmc --help`

## Project structure

```
SPIDA/
├── src/spida/
│   ├── config.py      # configuration
│   ├── io/            # all reads and writes of SPIDA data (layout, zarr store, AnnData)
│   ├── S/             # spatial processing
│   ├── P/             # preprocessing
│   ├── I/             # integration / annotation
│   ├── pl/            # plotting
│   ├── utilities/     # shared helpers
│   └── site.py        # QC figures for the SPIDA website
├── slurm_scripts/     # Anvil and AWS job-script generators and templates
├── docs/              # this documentation
└── pixi.toml          # environments and tasks
```
