# SPIDA Segmentation CLI Usage Examples

The `cli.py` script provides a command-line interface for the SPIDA segmentation pipeline using Python's `argparse` instead of `fire`.

## Installation and Setup

Make sure you have the required dependencies installed and your environment variables configured:
- `PROCESSED_ROOT_PATH`
- `SEGMENTATION_OUT_PATH` 
- `ZARR_STORAGE_PATH`

## Usage

The CLI provides three main commands:

### 1. Run Segmentation on a Single Region

```bash
python cli.py run <type> <exp_name> <reg_name> [options]
```

**Examples:**
```bash
# Run Cellpose segmentation
python cli.py run cellpose experiment_1 region_001

# Run VPT segmentation with custom config
python cli.py run vpt experiment_1 region_001 --config-path /path/to/config.json

# Run Proseg with custom input/output directories
python cli.py run proseg experiment_1 region_001 --input-dir /path/to/input --output-dir /path/to/output

# Run Mesmer segmentation with additional kwargs
python cli.py run mesmer experiment_1 region_001 --kwargs model_type=nuclear channels=1,2
```

### 2. Run Segmentation on All Regions in an Experiment

```bash
python cli.py experiment <type> <exp_name> [options]
```

**Examples:**
```bash
# Run Cellpose on all regions in an experiment
python cli.py experiment cellpose experiment_1

# Run VPT on all regions with custom directories
python cli.py experiment vpt experiment_1 --input-dir /path/to/input --output-dir /path/to/output

# Run Cellpose on all regions with additional parameters
python cli.py experiment cellpose experiment_1 --kwargs diameter=30 flow_threshold=0.4
```

### 3. Align Proseg Transcripts

```bash
python cli.py align <exp_name> <reg_name> [options]
```

**Examples:**
```bash
# Basic alignment
python cli.py align experiment_1 region_001

# Alignment with custom parameters
python cli.py align experiment_1 region_001 \
    --seed-prefix-name custom_seed \
    --prefix-name custom_proseg \
    --min-jaccard 0.5 \
    --min-prob 0.6 \
    --filter-blank

# Alignment with custom file names
python cli.py align experiment_1 region_001 \
    --cell-metadata-fname custom_metadata.csv \
    --cell-by-gene-fname custom_cell_gene.csv \
    --detected-transcripts-fname custom_transcripts.csv

# Alignment with additional custom parameters
python cli.py align experiment_1 region_001 \
    --kwargs max_distance=50 use_gpu=true custom_param=test_value
```

## Additional Keyword Arguments (--kwargs)

All commands now support the `--kwargs` parameter to pass additional keyword arguments to the downstream functions. This provides flexibility for algorithm-specific parameters without cluttering the main CLI interface.

**Format:**
```bash
--kwargs key1=value1 key2=value2 key3=value3
```

**Examples:**
```bash
# Cellpose with custom parameters
python cli.py run cellpose exp1 reg1 --kwargs diameter=25 flow_threshold=0.5 cellprob_threshold=0.0

# VPT with custom settings
python cli.py run vpt exp1 reg1 --kwargs num_workers=8 batch_size=32

# Proseg alignment with custom parameters
python cli.py align exp1 reg1 --kwargs max_iterations=100 convergence_threshold=0.001
```

## Available Segmentation Types

- `proseg`: ProSeg segmentation algorithm
- `vpt`: VPT (Voronoi Polygon Tessellation) segmentation
- `cellpose`: Cellpose segmentation algorithm
- `mesmer`: Mesmer segmentation algorithm

## Help

For detailed help on any command:

```bash
# General help
python cli.py --help

# Help for specific commands
python cli.py run --help
python cli.py experiment --help
python cli.py align --help
```

## Migration from Fire

If you were previously using the fire-based interface, here's how to migrate:

**Old (fire):**
```bash
python main.py run_segmentation --type=cellpose --exp_name=exp1 --reg_name=reg1
```

**New (argparse):**
```bash
python cli.py run cellpose exp1 reg1
```

**Old (fire):**
```bash
python main.py segment_experiment --type=vpt --exp_name=exp1 --config_path=/path/to/config.json
```

**New (argparse):**
```bash
python cli.py experiment vpt exp1 --config-path /path/to/config.json
```

The new CLI provides better help documentation, type checking, and clearer command structure while maintaining all the functionality of the original fire-based interface.
