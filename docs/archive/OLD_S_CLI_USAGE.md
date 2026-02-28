# SPIDA S module CLI Usage Examples

The `cli.py` script provides a command-line interface for the SPIDA spatial preprocessing pipeline.

## Description

Command line interface for SPIDA S (Spatial Processing) Module. This interface provides methods to ingest raw data, segment spatial data using a variety of algorithms, deconvolve large images for better segmentation, and align cell based segmentations to the nuclei segmentation they are rooted in to remove dubious artifacts.

## Installation and Setup

Make sure you have the required dependencies installed and your environment variables configured:
- `PROCESSED_ROOT_PATH`
- `SEGMENTATION_OUT_PATH` 
- `ZARR_STORAGE_PATH`

## Available Commands

The CLI provides the following main command categories:

### [Utilities]
- `decon-image` - Deconvolve large image files in tiles using DeconWolf algorithm

### [Segmentations] 
- `run` / `run-segmentation-region` - Run segmentation on a single region using specified algorithm
- `experiment` / `segment-experiment` - Run segmentation for all regions in an experiment using specified algorithm

### [Alignment]    
- `align` / `align-proseg` - Align cell-based segmentations to the nuclei segmentation

### [I/O]
- `ingest-region` - Ingest a specific region of an experiment into a spatialdata object
- `ingest-all` - Ingest all regions of an experiment into spatialdata objects
- `load-segmentation-region` - Load segmentation data for a specific region into a spatialdata object
- `load-segmentation-all` - Load segmentation data for all regions of an experiment into spatialdata objects

## Command Details

### 1. Deconvolve Large Images

```bash
python cli.py decon-image --image_path <path> --data_org_path <path> --channels <channels> [options]
```

**Examples:**

```bash
# Basic deconvolution with DAPI channel
python cli.py decon-image --image_path /path/to/image --data_org_path /path/to/data_org.txt --channels DAPI

# Deconvolution with multiple channels and custom tile size
python cli.py decon-image --image_path /path/to/image --data_org_path /path/to/data_org.txt --channels PolyT,DAPI --tile_size 2960 --overlap 100

# Deconvolution with visualization and GPU acceleration
python cli.py decon-image --image_path /path/to/image --data_org_path /path/to/data_org.txt --channels DAPI --visualize_grid --gpu true
```

**Options:**
- `--output_dir`: Output directory for tiles (default: tiles_output)
- `--tile_size`: Tile size in pixels (default: 2960)
- `--overlap`: Overlap between tiles in pixels (default: 100)
- `--visualize_grid`: Visualize the tiling grid
- `--z_step`: Axial(z) step size in micrometers (default: 1.5)
- `--filter`: Filter to apply to the image before segmentation
- `--filter_args`: Additional filter arguments
- `--gpu`: Use GPU (default: false)
- `--continue_stalled`: Continue processing if some tiles already processed (default: false)
- `--plot_thr`: Plot thresholding histogram (default: false)

### 2. Run Segmentation on a Single Region

```bash
python cli.py run <type> <exp_name> <reg_name> [options]
```

**Examples:**

```bash
# Run Cellpose segmentation
python cli.py run cellpose experiment_1 region_001

# Run VPT segmentation with custom config
python cli.py run vpt experiment_1 region_001 --config_path /path/to/config.json

# Run Proseg with custom input/output directories
python cli.py run proseg experiment_1 region_001 --input_dir /path/to/input --output_dir /path/to/output

# Run Mesmer segmentation with additional kwargs
python cli.py run mesmer experiment_1 region_001 --kwargs model_type=nuclear channels=1,2
```

**Available segmentation types:**
- `proseg`: ProSeg segmentation algorithm
- `vpt`: VPT (Vizgen Postprocessing Tool) segmentation
- `cellpose`: Cellpose segmentation algorithm
- `mesmer`: Mesmer segmentation algorithm

**Options:**
- `--input_dir`: Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)
- `--output_dir`: Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)
- `--config_path`: Configuration file path (for VPT segmentation)
- `--kwargs`: Additional keyword arguments to segmentation algorithms in key=value format

### 3. Run Segmentation on All Regions in an Experiment

```bash
python cli.py experiment <type> <exp_name> [options]
```

**Examples:**

```bash
# Run Cellpose on all regions in an experiment
python cli.py experiment cellpose experiment_1

# Run VPT on all regions with custom directories
python cli.py experiment vpt experiment_1 --input_dir /path/to/input --output_dir /path/to/output

# Run Cellpose on all regions with additional parameters
python cli.py experiment cellpose experiment_1 --kwargs diameter=30 flow_threshold=0.4
```

**Options:**
- `--input_dir`: Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)
- `--output_dir`: Directory to save the output data (default: uses SEGMENTATION_OUT_PATH env var)
- `--config_path`: Configuration file path (for VPT segmentation)
- `--kwargs`: Additional keyword arguments in key=value format

### 4. Align Proseg Transcripts

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

**Options:**
- `--seed-prefix-name`: Seed prefix name (default: default)
- `--prefix-name`: Prefix name (default: proseg)
- `--out-prefix-name`: Output prefix name (default: proseg_aligned)
- `--input-dir`: Directory containing the input data (default: uses PROCESSED_ROOT_PATH env var)
- `--seg-dir`: Segmentation directory (default: uses SEGMENTATION_OUT_PATH env var)
- `--x`: X coordinate column name (default: x)
- `--y`: Y coordinate column name (default: y)
- `--z`: Z coordinate column name (default: global_z)
- `--cell-column`: Cell column name (default: cell_id)
- `--barcode-column`: Barcode column name (default: barcode_id)
- `--gene-column`: Gene column name (default: gene)
- `--fov-column`: FOV column name (default: fov)
- `--cell-missing`: Cell missing value (default: -1)
- `--min-jaccard`: Minimum Jaccard index (default: 0.4)
- `--min-prob`: Minimum probability (default: 0.5)
- `--filter-blank`: Filter blank genes
- `--cell-metadata-fname`: Cell metadata filename (default: merged_cell_metadata.csv)
- `--cell-by-gene-fname`: Cell by gene filename (default: merged_cell_by_gene.csv)
- `--detected-transcripts-fname`: Detected transcripts filename (default: merged_transcript_metadata.csv)
- `--cell-polygons-fname`: Cell polygons filename (default: merged_cell_polygons.geojson)
- `--kwargs`: Additional keyword arguments in key=value format

### 5. Ingest Region Data

```bash
python cli.py ingest-region <exp_name> <reg_name> [options]
```

**Examples:**

```bash
# Basic ingestion
python cli.py ingest-region experiment_1 region_001

# Ingestion with custom type and prefix
python cli.py ingest-region experiment_1 region_001 --type merscope --prefix-name custom_prefix

# Ingestion with plotting
python cli.py ingest-region experiment_1 region_001 --plot
```

**Options:**
- `--type`: Type of the data to ingest (default: merscope)
- `--prefix-name`: Prefix for the keys in the spatialdata object (default: default)
- `--source`: Source of the data (default: machine)
- `--plot`: Plot results after ingestion

### 6. Ingest All Regions

```bash
python cli.py ingest-all <exp_name> [options]
```

**Examples:**

```bash
# Basic ingestion of all regions
python cli.py ingest-all experiment_1

# Ingestion with custom settings
python cli.py ingest-all experiment_1 --type merscope --prefix-name batch_001 --plot
```

**Options:**
- `--type`: Type of the data to ingest (default: merscope)
- `--prefix-name`: Prefix for the keys in the spatialdata object (default: default)
- `--source`: Source of the data (default: machine)
- `--plot`: Plot results after ingestion

### 7. Load Segmentation Data for Region

```bash
python cli.py load-segmentation-region <exp_name> <reg_name> <seg_dir> [options]
```

**Examples:**

```bash
# Basic loading
python cli.py load-segmentation-region experiment_1 region_001 /path/to/segmentation

# Loading with custom type and plotting
python cli.py load-segmentation-region experiment_1 region_001 /path/to/segmentation --type vpt --plot
```

**Options:**
- `--type`: Type of the segmentation data to load (default: vpt)
- `--prefix-name`: Prefix for the keys in the spatialdata object (default: default)
- `--plot`: Plot results after loading segmentation
- `--load_kwargs`: Additional keyword arguments for loading segmentation data

### 8. Load Segmentation Data for All Regions

```bash
python cli.py load-segmentation-all <exp_name> <seg_dir> [options]
```

**Examples:**

```bash
# Basic loading for all regions
python cli.py load-segmentation-all experiment_1 /path/to/segmentation

# Loading with custom settings
python cli.py load-segmentation-all experiment_1 /path/to/segmentation --type vpt --prefix-name batch_seg --plot
```

**Options:**
- `--type`: Type of the segmentation data to load (default: vpt)
- `--prefix-name`: Prefix for the keys in the spatialdata object (default: default)
- `--plot`: Plot results after loading segmentation
- `--load_kwargs`: Additional keyword arguments for loading segmentation data

## Additional Keyword Arguments (--kwargs)

All segmentation and alignment commands support the `--kwargs` parameter to pass additional keyword arguments to the downstream functions. This provides flexibility for algorithm-specific parameters without cluttering the main CLI interface.

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

## Version Information

```bash
# Check version
python cli.py -v
# or
python cli.py --version
```

## Help

For detailed help on any command:

```bash
# General help
python cli.py --help

# Help for specific commands
python cli.py run --help
python cli.py experiment --help
python cli.py align --help
python cli.py decon-image --help
python cli.py ingest-region --help
python cli.py ingest-all --help
python cli.py load-segmentation-region --help
python cli.py load-segmentation-all --help
```

## Environment Variables

The CLI uses the following environment variables as defaults:

- `PROCESSED_ROOT_PATH`: Default input directory for processed data
- `SEGMENTATION_OUT_PATH`: Default output directory for segmentation results
- `ZARR_STORAGE_PATH`: Path for zarr storage

## Command Aliases

Several commands have aliases for convenience:

- `run` can also be called as `run-segmentation-region`
- `experiment` can also be called as `segment-experiment`
- `align` can also be called as `align-proseg`

## Author

Amit Klein

---

*This documentation reflects the current implementation of the SPIDA S module CLI interface.*
