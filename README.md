#  SPIDA (Spatial Preprocessing Imputation Distrubition Analysis)

A package to host MERSCOPE experiment analysis. 


# Usage: 

## pixi

This package uses pixi to manage dependencies, which means that as of now it is composed of multiple features with conflicting dependencies, each installed as its own environment.

To run commands using the cli you need to first clone this repository, and then from the root of this repo run `pixi install -e {environment}`. 
"preprocessing" is the major environment that can handle most of the functions of this package, however cellpose segmentation runs from its own environment (`cellpose`), as does MapMyCells (`mmc`), due to conflicting dependencies. 

After installing a specific environment (stored in .pixi/envs), you can run 
```
pixi run -e {environment} {command}
```

The command is composed of multiple parts: 
```
python -m spida.{module}.cli --config {config.json} {command} --exp {experiment} --reg {region} {args}
```
(`spida-S`, `spida-P` and `spida-I` are pixi task shortcuts for `python -m spida.{module}`.)

## Config file

SPIDA writes data stores and intermediate files (segmentation outputs, h5ad tables, figures) to
your filesystem. Where they go, and where external binaries live, is set in a per-project JSON
config file passed with `--config` (environment variables of the same names are used for keys the
file does not set; `.env` files are no longer read). See the
[Configuration guide](./docs/configuration.md) for every key. Data is written as:

```
{store}/{experiment_name}/{region_name}/...
```

and multi-region projects under `{store}/_projects/{project_name}/` (see
[Project-level datasets](./docs/project_level.md)).

Some commands need external binaries: proseg (segmentation) and deconwolf (deconvolution). Point
your config at them if you use those steps. The vpt (vizgen-postprocessing-tool) binary is
**optional** — segmentation post-processing is reimplemented in pure Python and runs by default
(`--backend native`); vpt is only kept as a backwards-compatibility fallback (`--backend vpt`) and
is slower than native. Outputs match up to small library-version differences (e.g. GDAL
rasterization affects `sum_signals` by ~0.5–1%).

Main keys:
```
# Program binaries
RUST_BIN_PATH          directory with the proseg binary
DECONWOLF_CONFIG       deconwolf .config.ini (see below)
VPT_BIN_PATH           OPTIONAL; only for the legacy --backend vpt fallback

# Data stores
PROCESSED_ROOT_PATH    raw MERSCOPE output
ZARR_STORAGE_PATH      SpatialData .zarr stores
SEGMENTATION_OUT_PATH  segmentation outputs
ANNDATA_STORE_PATH     h5ad tables
ANNOTATION_STORE_PATH  annotation intermediates
IMAGE_STORE_PATH       figures
CUTOFFS_PATH           QC cutoffs JSON
```

## Documentation

- [Quick start](./docs/quick_start.md) · [Configuration](./docs/configuration.md) ·
  [Processing pipeline](./docs/pipeline.md) · [Project-level datasets](./docs/project_level.md) ·
  [Data layout](./docs/data_layout.md) · [Software setup](./docs/software_setup.md) ·
  [Deconvolution](./docs/deconvolution.md)

# Modules: 

SPIDA is composed of several modules meant to handle different tasks: 

## S (spatial processing): 

The S module is built to: 
- Setup the spatialdata objects (.zarr stores)
- Run any number of segmentation pipelines (cellpose / proseg / ...)
- Deconvolve images using deconwolf (see decon section below)

## P (Preprocessing): 

The P module is built to: 
- Filter cells and generate QC metrics
- Extract anndata objects into separate stores for easier downstream analysis 
- Calculate dimensionality reduction and clustering through a number of separate spatially aware methods
- Calculate and resolve doublets

## I (Integration / Imputation): 

The I module is built to: 
- annotate data using methods like MapMyCells, ALLCools, MOSCOT, ...
- Impute whole transcriptome RNA expression
- Impute methylation data onto the spatial context. 



# Deconvolution (deconwolf): 

To run deconvolution (highly recommended for high confidence nuclear segmentations), deconwolf needs to be installed on your machine. Then dw_bw needs to be used to generate PSFs for the specific imaging setup you have (the NA number, xy / z resolution, imaging light wavelength, etc...). Then a .config.ini file needs to be generated with the following structure: 

```
[Paths]
dw_path = {Path to dw binary}
dw_psf_path = {path to directory containing your generated psf files}

[Optics]
na = 1.4
ni = 1.518
pixel_size_nm = 106.0
z_step_nm = 1500
emission_nm_DAPI = 405
emission_nm_PolyT = 488
```

generated psfs are stored as follows: {wavelength}_z{z-step-size}_psf.tif

`[Optics]` is optional and holds the imaging parameters for the setup those PSFs were
built for, so filters needing NA / pixel size / wavelength read one value instead of
hardcoding their own. Missing keys fall back to `OPTICS_DEFAULTS` in
`spida/_constants.py`. See `docs/deconvolution.md`.
