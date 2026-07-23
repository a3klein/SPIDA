#  SPIDA (Spatial Preprocessing Imputation Distrubition Analysis)

A package to host MERSCOPE experiment analysis. 


# Usage: 

## pixi

This package uses pixi to manage dependencies, which means that as of now it is composed of multiple features with conflicting dependencies. Therefore to install a specific version of this package you need to ..... 

To run commands using the cli you need to first clone this repository, and then from the root of this repo run `pixi install -e {environment}`. 
"preprocessing" is the major environment that can handle most of the functions of this package, however some specific features like cellpose / deepcell (mesmer) segmentations need to run from their own environments. Additionally, mapmycells also has its own environment due to conflicting dependencies. 

After installing a speicifc environment (stored in .pixi/envs), you can run 
```
pixi run -e {environment} {command}. 
```

The command is composed of multiple parts: 
```
python -m spida.{module}.cli {command} {args}
```

## .env / config file

This package writes data stores and intermediate files (like segmentation specific files) into your local filesystem. In order to control where those files are written to, we use a configuration of paths and binary locations. This can be supplied either as a `.env` file or as a per-project config file (`.json`) — SPIDA is increasingly moving toward the config-file approach, but `.env` is still supported. Either way you specify the same keys (shown below); see the [Configuration guide](./docs/configuration.md) for the config-file format. For the most part, data is going to be written in the following form: 

```
{storage_location}/{experiment_name}/{region_name}/...
```

Additionally, some commands need external binaries: proseg (segmentation) and deconwolf (deconvolution). Point your config (`.env` / `.json`) at them if you use those steps. The vpt (vizgen-postprocessing-tool) binary is **optional** — segmentation post-processing is reimplemented in pure Python and runs by default (`--backend native`); vpt is only kept as a backwards-compatibility fallback (`--backend vpt`) and is slower than native. Outputs match up to small library-version differences (e.g. GDAL rasterization affects `sum_signals` by ~0.5–1%). 

Important environment variables to set: 
```
# Program Binaries
RUST_BIN_PATH={rust binary for running proseg / other rust packages}
VPT_BIN_PATH={OPTIONAL vpt binary; only needed for the legacy --backend vpt fallback}
DECONWOLF_CONFIG={where the deconwolf .config.ini file lives (see below)}

# Data Storage Path
ZARR_STORAGE_PATH={where the spatialdata .zarr files live}
PROCESSED_ROOT_PATH={where the raw MERSCOPE data lives}
SEGMENTATION_OUT_PATH={where intermediate files used for segmentation are stored}
ANNDATA_STORE_PATH={where adata objects are stored}
ANNOTATIONS_STORE_PATH={where intermediate files used for annotations are stored}

# Image Outputs
IMAGE_STORE_PATH={where generated figures will be stored}
```

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
```

generated psfs are stored as follows: {wavelength}_z{z-step-size}_psf.tif
