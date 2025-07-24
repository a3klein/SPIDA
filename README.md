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
python -m spida.{module} {command} {args}
```

## .env file

This package writes data stores and intermediate files (like segmentation specific files) into your local filesystem. In order to control where those files are written to, we utilize the .env file. In the .env file you will need to specify where you want specific forms of data stored. For the most parts, data is going to be written in the following form: 

```
{storage_location}/{experiment_name}/{region_name}/...
```

Additionally, there are times when specific binaries are needed to run specific commands, i.e. proseg, vpt (vizgen-postprocessing-tool), and deconwolf. For this there is a section in .env file where binaries need to be stored if those functions are to be used. 

Important environment variables to set: 
```
# Program Binaries
RUST_BIN_PATH={rust binary for running proseg / other rust packages}
VPT_BIN_PATH={vpt binary for converting between geometries mapping segmentations to a common structure}
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
