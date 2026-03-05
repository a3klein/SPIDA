# Quick Start Guide

This page will take you through setting up the basic QC on a merscope experiment.

## File Format

MRESCOPE data is assumed to be in the following folder structure: 
```
PROCESSED_ROOT_PATH
|__experiment_1
  |__out
    |__region_1
    |__region_2
    |__ ...
  |__raw
|__experiment_2
  |__out
  |__raw
|__experiment_3
|__ ... 
```

Where the `out` folder is the contents of `merscope_output` and the `raw` folder is the contents of `merscope_raw_data` from the standard vizgen MERSCOPE run Z drive. Technically for a simple preprocessing pipeline only the `dataorganization.csv` is needed from the raw data, and only for the deconvolution step. So, if you want to skip deconvolution you can ignore the raw folder. If you want use deconvolution but not have to worry about the raw folder, you can copy the `dataorganization.csv` file into the `out` folder and ignore the raw data. 

## Cloning Spida

SPIDA is currently not available as a pypi or conda package, you need to install it using git. 

```
git clone https://github.com/a3klein/SPIDA.git
export CONDA_OVERRIDE_CUDA=12.0 # This is necessary if attempting to install a package with dependency on a GPU from an environment with no GPU available (i.e. from a different slurm / aws instance).
cd SPIDA
pixi install -a
```

This will make sure you all the spida environments installed within the SPIDA pixi environment.

## configuring an empty config file

in order to run most steps there needs to at least exist a config file. SPIDA takes as an argument a --config option (.json or .env file) that holds within it default parameters such as paths to binaries, data paths, etc... . If the config file is not passed, SPIDA looks for a .env file in the SPIDA directory. I suggest at this point populating an empty .env file in the SPIDA directory and then adding to it later once you get a better understanding of how the package works. 

For now, run the following command: 
```
pixi run config setup-config
```

This will guide you through instructions on what each parameter in the config file represents, for now I will only setup the `PROCESSED_ROOT_PATH` parameter (this is where your MERSCOPE outputs live as described above), and the `ZARR_STORAGE_PATH` parameter (this is where the spatialdata output objects are going to live). It will also be nice to populate the `IMAGE_STORE_PATH` at this point, as this will allow you to automatically write out images from these initial steps. 

See the [Configurations](./configuration.md) page for more information on this. 

## Running test commands

The standard strucure for a SPIDA command is: 

```
pixi run -e {environment} SPIDA-{module} {command} {EXPERIMENT_NAME} {REGION_NAME} {PREFIX}
```

- `environment` specifies which pixi environment is used for a given command
- `module` specifies which module to run (S, P, I, ...)
- `command` specifies the precise command within that module
- `EXPERIMENT_NAME` is the name of the experiment for which to run the command (i.e. experiment_{1/2/3/...} in this example)
- `REGION_NAME` is the name of the region in the given experiment for which to run the command (i.e. region_{1/2/...} in this example)
- `PREFIX` (Optional) is the name of the segmentation on which to run that command (this only becomes relevant in segmentation dependent commands)


## Running initial QC runs
