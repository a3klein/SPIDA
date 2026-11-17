# Quick Start Guide

This page takes you through ingesting one MERSCOPE region and running the first QC steps.

## File format

MERSCOPE data is expected in this layout:
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
|__ ...
```

`out` holds the contents of `merscope_output` and `raw` the contents of `merscope_raw_data` from a
standard Vizgen MERSCOPE run. Only `raw/dataorganization.csv` is needed from the raw data, and only
for deconvolution; without deconvolution the `raw` folder can be left out.

## Installing SPIDA

SPIDA is not on PyPI or conda; install it from git:

```bash
git clone https://github.com/a3klein/SPIDA.git
cd SPIDA
export CONDA_OVERRIDE_CUDA=12.0  # only when installing GPU environments on a machine without a GPU
pixi install -a                  # or one environment: pixi install -e preprocessing
```

## A config file

Most commands need store paths from a JSON config file. Create one:

```bash
pixi run config setup-config --config_store_path /path/to/project/config.json
```

For a first run set at least `PROCESSED_ROOT_PATH` (where the MERSCOPE output lives, as above),
`ZARR_STORAGE_PATH` (where the SpatialData stores are written) and `IMAGE_STORE_PATH` (figures).
See [Configuration](./configuration.md) for every key.

## Command structure

```bash
pixi run -e {environment} spida-{module} --config {config.json} {command} --exp {EXPERIMENT} --reg {REGION} [PREFIX] [options]
```

- `environment`: the pixi environment for the command (`preprocessing` for most)
- `module`: `S` (spatial processing), `P` (preprocessing), `I` (integration / annotation)
- `--exp` / `--reg`: the experiment and region, e.g. `experiment_1` / `region_1`
- `PREFIX`: the segmentation a command works on (e.g. `cellpose_cell`), for segmentation steps

`pixi run -e preprocessing spida-S --help` lists a module's commands; `<command> --help` lists its
options.

## First QC run

```bash
C=/path/to/project/config.json
pixi run -e preprocessing spida-S --config $C ingest-region --exp experiment_1 --reg region_1 --plot
pixi run -e preprocessing spida-P --config $C tissue-mask --exp experiment_1 --reg region_1 --plot
pixi run -e preprocessing spida-P --config $C transcript-qc --exp experiment_1 --reg region_1 --hex_size 30 --min_transcripts 100 --plot
pixi run -e preprocessing spida-P --config $C cluster-hexes --exp experiment_1 --reg region_1 --hex_size 30 --min_transcripts 100 --leiden_resolution 0.7 --plot
```

This writes the region's zarr store (`ZARR_STORAGE_PATH/experiment_1/region_1/`), its tissue mask,
and QC figures under `IMAGE_STORE_PATH/experiment_1/region_1/`.

## Next

- [The processing pipeline](./pipeline.md): deconvolution, segmentation, filtering, annotation, and
  job scripts for SLURM.
- [Project-level datasets](./project_level.md): combining processed regions.
- [Data layout](./data_layout.md): what is written where.
