# SPIDA Configuration Guide

SPIDA reads its paths — data stores, binaries, cutoff files — from a JSON config file passed
with `--config` before the command:

```bash
pixi run -e preprocessing spida-P --config /path/to/project/config.json <command> ...
```

## Resolution order

1. Command-line options (e.g. `--zarr_store`, `--zarr_path`)
2. The `--config` JSON
3. Environment variables of the same name
4. Unset — a command that needs a missing path stops and names the key

`.env` files are no longer read; a config path that is not `.json` is rejected.

## Keys

| Key | Meaning |
|---|---|
| `PROCESSED_ROOT_PATH` | Raw MERSCOPE output: `{root}/{experiment}/out/{region}/` and `{root}/{experiment}/raw/`. |
| `ZARR_STORAGE_PATH` | SpatialData zarr stores. |
| `SEGMENTATION_OUT_PATH` | Segmenter output and segmentation-schema files. |
| `ANNDATA_STORE_PATH` | h5ad copies of tables, and project tables. |
| `ANNOTATION_STORE_PATH` | Annotation intermediates (`ANNOTATIONS_STORE_PATH` is accepted as an alias). |
| `IMAGE_STORE_PATH` | Figures. |
| `CUTOFFS_PATH` | QC cutoffs JSON for `filter-cells-region`. |
| `REGION_LABELS` | Optional JSON `{region name: label}` (e.g. donor); default is the region name without `region_`. |
| `RUST_BIN_PATH` | Directory with the `proseg` binary. |
| `DECONWOLF_CONFIG` | deconwolf `.config.ini` (see [Deconvolution](./deconvolution.md)). |
| `MMC_DIR` | MapMyCells reference store. |
| `CELLPOSE_MODEL_PATH` | Optional cellpose weights directory. |
| `GENE_PANEL_PATH` | Optional gene panel files. |
| `SPIDA_SITE_DIR` | Optional SPIDA website directory (QC figures). |
| `VPT_BIN_PATH` | **Optional** — only for the legacy `--backend vpt` post-processing fallback. |

See [Data layout](./data_layout.md) for what each store holds.

Example:

```json
{
    "PROCESSED_ROOT_PATH": "/data/merscope",
    "ZARR_STORAGE_PATH": "/project/data/zarr_store",
    "SEGMENTATION_OUT_PATH": "/project/data/segmented",
    "ANNDATA_STORE_PATH": "/project/data/aggregated",
    "ANNOTATION_STORE_PATH": "/project/data/annotated",
    "IMAGE_STORE_PATH": "/project/images",
    "CUTOFFS_PATH": "/project/config/filtering_cutoffs.json",
    "REGION_LABELS": "/project/config/region_labels.json",
    "RUST_BIN_PATH": "/home/user/.cargo/bin",
    "DECONWOLF_CONFIG": "/project/config/deconwolf.config.ini",
    "MMC_DIR": "/project/data/reference/mmc"
}
```

## Create a config file

```bash
pixi run config setup-config --config_store_path /path/to/project/config.json
```

Prompts for each key, defaulting to the environment variable of the same name; any key can be
given as an option instead (`--zarr_storage_path`, `--root_path`, ...). `.json` is appended if
missing, and an existing file is kept unless `--overwrite`. Keys it does not prompt for
(`REGION_LABELS`, `MMC_DIR`, ...) can be added to the file by hand.

## Display a config

```bash
pixi run config display-config /path/to/project/config.json
```

To see what a command will actually use (config file, environment and defaults combined):

```bash
pixi run -e preprocessing spida-P --config /path/to/project/config.json test
```
