# SPIDA Configuration Guide

SPIDA uses a per-project configuration file (`.env` or `.json`) to populate standard parameters such as data paths, zarr stores, binary locations, output directories, and other config files (e.g., QC/cutoffs).

## What the config controls

Common keys:
- `RUST_BIN_PATH`, `VPT_BIN_PATH`
- `DECONWOLF_CONFIG`
- `ZARR_STORAGE_PATH`
- `PROCESSED_ROOT_PATH`
- `SEGMENTATION_OUT_PATH`
- `ANNDATA_STORE_PATH`
- `ANNOTATION_STORE_PATH` (or `ANNOTATIONS_STORE_PATH`)
- `IMAGE_STORE_PATH`
- `CUTOFFS_PATH`

## Configuration parameters

- `RUST_BIN_PATH`: Path to the Rust binary used by SPIDA.
- `VPT_BIN_PATH`: Path to the VPT binary used by SPIDA.
- `DECONWOLF_CONFIG`: Path to the Deconwolf configuration file.
- `ZARR_STORAGE_PATH`: Root directory for Zarr stores.
- `PROCESSED_ROOT_PATH`: Root directory for processed outputs.
- `SEGMENTATION_OUT_PATH`: Directory for segmentation outputs.
- `ANNDATA_STORE_PATH`: Directory for AnnData (.h5ad) storage.
- `ANNOTATION_STORE_PATH`: Directory for annotation artifacts/labels.
- `IMAGE_STORE_PATH`: Directory for image outputs/exports.
- `CUTOFFS_PATH`: Path to QC/cutoffs configuration file.

## Create a config file

Use the CLI to generate a `.env` or `.json` file:

```bash
spida config setup-config --config_store_path /path/to/project/.env --ext_type env
```

```bash
spida config setup-config --config_store_path /path/to/project/config --ext_type json
```

Notes:
- `.env` or `.json` is appended automatically if missing.
- Existing files are preserved unless `--overwrite` is used.
- Defaults come from environment variables or built-in defaults.

## Display an existing config

```bash
spida config display-config /path/to/project/.env
```

```bash
spida config display-config /path/to/project/config.json
```

## Loading behavior (summary)

Resolution order is:
1. CLI options
2. Config file
3. Environment variables
4. Built-in defaults

The `ANNOTATIONS_STORE_PATH` env var (or key) is accepted as an alias for `ANNOTATION_STORE_PATH`.
