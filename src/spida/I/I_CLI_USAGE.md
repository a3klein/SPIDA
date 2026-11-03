# SPIDA I module — command-line usage

The **I module** annotates cell types on a segmentation's table with MapMyCells (Allen
`cell-type-mapper`) or ALLCools integration against a single-cell reference, and copies the labels
back into the region's zarr store.

I works on **AnnData only**: it runs in the `mmc` pixi environment and does not need spatialdata.
Labels are written to the table's AnnData copy as `mmc_<reference>_*` / `allcools_<reference>_*`
columns; `sync-annotations` then copies them into the zarr table's obs.

```bash
pixi run -e mmc python -m spida.I --help        # or: pixi run -e mmc spida-I --help
pixi run -e mmc python -m spida.I --config run.json <command> ...
```

## Where things are read and written

| What | Location |
|---|---|
| The table being annotated | `ANNDATA_STORE_PATH/<exp>/<reg>/<prefix>_table<suffix>.h5ad` — copied from the zarr table the first time, reused after (`--refresh_table` re-copies it and drops its annotation columns) |
| Or any AnnData file | `--adata_file path.h5ad` (or a zarr-encoded AnnData) — labels are written back to that file |
| MapMyCells intermediates | `ANNOTATION_STORE_PATH/<exp>/<reg>/<table>/mmc/<reference>/` |
| ALLCools integrator, combined AnnData | `ANNOTATION_STORE_PATH/<exp>/<reg>/<table>/allcools/<reference>/` |
| MapMyCells references | `MMC_DIR/<identifier>/` (built by `mmc-setup`) |
| Figures | `IMAGE_STORE_PATH/<exp>/<reg>/<prefix>/` |

Config keys (`--config` JSON, else environment): `ZARR_STORAGE_PATH`, `ANNDATA_STORE_PATH`,
`ANNOTATION_STORE_PATH`, `IMAGE_STORE_PATH`, `MMC_DIR`, `GENE_PANEL_PATH`. Precedence is CLI option >
`--config` > environment.

Region commands take `--exp` and `--reg`. `--zarr_path`, `--anndata_path`, `--annotations_path`,
`--figures_path` (and `--raw_path`, `--seg_path`) override one region's folders; `--zarr_store`,
`--anndata_store`, `--annotation_store`, `--image_store` override the store roots.

## Commands

| Command | Arguments | What it does |
|---|---|---|
| `mmc-setup` | `REF_PATH IDENTIFIER` | Precompute a MapMyCells reference (stats + markers) under `MMC_DIR/IDENTIFIER` |
| `mmc-annotation-region` | `PREFIX_NAME IDENTIFIER` | Annotate one region's table (or `--adata_file`) |
| `mmc-annotation-experiment` | `PREFIX_NAME IDENTIFIER` | `mmc-annotation-region` for every region of `--exp` |
| `allcools-integration-region` | `PREFIX_NAME REF_PATH` | ALLCools integration of one region's table (or `--adata_file`) |
| `allcools-integration-experiment` | `PREFIX_NAME REF_PATH` | The same for every region of `--exp` |
| `plot-allcools-integration` | `PREFIX_NAME REFERENCE_NAME` | Plot one finished ALLCools integration |
| `sync-annotations` | `PREFIX_NAME` | Copy annotation columns from the AnnData file into the zarr table's obs |
| `sync-annotations-experiment` | `PREFIX_NAME` | `sync-annotations` for every region of `--exp` |
| `moscot-integration` | | Not implemented |

`PREFIX_NAME` is the segmentation name (e.g. `cellpose_cell`); `--suffix` picks the table
(default `_filt`, i.e. `cellpose_cell_table_filt`).

### `mmc-setup REF_PATH IDENTIFIER`

| Option | Default | |
|---|---|---|
| `--hierarchy_list` | `Class,Subclass,Group` | Comma-separated hierarchy levels in the reference obs |
| `--gene_names_path` | searched in `GENE_PANEL_PATH` by identifier | `gene_names.txt`, one gene per line |
| `--gene_name_mapping_path` | searched in `GENE_PANEL_PATH` by identifier | `gene_name_mapping.json`, reference -> panel names |
| `--mmc_store_path` / `--gene_panel_path` | `MMC_DIR` / `GENE_PANEL_PATH` | |
| `--ref_norm` | `log2CPM` | Normalization of the reference `X` |
| `--n_cpu`, `--n_valid`, `--n_per_utility` | 8, 10, 10 | |

### `mmc-annotation-region PREFIX_NAME IDENTIFIER`

| Option | Default | |
|---|---|---|
| `--suffix` | `_filt` | |
| `--adata_file` | | Annotate this file instead; labels are written back to it |
| `--refresh_table` | off | Re-copy the table from the zarr store (drops the file's annotations; the dropped columns are logged) |
| `--mmc_store_path` | `MMC_DIR` | |
| `--n_cpu` | 1 | |
| `--bootstrap_factor`, `--bootstrap_iterations`, `--rng_seed` | 0.8, 100, 13 | |
| `--filter_annotations` | off | Filter calls by score and count |
| `--plot`, `--plot_name`, `--palette_path` | off, `mmc_annotation_plots` | |

### `allcools-integration-region PREFIX_NAME REF_PATH`

`REF_PATH` is the reference RNA AnnData. Main options (all are listed by `--help`):

| Option | Default | |
|---|---|---|
| `--reference_name` | `REF_PATH` stem | Name used in column names and paths |
| `--rna_cell_type_column` | `supercluster_name` | Reference label column |
| `--qry_cluster_column` | `leiden` | Query cluster column (e.g. `base_leiden`) |
| `--gene_rename_dict` | | JSON string or file mapping gene names |
| `--top_deg_genes` | -1 (off) | Top DEGs per cluster used as features |
| `--label_transfer_k` | 100 | |
| `--max_cells_per_cluster`, `--min_cells_per_cluster` | 3000, 0 | Reference downsampling per cluster |
| `--run_joint_embeddings`, `--joint_embedding_leiden_res` | off, 1.5 | |
| `--run_clust_label_transfer` | on | Cluster-to-cluster transfer (`confusion_matrix_*`, `*_cluster_threshold` options) |
| `--save_integrator`, `--save_adata_comb` | on | |
| `--plot` | off | |

Extra `--key=value` options are passed through to the integration (e.g. `--deg_type=cef`).

### `sync-annotations PREFIX_NAME`

| Option | Default | |
|---|---|---|
| `--reference` | every annotation column | Only `mmc_<reference>_*` / `allcools_<reference>_*` |
| `--columns` | | Exactly these obs columns (repeatable) |
| `--adata_file` | the region's AnnData file | |
| `--overwrite` | off | Replace columns already in the zarr table |

Only obs is written; X, layers, obsm and uns in the zarr store are untouched. Without
`--overwrite`, columns already present stop the command before anything is written. The
experiment version skips regions whose columns are already synced.

## Typical workflow

```bash
I="pixi run -e mmc python -m spida.I --config run.json"

# once per reference
$I mmc-setup /refs/hippocampus.h5ad hippocampus-D2P2484 --hierarchy_list Class,Subclass,Group --n_cpu 8

# per region (or *-experiment with --exp only)
$I mmc-annotation-region cellpose_cell hippocampus-D2P2484 --exp EXP --reg region_UCI-5224 --n_cpu 8 --plot
$I allcools-integration-region cellpose_cell /refs/AIT_PU.h5ad --exp EXP --reg region_UCI5224 \
    --rna_cell_type_column Subclass --qry_cluster_column base_leiden --plot

# labels into the zarr table
$I sync-annotations cellpose_cell --exp EXP --reg region_UCI5224
```
