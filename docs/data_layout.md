# Data layout

Where SPIDA reads and writes, and what things are called. Every store root comes from the
[config](./configuration.md); each region lives at `{root}/{experiment}/{region}/`.

## Stores

Example: experiment `202501011200_STUDY-R1-E-01_VMSC00000`, region `region_D1`.

| Config key | Path | Contents |
|---|---|---|
| `PROCESSED_ROOT_PATH` | `{root}/{exp}/out/{reg}/` | MERSCOPE output (input). SPIDA adds `images/tissue_mask.tif` and `decon-images/`. `{root}/{exp}/raw/dataorganization.csv` is read for deconvolution. |
| `ZARR_STORAGE_PATH` | `{root}/{exp}/{reg}/` | The region's SpatialData zarr store (below). |
| `SEGMENTATION_OUT_PATH` | `{root}/{exp}/{reg}/{label}/` | Raw segmenter output and the [segmentation schema](#segmentation-schema) files, one folder per segmentation (`cellpose_cell`, `proseg_cell`, ...). |
| `ANNDATA_STORE_PATH` | `{root}/{exp}/{reg}/{table}.h5ad` | h5ad copies of the store's tables, e.g. `cellpose_cell_table_filt.h5ad`. Annotation steps work on these. |
| `ANNOTATION_STORE_PATH` | `{root}/{exp}/{reg}/{table}/{method}/{reference}/` | Annotation intermediates (MapMyCells / ALLCools inputs and raw results). |
| `IMAGE_STORE_PATH` | `{root}/{exp}/{reg}/{label}/` | Figures; `default/` for load QC, one folder per segmentation. |

Projects use the same roots under `_projects/<project_name>/` (see
[project-level datasets](./project_level.md)).

Any one region's folders can be overridden on a command with `--zarr_path`, `--raw_path`,
`--seg_path`, `--figures_path`, `--anndata_path` or `--annotations_path`.

## Inside a region's zarr store

Element names carry the experiment and region (`<exp>_<reg>`, the store's *namespace*, recorded in
the store's attrs at ingest), except tables:

| Element | Name |
|---|---|
| shapes | `<prefix>_<exp>_<reg>_polygons` (`default` = MERSCOPE's own segmentation); `qc_mask` (tissue polygons) |
| points | `<prefix>_<exp>_<reg>_transcripts` |
| tables | `<prefix>_table`, and `<prefix>_table_filt` after QC filtering |
| images | `default_<exp>_<reg>_z3` |
| labels | `tissue_mask` |

`<prefix>` is the segmentation name, e.g. `cellpose_cell` or `proseg_cell`. A table holds every cell
of its segmentation with QC metrics, `pass_qc` and the applied cutoffs in `uns["cutoffs"]`; the
`_filt` table holds the passing cells with normalisation, embeddings and clusters.

### Labels per region

obs columns such as `donor` are taken from `REGION_LABELS`, a JSON mapping region names to labels
(`{"region_D1": "D1"}`). Without it, the label is the region name without `region_`.

## Segmentation schema

Every segmentation backend (cellpose, proseg, mesmer, or a user-provided segmentation) is
normalised to the same files in its segmentation folder, keyed on `EntityID`, before it is loaded
into the store: `boundaries_micron.parquet`, `cell_by_gene.csv`, `detected_transcripts.csv`,
`cell_metadata.csv` (cell centres as `center_x` / `center_y`) and `sum_signals.csv`.

## For developers: the `spida.io` layer

All reads and writes of SPIDA's own data go through `spida.io`:

| Module | Use |
|---|---|
| `io.layout` | Where a dataset lives: `dataset_from(config, exp, reg, ...)` / `TreeLayout`, element names (`element_keys`, `table_key`), region labels. |
| `io.store` | The zarr store: `open_store(path, selection=...)`, `read_table` (anndata only), `write_element` / `write_table` (replace an element by temp-swap), `update_obs` / `sync_obs` (rewrite only obs), `write_store` (a whole new store, temp-swap), and up-front checks `is_writeable` / `require`. |
| `io.anndata` | SPIDA's AnnData files: `read_adata` (.h5ad or zarr-encoded directory), `read_obs`, `write_adata` (written beside the target, then swapped in). anndata only — usable from the mmc env. |
| `io.load` | Building store elements from files (MERSCOPE ingest, segmentation load, decon images). |
| `io.images` | TIFF reading (windowed, downscaled, masked). |
| `io.formats` / `io.schema` | Files SPIDA reads but does not define (MERSCOPE, segmenter output) / SPIDA's own formats (segmentation schema, tissue mask, decon plane names). |

Files SPIDA does not own — reference atlases, third-party tool scratch files — are read and written
by the code that uses them. Writing steps check all their outputs up front and refuse to replace
existing ones without `--overwrite`.
