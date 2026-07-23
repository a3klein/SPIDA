# SPIDA S module CLI Usage

The `spida-S` command is the command-line interface for the SPIDA spatial-processing
(`S`) module: it ingests raw MERSCOPE data into SpatialData, (optionally) deconvolves
images, segments the tissue with a choice of algorithms, post-processes any segmentation
into a standard on-disk **segmentation schema**, and loads that schema back into the
SpatialData object.

Invoke it inside the appropriate pixi environment:

```bash
pixi run -e preprocessing spida-S --help
pixi run -e cellpose     spida-S segment-region cellpose EXP REGION
```

---

## Pipeline at a glance

The commands follow the processing pipeline, in order:

```
ingest-region                 raw MERSCOPE  ->  SpatialData (.zarr)
   │
   ▼ (optional) decon-image   deconvolve stain images for better segmentation
   │
   ▼ segment-region           run a segmentation backend  ->  raw boundaries
   │
   ▼ process-segmentation-region   raw boundaries  ->  segmentation schema
   │                               (native by default; VPT optional)
   ▼ load-segmentation-region  segmentation schema  ->  back into SpatialData
```

Two things to know up front:

- **Segmentation is two commands in two environments.** `segment-region` runs the backend
  in its own env (cellpose → `cellpose`, mesmer → `deepcell`, proseg → `preprocessing`);
  everything else runs in `preprocessing`. This split exists because no single env has both
  the deep-learning segmenters and rasterio/spatialdata.
- **Post-processing is native (pure-Python) by default; VPT is optional.** You do not need
  the VPT binary for normal use. `--backend vpt` is kept only for backwards compatibility
  and is slower (details under *Post-process*, below).

There are also two side paths: **`process-custom-segmentation`** (bring your own
segmentation files) and **`align-segmentation`** (reconcile two segmentations).

## Setup

Configure these via your project config (`.env` or `.json`):
- `PROCESSED_ROOT_PATH` — raw MERSCOPE root (`{root}/{exp}/out/{region}`)
- `SEGMENTATION_OUT_PATH` — segmentation output root (`{root}/{exp}/{method}/{region}`)
- `ZARR_STORAGE_PATH` — SpatialData zarr store
- `VPT_BIN_PATH` — **optional**, only for `--backend vpt`

---

## Commands (in pipeline order)

### 1. Ingest raw data → SpatialData (`ingest-region` / `ingest-all`)

Build a SpatialData object (`.zarr`) from raw MERSCOPE output.

```bash
pixi run -e preprocessing spida-S ingest-region <exp_name> <reg_name> [options]
pixi run -e preprocessing spida-S ingest-all    <exp_name>            [options]
```

**Examples:**

```bash
spida-S ingest-region experiment_1 region_001
spida-S ingest-region experiment_1 region_001 --plot
spida-S ingest-all    experiment_1 --prefix-name batch_001 --plot
```

**Options:** `--type` (default `merscope`), `--prefix-name` (default `default`),
`--source` (default `machine`), `--plot`.

### 2. Deconvolve images (`decon-image`) — optional

Deconvolve large stain images in tiles (DeconWolf) before segmentation. Optional but
improves segmentation for low-boundary-signal tissue.

```bash
pixi run -e preprocessing spida-S decon-image --image_path <path> --data_org_path <path> --channels <channels> [options]
```

**Examples:**

```bash
spida-S decon-image --image_path /path/to/image --data_org_path /path/to/data_org.txt --channels DAPI
spida-S decon-image --image_path /path/to/image --data_org_path /path/to/data_org.txt --channels PolyT,DAPI --tile_size 2960 --overlap 100
```

**Options:** `--output_dir` (default `tiles_output`), `--tile_size` (2960), `--overlap`
(100), `--visualize_grid`, `--z_step` (1.5), `--filter`, `--filter_args`, `--gpu`
(false), `--continue_stalled`, `--plot_thr`.

### 3. Segment a region (`segment-region`)

Run a segmentation *backend* only (produces raw boundaries). Runs in the backend's env.

```bash
pixi run -e <backend_env> spida-S segment-region <method> <exp_name> <reg_name> [options]
```

**Examples:**

```bash
# cellpose -> cellpose env;  mesmer -> deepcell env;  proseg -> preprocessing env
pixi run -e cellpose      spida-S segment-region cellpose experiment_1 region_001
pixi run -e deepcell      spida-S segment-region mesmer   experiment_1 region_001
pixi run -e preprocessing spida-S segment-region proseg   experiment_1 region_001
```

**Methods:** `cellpose`, `mesmer`, `proseg`.

**Options:** `--root_path` (default `PROCESSED_ROOT_PATH`), `--segmentation_store`
(default `SEGMENTATION_OUT_PATH`), `--rust_bin_path`; extra `key=value` args are forwarded
to the backend. (No `--version`: the proseg binary version is fixed at install time —
`--version` applies only to `process-segmentation-region`.)

### 4. Post-process into the segmentation schema (`process-segmentation-region`)

Normalize + post-process the raw output into the segmentation schema. Always runs in
`preprocessing`.

```bash
pixi run -e preprocessing spida-S process-segmentation-region <method> <exp_name> <reg_name> [options]
```

**Examples:**

```bash
# Native (default) — pure-Python, no VPT binary
spida-S process-segmentation-region cellpose experiment_1 region_001

# Legacy VPT fallback (see note below)
spida-S process-segmentation-region cellpose experiment_1 region_001 --backend vpt
```

Writes the schema to `{SEGMENTATION_OUT_PATH}/{exp}/{method}/{region}`:
`boundaries_micron.parquet`, `cell_by_gene.csv`, `detected_transcripts.csv`,
`cell_metadata.csv`, `sum_signals.csv`.

**Options:** `--backend` (`native` default | `vpt`), `--version` (proseg 2/3),
`--root_path`, `--segmentation_store`, `--micron_per_z` (1.5), `--n_z_planes` (7),
`--n_jobs` (7; sum-signals is IO-bound, so 6–8 is the sweet spot — more workers are slower
**and** cost more SUs), `--vpt_bin_path` (only for `--backend vpt`).

> #### VPT backend (`--backend vpt`) — optional, backwards-compatibility only
> The native path is the default and expected one; `--backend vpt` exists only for
> pre-redesign workflows/data. It requires the external `vpt` binary (`VPT_BIN_PATH`) and
> is **slower** than native. Outputs are equivalent up to small, library-version-driven
> numeric differences (not bugs): `sum_signals` intensities differ by ~0.5–1% (GDAL
> `rasterize(all_touched=True)` boundary pixels differ across GDAL versions), `anisotropy`
> differs slightly (shapely `minimum_rotated_rectangle` version), and a tiny fraction
> (~0.16%) of tile-seam boundary slivers may be resolved differently. Cell counts,
> transcript assignments, and gene totals match.

### 5. Process a user-provided segmentation (`process-custom-segmentation`)

Bring your own boundaries (+ optional transcripts / stain images) with whatever column
names you use, and produce the segmentation schema (native; no VPT). Steps are
auto-selected by what you provide (transcripts → partition; images → sum-signals).

```bash
pixi run -e preprocessing spida-S process-custom-segmentation <boundaries_path> <output_dir> [options]
```

**Examples:**

```bash
# 2D boundaries (micron space) + transcripts, columns named yourself
spida-S process-custom-segmentation my_cells.parquet out/ \
    --cell_id_col my_cell \
    --transcripts_path my_transcripts.csv \
    --transcript_z_col z --gene_col target --barcode_col code \
    --transcript_x_col x_um --transcript_y_col y_um

# 3D boundaries (pixel space) + a stain-image stack
spida-S process-custom-segmentation cells.geojson out/ \
    --boundaries_space pixel --micron_to_mosaic_path transform.csv \
    --cell_id_col cell --boundary_z_col z_plane --n_z_planes 7 \
    --images_dir images/ --z_spacing 1.5
```

**Key options:**
- `--boundaries_space` — `micron` (default) or `pixel`
- `--micron_to_mosaic_path` — path to `micron_to_mosaic_pixel_transform.csv` (required for
  `pixel` space and/or when `--images_dir` is given)
- `--z_spacing` (1.5); `--n_z_planes` — 1 ⇒ 2D (polygons expanded into cylinders across the
  transcript/image planes for partitioning; sum-signals uses one plane), >1 ⇒ 3D (boundary
  and transcript plane counts must match)
- `--cell_id_col` — your cell-identifier column (grouping key; preserved in the output; a
  fresh `EntityID` is always assigned). Required for 3D and for metadata merge.
- `--boundary_z_col` — integer z-plane column for 3D input
- `--transcripts_path`, `--transcript_z_col`, `--transcript_z_in_microns`, `--gene_col`,
  `--barcode_col`, `--transcript_x_col`, `--transcript_y_col`
- `--images_dir` — dir of `mosaic_{stain}_z{N}.tif` (enables sum-signals)
- `--segmentation_z_index` — for 2D, the single image plane sum-signals uses (default: the
  middle image plane)
- `--metadata_path` / `--metadata_cell_id_col` — user cell-metadata merged onto the derived
  metadata on your cell-id
- `--n_jobs` (7)

> **Z-indexing note:** transcripts are placed on a z-plane via `--transcript_z_col` (integer
> plane, or micron `ZLevel` with `--transcript_z_in_microns` + `--z_spacing`), while images
> are indexed by the integer `ZIndex` in their `mosaic_*_z{N}.tif` filenames — two
> independent conventions, as in VPT.

### 6. Align two segmentations (`align-segmentation`)

Reconcile two segmentations (e.g. `cellpose_nuc` vs `cellpose_cell`) into one boundary set
by spatial overlap.

```bash
pixi run -e preprocessing spida-S align-segmentation <exp_name> <reg_name> [options]
```

> ⚠️ **Caveat:** not yet migrated to the segmentation-schema redesign and effectively
> **2D-only** (its geometry-conversion call does not thread the z-spacing / 3D flag). Treat
> 3D results with caution.

**Options:** `--prefix1` (default `cellpose_nuc`), `--prefix2` (default `cellpose_cell`),
`--geometry_mode` (`larger`|`prefix1`|`prefix2`|`intersection`), `--cell_id` (default
`EntityID`), `--min_intersection_area`, `--coordinate_system` (default `global`),
`--out_dir_name` (default `align`), `--zarr_store`, `--segmentation_store`, `--root_path`,
`--vpt_bin_path`.

### 7. Load a segmentation into SpatialData (`load-segmentation-region` / `-all`)

Load segmentation-schema output back into the SpatialData object.

```bash
pixi run -e preprocessing spida-S load-segmentation-region <exp_name> <reg_name> <seg_dir> [options]
pixi run -e preprocessing spida-S load-segmentation-all    <exp_name> <seg_dir>            [options]
```

**Examples:**

```bash
spida-S load-segmentation-region experiment_1 region_001 /path/to/segmentation --type cellpose --plot
spida-S load-segmentation-region experiment_1 region_001 /path/to/segmentation --type proseg
```

The loader is backend-agnostic: it resolves the boundary file whether it was written under
the canonical name (`boundaries_micron.parquet`, native) or a legacy name
(`cellpose_micron_space.parquet` / `proseg_polygons.parquet`, e.g. from `--backend vpt` or
older data).

**Options:** `--type` (segmentation method `cellpose`/`proseg`; legacy `vpt` maps to
cellpose), `--prefix-name` (default `default`), `--plot`, `--load_kwargs`.

---

## Reference

### Additional keyword arguments

Segmentation commands forward extra `key=value` arguments to the underlying functions:

```bash
spida-S segment-region cellpose exp1 reg1 diameter=25 flow_threshold=0.5
```

### Help & version

```bash
spida-S --help
spida-S <command> --help
spida-S --version
```

### Command aliases

- `segment-region` → `segment`
- `process-segmentation-region` → `process-segmentation`
- `process-custom-segmentation` → `process-custom`
- `align-segmentation` → `align_segmentation`
- `run-segmentation-region` → `run` *(deprecated — see below)*

### Deprecated / removed

- `run-segmentation-region` / `run` — the old bundled single-step command. **No longer
  runs**; it prints migration guidance pointing at `segment-region` +
  `process-segmentation-region`.
- `segment-experiment` and `align-proseg` were **removed** (batch-over-regions is now done
  per-region via SLURM; `align-proseg` was obsolete under proseg v3's native assignment).

## Author

Amit Klein

---

*This documentation reflects the post-VPT-removal S module: native pure-Python
post-processing by default, with an optional VPT backend for backwards compatibility.*
