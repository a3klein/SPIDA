# Project-level datasets

A **project** combines many processed regions (samples) from one or several experiments into one dataset for joint analysis. It is described by a JSON file
and built in three commands:

```bash
pixi run -e preprocessing spida-P --config config.json combine-datasets project.json
pixi run -e preprocessing spida-P --config config.json setup-dataset _projects/<project_name>/<prefix>_table<suffix> --batch_key donor
pixi run -e preprocessing spida-P --config config.json sync-project project.json --columns leiden
```

Every sample must already have its zarr store with the segmentation's shapes, transcripts and
table (i.e. steps 1–3 of the [pipeline](./pipeline.md) have run for that region).

## The project JSON

```json
{
  "project_name": "STUDY_v1",
  "prefix": "cellpose_cell",
  "suffix": "_filt",
  "sample_id": "{brain_region}_{donor}",
  "obs_defaults": {"replicate": "lab1", "major_structure": "STRUCT"},
  "samples": [
    {"exp": "202501011200_STUDY-R1-E-01_VMSC00000", "reg": "region_D1", "obs": {"brain_region": "R1", "donor": "D1"}},
    {"exp": "202501011200_STUDY-R1-E-01_VMSC00000", "reg": "region_D2", "obs": {"brain_region": "R1", "donor": "D2"}},
    {"exp": "202502021300_STUDY-R2-E-01_VMSC00000", "reg": "region_D1", "obs": {"brain_region": "R2", "donor": "D1"}},
    {"exp": "202502021300_STUDY-R2-E-01_VMSC00000", "reg": "region_D3", "obs": {"brain_region": "R2", "donor": "D3", "replicate": "lab2"}}
  ]
}
```

| Key | Required | Meaning |
|---|---|---|
| `project_name` | yes | Output folder name (letters, digits, `_`, `.`, `-`). |
| `prefix` | yes | Segmentation to combine, e.g. `cellpose_cell` or `proseg_cell`. |
| `suffix` | no (`""`) | Table suffix; `_filt` combines the QC-filtered tables. |
| `sample_id` | no (`"{exp}_{reg}"`) | Template for each sample's short name, filled from `exp`, `reg` and the sample's obs. Used in element names and as the sample's coordinate system, so keep it short. |
| `obs_defaults` | no | obs columns set on every sample. |
| `required_obs` | no | obs columns every sample must end up with (after defaults); checked before any data is read. |
| `samples` | yes | One entry per region: `exp`, `reg`, and optional `obs`, which override `obs_defaults`. |

The file is validated as a whole before any data is read: unknown keys, duplicate samples, a
`sample_id` that is not unique or uses a column a sample does not set, and missing inputs are all
reported together.

## combine-datasets

Builds two outputs from each sample's zarr store (`prefix` shapes, transcripts and
`<prefix>_table<suffix>`; images are left out):

| Output | Path | Contents |
|---|---|---|
| Project zarr | `ZARR_STORAGE_PATH/_projects/<project_name>/` | Per sample: `<prefix>_<sample_id>_polygons`, `<prefix>_<sample_id>_transcripts`, `<prefix>_<sample_id>_table<suffix>`. Each sample's `global` coordinate system is renamed to its `sample_id`. |
| Combined h5ad | `ANNDATA_STORE_PATH/_projects/<project_name>/<prefix>_table<suffix>.h5ad` | All samples' cells in one table. |

The combined table:

- `obs`: each sample's obs plus the project obs (`obs_defaults` / sample `obs`) and `sample_id`.
  The region key (`cells_region`) points each cell at its sample's polygons in the project zarr.
- `X` and `layers`: concatenated (genes shared by all samples; a warning lists any that differ).
- `obsm`: keys every sample has with the same width (e.g. `spatial`); per-sample embeddings of
  different widths (e.g. `base_pca`) are left out with a warning.
- `uns["samples"][<sample_id>]`: each sample's own `uns`; `uns["spida_project"]`: the project JSON.

Annotation columns (`mmc_*`, `allcools_*`) are taken from each sample's **zarr** table. If a
sample's h5ad copy has annotation columns its zarr table lacks, a warning names them: run
`spida-I sync-annotations` for that region first.

An existing project is only replaced with `--overwrite`; both outputs are written beside the
target and swapped in, so a failed run never leaves a partial project.

## setup-dataset

Joint normalisation, embedding and clustering on the combined table:

```bash
spida-P setup-dataset _projects/STUDY_v1/cellpose_cell_table_filt --batch_key donor
```

- Re-normalises from `layers["raw"]` (log1p), then PCA, optional Harmony on `--batch_key`,
  neighbours, UMAP and Leiden (consensus clustering by default).
- Writes `X_pca`, `X_pca_harmony` (with `--batch_key`), `X_umap` and `leiden` (prefix them with
  `--key_added`).
- Options: `--leiden_res`, `--knn`, `--min_dist`, `--do_tsne`, `--consensus/--no-consensus`,
  `--n_jobs`.

The dataset name is relative to `ANNDATA_STORE_PATH`; the h5ad is updated in place.

## sync-project

Copies obs columns from the combined h5ad back into the per-sample tables of the project zarr,
splitting cells by the region key — for example joint cluster labels, so they can be plotted on
each sample's polygons:

```bash
spida-P sync-project project.json --columns leiden --columns cell_type
```

Only obs is rewritten. A column that already exists in a project table is refused unless
`--overwrite`; every table is checked before any is written.

## Reading a project

```python
import spatialdata as sd
from spida.P.project_level import read_combined

adata = read_combined(".../anndata_store/_projects/STUDY_v1/cellpose_cell_table_filt.h5ad")
sdata = sd.read_zarr(".../zarr_store/_projects/STUDY_v1")
```

`read_combined` restores the table's region list (h5ad stores it as an array), so the combined
table can be paired with the project zarr's shapes.

Figures for a project go to `IMAGE_STORE_PATH/_projects/<project_name>/` by convention.
