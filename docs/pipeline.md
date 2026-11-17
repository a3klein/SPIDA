# The processing pipeline

One MERSCOPE region goes from raw output to an annotated table in five steps. Each step is a
handful of SPIDA commands; on a SLURM cluster they are run as a chain of jobs generated from
templates (below).

```
[1] ingest + tissue mask (cpu) ──┬──> [3a] cellpose (gpu) ─> [3b] cellpose post-processing (cpu)
[2] deconvolution DAPI, PolyT ───┘                                  │
    (gpu; needs step 1's mask)                                      v
                                               [4] proseg (cpu) ─> [5] MapMyCells (cpu, optional)
```

| Step | Env | Commands | Writes |
|---|---|---|---|
| 1 start | preprocessing | `spida-S ingest-region`, `spida-P tissue-mask`, `spida-site generate-load-qc-figs` / `generate-seg-qc-figs default` | the region's zarr store (images, transcripts, default segmentation), `tissue_mask.tif`, load-QC figures |
| 2 decon | preprocessing-gpu | `spida-S transform-image --transform deconwolf` (once per stain) | `out/<region>/decon-images/` |
| 3a segment | cellpose (gpu) | `spida-S segment-region cellpose` | raw cellpose output in the segmentation store |
| 3b cellpose | preprocessing | `process-segmentation-region cellpose`, `load-segmentation-region`, `spida-P filter-cells-region`, `setup-adata-region`, seg-QC figures | `cellpose_cell` shapes / transcripts / table (+ `_filt`) in the zarr store |
| 4 proseg | preprocessing | `segment-region proseg` (cellpose as prior), then as 3b | `proseg_cell` elements |
| 5 annotate | mmc | `spida-I mmc-annotation-region`, `sync-annotations` | `mmc_<reference>_*` columns on the tables |

Notes:
- Step 2 selects deconvolution tiles from step 1's tissue mask, so `chain.sh` submits it after
  step 1. `--no_tissue_mask` selects tiles by an Otsu threshold instead (and also turns off the
  mask gate in step 3a, which then needs an explicit `--gpu_tile_size` or `none`).
- Step 3a segments the deconvolved planes (`3d_true`, `3d_stitched`, or `2d` on the middle
  plane), gated to the tissue mask and tiled across the GPU (`--gpu_tile_size auto`).
- Step 4 uses the cellpose result as its nuclear prior, so it runs after 3b.
- Commands take `--exp <experiment> --reg <region>`; paths come from the [config](./configuration.md).
  Run the commands in the environment shown (`pixi run -e <env> ...`).

## Generating job scripts

`slurm_scripts/script.py` (Anvil) and `slurm_scripts/script_aws.py` (AWS ParallelCluster) fill
the templates in `slurm_scripts/*_templates/` for one region (or every `region_*` of an
experiment) and write a `chain.sh` that submits them with `--dependency=afterok`:

```bash
pixi run -e preprocessing python slurm_scripts/script.py 202501011200_STUDY-R1-E-01_VMSC00000 \
    --region_name region_D1 --exp_n R1 --reg_n D1 -b R1 \
    --config_path /path/to/project/config/config.json \
    --output_dir /path/to/project/image_processing \
    --seg_mode 3d_true --cellpose_model cpsam_v2 \
    --decon_iter 50 --decon_tilesize 0 --decon_psigma 1.5 \
    --mmc_reference <reference> \
    --cellpose_gpu_time 9:00:00

bash /path/to/project/image_processing/R1/D1/chain.sh
```

Scripts land in `{output_dir}/{EXP_N}/{REG_N}/`, logs in `{output_dir}/logs/{EXP_N}/`. The config
JSON must set `PROCESSED_ROOT_PATH`, `SEGMENTATION_OUT_PATH` and `CUTOFFS_PATH`; `MMC_DIR` is
needed for step 5. Each job copies the SPIDA checkout when it starts, so later source edits do not
reach a running job (they do reach jobs still in the queue).

Main options (`--help` lists all):

| Option | Default | Meaning |
|---|---|---|
| `--steps` | `all` | Steps written into `chain.sh`, e.g. `1,2` or `3,4,5`. |
| `--lab salk/ucsd`, `--exp_n`, `--reg_n` | full names | Short names for directories and job names. |
| `-b/--brain_region` | required | Label for the site figures. |
| `--seg_mode` | `3d_stitched` | `2d`, `3d_stitched` or `3d_true`. |
| `--cellpose_model` | `cpsam` | Cellpose model. |
| `--tissue_mask/--no_tissue_mask` | on | Gate decon tiles and step 3a to `tissue_mask.tif`. |
| `--tissue_rho_min` | 5.0 | Tissue-mask density threshold; lower keeps sparser tissue. |
| `--decon_iter`, `--decon_tilesize`, `--decon_psigma` | 100, 1000, – | deconwolf settings (`--decon_psigma` needs `--decon_tilesize 0`). |
| `--mmc_reference` | – | MapMyCells reference; enables step 5. |
| `--proseg_cutoffs_path` | `CUTOFFS_PATH` | QC cutoffs for the proseg tables. |
| `--overwrite` | off | Re-run: every writing command gets `--overwrite`. |
| `--*_time` | measured defaults | Wall time per step. |

### Re-running part of a chain

`chain.sh` has one `DO_<step>=true/false` switch per step; edit them and run it again. Skipped
steps add no dependencies. To make step 3a wait for jobs that are already queued (e.g. this
region's decon jobs from an earlier submission):

```bash
PREP_JOBS=<decon_DAPI_jobid>:<decon_PolyT_jobid> bash chain.sh
```

`STORE_JOBS` does the same for step 3b, when a second chain writes to the same zarr store.
Submitted job IDs are written to `jobs.env` next to `chain.sh`.

### Sizing wall times

Ask for about 1.1–1.2x a measured runtime on GPU steps and 1.5–2x on CPU steps; over-asking only
costs queue time. Step 3a dominates. With `--gpu_tile_size auto` its runtime is close to linear in
the GPU work: in our runs (3d_true, cpsam_v2, two tiles at a time on one A100)

```
hours ~ 0.23 + 16.9 s/Mpx * ceil(n_tiles / 2) * tile_Mpx / 3600
```

where `n_tiles` and the tile size are the plan step 3a logs (`tiling plan: ...`). Re-fit these two
numbers on your own hardware from a few completed runs.

## After the pipeline

Combine regions into a project with [project-level datasets](./project_level.md).
