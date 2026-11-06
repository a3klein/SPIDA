# AWS template changes to relay to Anvil

Changes made to `aws_templates/aws_s3_template/` and `script_aws.py` after `io-redesign-phase1`
(`6b71263`). Each entry says whether the Anvil templates (`anvil_templates/templates/`,
`script.py`) need the same change.

## 2026-11-03

### 1. Output guards against S3 state
Every AWS job wipes `/scratch` first, so SPIDA's `--overwrite` checks only see what the job
synced down. Steps that never sync their outputs down could silently merge into an existing S3
copy.

| Step | Output checked on S3 before compute | On `--overwrite`, upload uses |
|---|---|---|
| 1 | `zarr_store/{EXP}/{REG}/zarr.json`, `images/tissue_mask.tif` (if raster) | `--delete` on the region store |
| 2 | `decon-images/mosaic_{STAIN}_*` | plain sync (same file names) |
| 3a | `segmentation/{EXP}/{REG}/cellpose_cell/` | `--delete` on `cellpose_cell/` |
| 4 | `segmentation/{EXP}/{REG}/proseg_cell/` | `--delete` on `proseg_cell/` |
| 3b, 5 | none needed: store and seg dir are synced down, SPIDA guards fire | unchanged |

Command lines use `$OW` instead of `{OVERWRITE_ARG}`.

**Anvil:** probably not needed. Anvil works on persistent storage, so SPIDA's own guards see
existing outputs. Check that no Anvil step stages to node-local scratch and copies back.

### 2. Requeue-aware overwrite
Every step has `#SBATCH --requeue`. On a requeued attempt (`SLURM_RESTART_COUNT > 0`), `$OW`
becomes `--overwrite`, so a preempted job replaces its own partial output instead of failing on
it. That is safe only because the first attempt's guard passed.

**Anvil:** needed only if Anvil jobs are requeued (preemptible partitions).

### 3. Partition fallback (AWS-specific)
- CPU steps (1, 3b, 4, 5): `{CPU_PARTITION}` / `{CPU_CONSTRAINT}`, defaults
  `cpu-spot,cpu-ondemand` and `cpu-32vcpu|cpu-48vcpu|cpu-64vcpu`. New `script_aws.py` options
  `--cpu_partition`, `--cpu_constraint`.
- Decon (2a, 2b) now use `{GPU_PARTITION}` / `{GPU_CONSTRAINT}` like 3a. Previously pinned to
  `gpu-l40s-ondemand` / `g6e.4xlarge`.
- Defaults are on-demand only (`cpu-ondemand`, `gpu-l40s-ondemand`). Spot-first was tried
  and reverted the same day: spot CPU jobs were killed without being requeued.

**Anvil:** no. Partition names are AWS-specific.

### 4. SPIDA code refreshed every job
`rsync -a --delete --exclude='.pixi'` runs unconditionally. Previously it ran only
`if [ ! -d /scratch/SPIDA ]`, so a reused node kept old code and an old `pixi.lock`.

**Anvil:** check whether Anvil copies SPIDA to node-local storage. If it runs from the shared
checkout, it is not affected.

### 5. Decon input include narrowed
Step 2 syncs `mosaic_{STAIN}_z?.tif` (raw planes only), not `mosaic_{STAIN}_*`. The old pattern
also pulled any `*.decon.tif` left in `images/` by earlier runs (THM1 has them).

**Anvil:** check the equivalent glob or copy in the Anvil decon template.

### 6. 3b downloads the raw `detected_transcripts.parquet`
`process-segmentation-region cellpose` (native backend) runs `partition_transcripts` on
`out/{REGION}/detected_transcripts.parquet`, because cellpose assigns no transcripts. 3b only
synced `images/`, so it failed at `require()`. This gap has been in the template since `eded125`
(2026-07-30) and first ran on THM4.

**Anvil:** probably not affected if 3b reads the raw region directory in place. Check it if Anvil
stages raw data per step.

### 7. Segmentation names are placeholders; `STORE_JOBS`
`cellpose_cell` / `proseg_cell` are no longer hardcoded. Templates use `{CELLPOSE_PREFIX}` /
`{PROSEG_PREFIX}`, set by `--cellpose_prefix` / `--proseg_prefix` (same defaults). That lets two
chains of one region coexist, e.g. HIPc `cellpose_cell_3d` and `cellpose_cell_2d`. `chain.sh` takes
`STORE_JOBS` (adds dependencies to 3b), so a second chain's store writes wait for the first's:
3b–5 download the zarr store, add a table and upload it, and two at once lose a table (root
`zarr.json` holds the consolidated metadata). The generator now errors on
`--force_2d_segmentation --seg_mode 3d_true`.

**Anvil:** the names apply there too. The store race applies only if Anvil chains copy the store.

### 8. `data/anndata` upload only if the folder exists (3b, 4)
No per-region command writes an h5ad: `setup-adata-region` writes only zarr tables, and the h5ad
export is the separate `write-adata`. The unconditional `aws s3 sync /scratch/data/anndata/ ...`
failed on the missing folder ("The user-provided path ... does not exist"). Under `set -e` that
failed 3b after its store and segmentation uploads, before the figure upload. Step 5 is not
affected, because it creates the folder itself.

**Anvil:** check whether its 3b/4 copy `anndata/` out. If so, the same guard applies.

### 9. Node-local pixi cache, `pixi install --frozen`
Every template runs `export PIXI_CACHE_DIR=/scratch/SPIDA/.pixi/cache` and then, unconditionally,
`pixi install --frozen -e <env>`. Before, the `if ! pixi env list | grep` guard was used, with no
`--frozen`. The default cache (`~/.cache/rattler`, on NFS) has one global lock. On 2026-11-03,
fourteen fresh nodes queued on it: env installs took 13–55 min instead of ~4, and the lock holder
copied libtorch headers at ~90 KB/s. The `grep` guard also matched `preprocessing` inside
`preprocessing-gpu`, and an install without `--frozen` can rewrite the node's `pixi.lock`.

**Anvil:** check where its pixi cache lives. A shared filesystem cache has the same lock.
