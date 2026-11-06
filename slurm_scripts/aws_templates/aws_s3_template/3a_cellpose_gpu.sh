#!/bin/bash
# FILENAME: 3a_cellpose_gpu.sh
# DESCRIPTION: GPU step: run ONLY the cellpose segmentation backend (raw output).
#              Post-processing into the segmentation schema is CPU/IO-bound and
#              lives in 3b, so the GPU isn't held during it.

#SBATCH -J cellpose_gpu_{EXP_N}_{REG_N}
#SBATCH --partition={GPU_PARTITION}
#SBATCH --constraint={GPU_CONSTRAINT}
#SBATCH --ntasks-per-node=16
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH --exclusive
#SBATCH --requeue

nvidia-smi || true
nvcc --version || true

# Fail fast: any non-zero exit (incl. inside pipelines) aborts the script and the
# SLURM job, so --dependency=afterok in chain.sh stops cascading broken state.
# Diagnostic commands above are non-fatal (|| true) so a missing nvcc doesn't kill the job.
set -euo pipefail

# Use EC2 instance role — bypass any SSO credentials inherited via shared /home
unset AWS_PROFILE
unset AWS_DEFAULT_PROFILE
unset AWS_ACCESS_KEY_ID
unset AWS_SECRET_ACCESS_KEY
unset AWS_SESSION_TOKEN
export AWS_CONFIG_FILE=/dev/null
export AWS_SHARED_CREDENTIALS_FILE=/dev/null

# --- Output guard ---
# /scratch starts empty, so SPIDA's own --overwrite check cannot see outputs that exist
# only on S3. A requeued attempt (spot preemption) replaces its own partial output.
OW="{OVERWRITE_ARG}"
if [[ "${{SLURM_RESTART_COUNT:-0}}" -gt 0 ]]; then OW="--overwrite"; fi
s3_guard() {{
    if [[ -z "$OW" ]] && aws s3 ls "$1" > /dev/null 2>&1; then
        echo "ERROR: $1 exists on S3. Regenerate with --overwrite to replace it." >&2
        exit 1
    fi
}}
s3_guard s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/

# Clean state left by previous jobs on this reused instance. Remove ALL
# experiment data + SPIDA outputs (not just our own — a prior job may have
# been a different experiment that filled /scratch with its own data).
# Preserve /scratch/SPIDA (the pixi env, ~13 GB) so we don't pay the rsync
# + pixi-install cost on every job.
find /scratch -mindepth 1 -maxdepth 1 \
    -not -name 'SPIDA' \
    -not -name 'lost+found' \
    -exec rm -rf {{}} + 2>/dev/null || true

# --- Sync from S3 ---
echo -e "\nSyncing deconvolved images, tissue mask and transform from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images {DECON_DIR}
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/decon-images/ {DECON_DIR}/ --only-show-errors
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ \
    --exclude "*" --include "tissue_mask.tif" --include "micron_to_mosaic_pixel_transform.csv" --only-show-errors

tree -L 5 {ROOT_DIR}/{EXPERIMENT}

# --- SPIDA Setup ---
# Refreshed every job: a reused node otherwise keeps an older checkout (and lock file).
rsync -a --delete --exclude='.pixi' /home/ubuntu/aklein/SPIDA/ /scratch/SPIDA/
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
# Node-local package cache: the default (~/.cache on NFS) is one lock shared by every node, so
# concurrent installs serialize. Under .pixi, which the scratch wipe and the rsync both keep.
export PIXI_CACHE_DIR=/scratch/SPIDA/.pixi/cache
pixi install --frozen -e cellpose
mkdir -p /scratch/images

# --- Compute ---
# IMAGE_EXT is set by script_aws.py (--force_2d_segmentation flag).
# .decon.tif = all deconvolved planes; _z3.decon.tif = the middle z-slice (2D segmentation).
IMAGE_EXT="{IMAGE_EXT}"

# Tiling and any other method kwargs script_aws.py chose. Empty means the whole
# region goes to the GPU in one model.eval call. Deliberately unquoted below so
# it word-splits into separate arguments (and to nothing when empty).
SEG_EXTRA_ARGS="{SEG_EXTRA_ARGS}"

echo -e "\nRunning Cellpose segmentation (backend only) on region {REGION} of experiment {EXPERIMENT}\n"
echo "seg_mode={SEG_MODE}  model={CELLPOSE_MODEL}  extra: ${{SEG_EXTRA_ARGS:-<none>}}"
# Raw boundaries land in the segmentation store under {EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}.
# Stain images come from the deconvolved planes; tissue mask and transform from images/.
pixi run --frozen -e cellpose \
    python -m spida.S.cli --config {CONFIG_PATH} \
    segment-region \
    cellpose \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --prefix_name {CELLPOSE_PREFIX} \
    --image_dir {DECON_DIR} \
    {TISSUE_MASK_ARG}$OW \
    --seg_mode={SEG_MODE} \
    --model={CELLPOSE_MODEL} \
    --scale=4 \
    --image_ext=$IMAGE_EXT \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-2 \
    --tile_norm_blocksize=2960 \
    $SEG_EXTRA_ARGS

# --- Sync to S3 ---
echo -e "\nSyncing cellpose segmentation to S3...\n"
# On overwrite, 3b's schema files from the previous run are stale: replace the whole prefix.
DELETE=""
if [[ -n "$OW" ]]; then DELETE="--delete"; fi
aws s3 sync {SEGMENTATION_DIR}/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ --only-show-errors $DELETE
