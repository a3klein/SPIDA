#!/bin/bash
# FILENAME: 2_decon_DAPI.sh
# DESCRIPTION: Apply deconwolf deconvolution on the DAPI channel.

#SBATCH -J decon_DAPI_{EXP_N}_{REG_N}
#SBATCH --partition={GPU_PARTITION}
#SBATCH --constraint={GPU_CONSTRAINT}
#SBATCH --requeue
#SBATCH --ntasks-per-node=16
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_DAPI_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_DAPI_{EXP_N}_{REG_N}.out
#SBATCH --exclusive

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
# /scratch starts empty, so nothing local can see outputs that exist only on S3.
# A requeued attempt (spot preemption) replaces its own partial output.
OW="{OVERWRITE_ARG}"
if [[ "${{SLURM_RESTART_COUNT:-0}}" -gt 0 ]]; then OW="--overwrite"; fi
s3_guard() {{
    if [[ -z "$OW" ]] && aws s3 ls "$1" > /dev/null 2>&1; then
        echo "ERROR: $1 exists on S3. Regenerate with --overwrite to replace it." >&2
        exit 1
    fi
}}
s3_guard s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/decon-images/mosaic_DAPI_

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
echo -e "\nSyncing raw images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/ \
    --exclude "region_*" --no-progress
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ \
    --exclude "*" --include "mosaic_DAPI_z?.tif" --no-progress

tree -L 5 {ROOT_DIR}/{EXPERIMENT}

# --- SPIDA Setup ---
# Refreshed every job: a reused node otherwise keeps an older checkout (and lock file).
rsync -a --delete --exclude='.pixi' /home/ubuntu/aklein/SPIDA/ /scratch/SPIDA/
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
# Node-local package cache: the default (~/.cache on NFS) is one lock shared by every node, so
# concurrent installs serialize. Under .pixi, which the scratch wipe and the rsync both keep.
export PIXI_CACHE_DIR=/scratch/SPIDA/.pixi/cache
pixi install --frozen -e preprocessing-gpu
mkdir -p /scratch/images

# --- Compute ---
# Deconvolved planes go to {DECON_DIR} (next to the raw images/, which stay untouched).
echo -e "\nRunning whole image deconvolution - DAPI - {REG_N} - {EXP_N}\n"
mkdir -p {DECON_DIR}
pixi run --frozen -e preprocessing-gpu \
    python -m spida.S.cli --config {CONFIG_PATH} \
    transform-image \
    -i {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images \
    --data_org_path {ROOT_DIR}/{EXPERIMENT}/out/dataorganization.csv \
    --tile_dir {ROOT_DIR}/{EXPERIMENT}/{REGION}/tile_images_DAPI \
    --stitched_dir {DECON_DIR} \
    --channels DAPI \
    --tile_size 2960 \
    --overlap 400 \
    --z_step 1.5 \
    --transform deconwolf \
    {DECON_FILTER_ARGS}\
    --gpu \
    --plot_thr

# --- Sync to S3 ---
echo -e "\nSyncing deconvolved images to S3...\n"
aws s3 sync {DECON_DIR}/ s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/decon-images/ --only-show-errors
