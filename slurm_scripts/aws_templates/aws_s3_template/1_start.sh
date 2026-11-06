#!/bin/bash
# FILENAME: 1_start.sh
# DESCRIPTION: Ingest region into spatialdata zarr store and build the tissue mask.
#              The KDE tissue mask replaces the old transcript-QC / hex-clustering
#              pass; step 3b consumes it as the qc_mask pass region.

#SBATCH -J start_{EXP_N}_{REG_N}
#SBATCH --partition={CPU_PARTITION}
#SBATCH --constraint={CPU_CONSTRAINT}
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --requeue
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out

# Fail fast: any non-zero exit (incl. inside pipelines) aborts the script and the
# SLURM job, so --dependency=afterok in chain.sh stops cascading broken state.
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
s3_guard s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/zarr.json
if [[ -z "{TISSUE_RASTER_ARG}" ]]; then
    s3_guard s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/tissue_mask.tif
fi

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
echo -e "\nSyncing raw data from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/ \
    --exclude "region_*" --no-progress
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/ --no-progress

tree -L 5 {ROOT_DIR}/{EXPERIMENT}

# --- SPIDA Setup ---
# Refreshed every job: a reused node otherwise keeps an older checkout (and lock file).
rsync -a --delete --exclude='.pixi' /home/ubuntu/aklein/SPIDA/ /scratch/SPIDA/
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
# Node-local package cache: the default (~/.cache on NFS) is one lock shared by every node, so
# concurrent installs serialize. Under .pixi, which the scratch wipe and the rsync both keep.
export PIXI_CACHE_DIR=/scratch/SPIDA/.pixi/cache
pixi install --frozen -e preprocessing
mkdir -p /scratch/images

# --- Compute ---
echo -e "\nIngesting region {REGION} of experiment {EXPERIMENT}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    ingest-region \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    $OW

# TISSUE_RASTER_ARG is --no_raster when another chain of this region writes tissue_mask.tif
# (the store still gets its qc_mask shapes / tissue_mask labels).
echo -e "\nComputing tissue mask - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli --config {CONFIG_PATH} \
    tissue-mask \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --rho_min {TISSUE_RHO_MIN} \
    {TISSUE_RASTER_ARG}$OW \
    --plot

echo -e "\nGenerating Load QC Figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site --config {CONFIG_PATH} \
    generate-load-qc-figs \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --brain-region {BR} \
    --lab salk {NAMING_MAP_ARG}

echo -e "\nGenerating default Segmentation QC figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site --config {CONFIG_PATH} \
    generate-seg-qc-figs \
    default \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --brain-region {BR} \
    --lab salk {NAMING_MAP_ARG}

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
if [ -f {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/tissue_mask.tif ]; then
    aws s3 cp {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/tissue_mask.tif \
        s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/tissue_mask.tif --only-show-errors
fi
# ingest-region --overwrite rebuilds the whole store, so the S3 copy is replaced to match
# (--delete, scoped to this region's store); otherwise stale elements would survive.
DELETE=""
if [[ -n "$OW" ]]; then DELETE="--delete"; fi
aws s3 sync {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors $DELETE
aws s3 sync {ROOT_DIR}/images/ s3://{S3_BUCKET}/spida_outputs/images/ --only-show-errors
