#!/bin/bash
# FILENAME: 3b_cellpose_cpu.sh
# DESCRIPTION: CPU step (after 3a): post-process the raw cellpose output into the
#              segmentation schema (native), load it into the zarr store, filter
#              cells, and set up the AnnData.

#SBATCH -J cellpose_cpu_{EXP_N}_{REG_N}
#SBATCH --partition={CPU_PARTITION}
#SBATCH --constraint={CPU_CONSTRAINT}
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --requeue
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out

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

# The store and segmentation dir are synced down, so SPIDA's own guards see existing outputs.
# A requeued attempt (spot preemption) replaces its own partial output.
OW="{OVERWRITE_ARG}"
if [[ "${{SLURM_RESTART_COUNT:-0}}" -gt 0 ]]; then OW="--overwrite"; fi

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
echo -e "\nSyncing zarr store, segmentation, and images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
mkdir -p {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}
mkdir -p {SEGMENTATION_DIR}/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ --only-show-errors
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ {SEGMENTATION_DIR}/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ --only-show-errors
# cellpose assigns no transcripts: process-segmentation-region builds cell-by-gene from the raw ones.
aws s3 cp s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/detected_transcripts.parquet {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/ --only-show-errors

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
echo -e "\nPost-processing cellpose output into the segmentation schema (native) - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    process-segmentation-region \
    cellpose \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --prefix_name {CELLPOSE_PREFIX} \
    $OW

echo -e "\nLoading cellpose segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    load-segmentation-region \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --type cellpose \
    --prefix_name {CELLPOSE_PREFIX} \
    --transcript-qc \
    $OW

echo -e "\nFiltering cells - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli --config {CONFIG_PATH} \
    filter-cells-region \
    {CELLPOSE_PREFIX} \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --cutoffs_path {CUTOFFS_PATH} \
    $OW

echo -e "\nSetting up AnnData - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli --config {CONFIG_PATH} \
    setup-adata-region \
    {CELLPOSE_PREFIX} \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --suffix _filt \
    $OW

echo -e "\nGenerating segmentation QC figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site --config {CONFIG_PATH} \
    generate-seg-qc-figs \
    {CELLPOSE_PREFIX} \
    --exp {EXPERIMENT} \
    --reg {REGION} \
    --brain-region {BR} \
    --lab salk {NAMING_MAP_ARG}

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
# The segmentation dir goes back too: process-segmentation-region wrote the schema files there
# (and step 4's proseg starts from its detected_transcripts.csv); /scratch is wiped by the next job.
aws s3 sync {SEGMENTATION_DIR}/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/{REGION}/{CELLPOSE_PREFIX}/ --only-show-errors
aws s3 sync {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
# No step here writes an h5ad (write-adata does); sync only if something did.
if [ -d {ROOT_DIR}/data/anndata ]; then
    aws s3 sync {ROOT_DIR}/data/anndata/ s3://{S3_BUCKET}/spida_outputs/data/anndata/ --only-show-errors
fi
aws s3 sync {ROOT_DIR}/images/ s3://{S3_BUCKET}/spida_outputs/images/ --only-show-errors
