#!/bin/bash
# FILENAME: 5_mmc_annot.sh
# DESCRIPTION: OPTIONAL (written only when script_aws.py gets --mmc_reference): MapMyCells
#              annotation of {CELLPOSE_PREFIX} and {PROSEG_PREFIX} against {MMC_REFERENCE}, then the labels
#              are copied into the zarr tables (sync-annotations).

#SBATCH -J mmc_{EXP_N}_{REG_N}
#SBATCH --partition={CPU_PARTITION}
#SBATCH --constraint={CPU_CONSTRAINT}
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --requeue
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/5_mmc_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/5_mmc_{EXP_N}_{REG_N}.out

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

# The store is synced down, so SPIDA's own guards see existing outputs.
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
echo -e "\nSyncing zarr store and AnnData copies from S3...\n"
mkdir -p {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION} {ROOT_DIR}/data/anndata
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/anndata/{EXPERIMENT}/{REGION}/ {ROOT_DIR}/data/anndata/{EXPERIMENT}/{REGION}/ --only-show-errors

# --- SPIDA Setup ---
# Refreshed every job: a reused node otherwise keeps an older checkout (and lock file).
rsync -a --delete --exclude='.pixi' /home/ubuntu/aklein/SPIDA/ /scratch/SPIDA/
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
# Node-local package cache: the default (~/.cache on NFS) is one lock shared by every node, so
# concurrent installs serialize. Under .pixi, which the scratch wipe and the rsync both keep.
export PIXI_CACHE_DIR=/scratch/SPIDA/.pixi/cache
pixi install --frozen -e mmc
mkdir -p /scratch/images

# --- Compute ---
for PREFIX in {CELLPOSE_PREFIX} {PROSEG_PREFIX}; do
    echo -e "\nMapMyCells {MMC_REFERENCE} - ${{PREFIX}} - {REG_N} - {EXP_N}\n"
    pixi run --frozen -e mmc \
        python -m spida.I --config {CONFIG_PATH} \
        mmc-annotation-region \
        ${{PREFIX}} \
        {MMC_REFERENCE} \
        --exp {EXPERIMENT} \
        --reg {REGION} \
        --suffix _filt \
        --mmc_store_path {MMC_STORE_PATH} \
        --plot \
        {MMC_ARGS}

    echo -e "\nCopying {MMC_REFERENCE} labels into the zarr table - ${{PREFIX}}\n"
    pixi run --frozen -e mmc \
        python -m spida.I --config {CONFIG_PATH} \
        sync-annotations \
        ${{PREFIX}} \
        --exp {EXPERIMENT} \
        --reg {REGION} \
        --reference {MMC_REFERENCE} \
        $OW
done

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
aws s3 sync {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
aws s3 sync {ROOT_DIR}/data/anndata/ s3://{S3_BUCKET}/spida_outputs/data/anndata/ --only-show-errors
aws s3 sync {ROOT_DIR}/data/annotation/ s3://{S3_BUCKET}/spida_outputs/data/annotation/ --only-show-errors
aws s3 sync {ROOT_DIR}/images/ s3://{S3_BUCKET}/spida_outputs/images/ --only-show-errors
