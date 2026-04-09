#!/bin/bash
# FILENAME: 1_start.sh
# DESCRIPTION: Ingest region into spatialdata zarr store, run transcript QC and hexagon clustering.

#SBATCH -J start_{EXP_N}_{REG_N}
#SBATCH --partition=cpu-ondemand
#SBATCH --constraint=cpu-32vcpu
#SBATCH --mem=128G
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out

# Use EC2 instance role — bypass any SSO credentials inherited via shared /home
unset AWS_PROFILE
unset AWS_DEFAULT_PROFILE
unset AWS_ACCESS_KEY_ID
unset AWS_SECRET_ACCESS_KEY
unset AWS_SESSION_TOKEN
export AWS_CONFIG_FILE=/dev/null
export AWS_SHARED_CREDENTIALS_FILE=/dev/null

# --- Sync from S3 ---
echo -e "\nSyncing raw data from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/ \
    --exclude "region_*" --no-progress
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/ --no-progress

tree -L 5 {ROOT_DIR}/{EXPERIMENT}

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "preprocessing"; then
    pixi install -e preprocessing
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env
mkdir /scratch/images

# --- Compute ---
echo -e "\nIngesting region {REGION} of experiment {EXPERIMENT}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    ingest-region \
    {EXPERIMENT} \
    {REGION}

echo -e "\nRunning Transcript QC\n"
pixi run --frozen -e preprocessing \
    spida-P --config {CONFIG_PATH} \
    transcript-qc \
    {EXPERIMENT} \
    {REGION} \
    --hex_size 30 \
    --min_transcripts 100

echo -e "\nClustering Hexagons\n"
pixi run --frozen -e preprocessing \
    spida-P --config {CONFIG_PATH} \
    cluster-hexes \
    {EXPERIMENT} \
    {REGION} \
    --hex_size 30 \
    --leiden_resolution 0.6 \
    --min_transcripts 100

echo -e "\nGenerating Load QC Figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site --config {CONFIG_PATH} \
    generate-load-qc-figs \
    {EXPERIMENT} \
    {REGION} \
    --brain-region {BR} \
    --lab salk

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
aws s3 sync {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
aws s3 sync {ROOT_DIR}/images/ s3://{S3_BUCKET}/spida_outputs/images/ --only-show-errors
