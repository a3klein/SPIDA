#!/bin/bash
# FILENAME: 1_start.sh
# DESCRIPTION: Ingest region into spatialdata zarr store, run transcript QC and hexagon clustering.

#SBATCH -J start_{EXP_N}_{REG_N}
#SBATCH --partition=cpu-ondemand
#SBATCH --constraint=c6id.12xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/1_start_{EXP_N}_{REG_N}.out

# --- Sync from S3 ---
echo -e "\nSyncing raw data from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}
rsync -av --exclude "region_*" /s3-data/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/
rsync -av /s3-data/spatial_data/{EXPERIMENT}/out/{REGION}/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "preprocessing"; then
    pixi install -e preprocessing
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env

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
rsync -av {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ /s3-data/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/
rsync -av {ROOT_DIR}/images/ /s3-data/spida_outputs/images/
