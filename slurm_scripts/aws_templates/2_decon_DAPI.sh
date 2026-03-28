#!/bin/bash
# FILENAME: 2_decon_DAPI.sh
# DESCRIPTION: Apply deconwolf deconvolution on the DAPI channel.

#SBATCH -J decon_DAPI_{EXP_N}_{REG_N}
#SBATCH --partition=gpu-l40s-ondemand
#SBATCH --constraint=g6e.4xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_DAPI_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_DAPI_{EXP_N}_{REG_N}.out

nvidia-smi
nvcc --version

# --- Sync from S3 ---
echo -e "\nSyncing raw images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
rsync -av --exclude "region_*" /s3-data/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/
rsync -av /s3-data/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "preprocessing-gpu"; then
    pixi install -e preprocessing-gpu
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env

# --- Compute ---
echo -e "\nRunning whole image deconvolution - DAPI - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing-gpu \
    spida-S --config {CONFIG_PATH} \
    decon_image \
    -i {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images \
    --data_org_path {ROOT_DIR}/{EXPERIMENT}/out/dataorganization.csv \
    -o {ROOT_DIR}/{EXPERIMENT}/{REGION}/tile_images \
    --channels DAPI \
    --tile_size 2960 \
    --overlap 400 \
    --z_step 1.5 \
    --filter deconwolf \
    --filter_args tilesize 1000 \
    --gpu \
    --plot_thr

# --- Sync to S3 ---
# Syncs the images directory which now contains the output .decon.tif files
echo -e "\nSyncing deconvoluted images to S3...\n"
rsync -av {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ /s3-data/spatial_data/{EXPERIMENT}/out/{REGION}/images/
