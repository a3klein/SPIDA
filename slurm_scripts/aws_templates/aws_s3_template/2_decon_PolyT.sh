#!/bin/bash
# FILENAME: 2_decon_PolyT.sh
# DESCRIPTION: Apply deconwolf deconvolution on the PolyT channel.

#SBATCH -J decon_PolyT_{EXP_N}_{REG_N}
#SBATCH --partition=gpu-l40s-ondemand
#SBATCH --constraint=g6e.4xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_PolyT_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/2_decon_PolyT_{EXP_N}_{REG_N}.out

nvidia-smi
nvcc --version

# Use EC2 instance role — bypass any SSO credentials inherited via shared /home
unset AWS_PROFILE
unset AWS_DEFAULT_PROFILE
unset AWS_ACCESS_KEY_ID
unset AWS_SECRET_ACCESS_KEY
unset AWS_SESSION_TOKEN
export AWS_CONFIG_FILE=/dev/null
export AWS_SHARED_CREDENTIALS_FILE=/dev/null

# --- Sync from S3 ---
echo -e "\nSyncing raw images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/ {ROOT_DIR}/{EXPERIMENT}/out/ \
    --exclude "region_*" --only-show-errors
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ --only-show-errors

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "preprocessing-gpu"; then
    pixi install -e preprocessing-gpu
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env
mkdir /scratch/images

# --- Compute ---
echo -e "\nRunning whole image deconvolution - PolyT - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing-gpu \
    spida-S --config {CONFIG_PATH} \
    decon_image \
    -i {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images \
    --data_org_path {ROOT_DIR}/{EXPERIMENT}/out/dataorganization.csv \
    -o {ROOT_DIR}/{EXPERIMENT}/{REGION}/tile_images \
    --channels PolyT \
    --tile_size 2500 \
    --overlap 400 \
    --z_step 1.5 \
    --filter deconwolf \
    --filter_args tilesize 1000 \
    --gpu \
    --plot_thr

# --- Sync to S3 ---
# Syncs the images directory which now contains the output .decon.tif files
echo -e "\nSyncing deconvoluted images to S3...\n"
aws s3 sync {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/ s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/images/ --only-show-errors
