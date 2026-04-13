#!/bin/bash
# FILENAME: 3a_cellpose_gpu.sh
# DESCRIPTION: Run Cellpose segmentation on deconvoluted images.

#SBATCH -J cellpose_gpu_{EXP_N}_{REG_N}
#SBATCH --partition=gpu-l40s-ondemand
#SBATCH --constraint=g6e.4xlarge
#SBATCH --ntasks-per-node=16
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH --exclusive

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
echo -e "\nSyncing deconvoluted images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
aws s3 sync s3://{S3_BUCKET}/spatial_data/{EXPERIMENT}/out/{REGION}/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/ --only-show-errors

tree -L 5 {ROOT_DIR}/{EXPERIMENT}

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
echo -e "\nInstalling pixi environments...\n"
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "cellpose"; then
    pixi install -e cellpose
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env
mkdir /scratch/images

# --- Compute ---
echo -e "\nRunning Cellpose segmentation on region {REGION} of experiment {EXPERIMENT}\n"
pixi run --frozen -e cellpose \
    python -m spida.S.cli --config {CONFIG_PATH} \
    run-segmentation-region \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --scale=4 \
    --image_ext=.decon.tif \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-2 \
    --tile_norm_blocksize=2960

# --- Sync to S3 ---
echo -e "\nSyncing cellpose segmentation to S3...\n"
aws s3 sync {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell/ s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/cellpose_cell/ --only-show-errors
