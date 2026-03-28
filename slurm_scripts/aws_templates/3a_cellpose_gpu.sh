#!/bin/bash
# FILENAME: 3a_cellpose_gpu.sh
# DESCRIPTION: Run Cellpose segmentation on deconvoluted images.

#SBATCH -J cellpose_gpu_{EXP_N}_{REG_N}
#SBATCH --partition=gpu-l40s-ondemand
#SBATCH --constraint=g6e.4xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out

nvidia-smi
nvcc --version

# --- Sync from S3 ---
echo -e "\nSyncing deconvoluted images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
rsync -av /s3-data/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/

# --- SPIDA Setup ---
if [ ! -d /scratch/SPIDA ]; then
    git clone https://github.com/a3klein/SPIDA.git /scratch/SPIDA
fi
cd /scratch/SPIDA
if ! pixi env list 2>/dev/null | grep -q "cellpose"; then
    pixi install -e cellpose
fi
cp /home/ubuntu/aklein/SPIDA/.env /scratch/SPIDA/.env

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
rsync -av {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell/ /s3-data/spida_outputs/data/segmentation/{EXPERIMENT}/cellpose_cell/
