#!/bin/bash
# FILENAME: 4_proseg.sh
# DESCRIPTION: Run ProSeg segmentation, load results, filter cells, and set up AnnData.

#SBATCH -J proseg_{EXP_N}_{REG_N}
#SBATCH --partition=cpu-ondemand
#SBATCH --constraint=cpu-32vcpu
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128G
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/4_proseg_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/4_proseg_{EXP_N}_{REG_N}.out

# Use EC2 instance role — bypass any SSO credentials inherited via shared /home
unset AWS_PROFILE
unset AWS_DEFAULT_PROFILE
unset AWS_ACCESS_KEY_ID
unset AWS_SECRET_ACCESS_KEY
unset AWS_SESSION_TOKEN
export AWS_CONFIG_FILE=/dev/null
export AWS_SHARED_CREDENTIALS_FILE=/dev/null

# --- Sync from S3 ---
echo -e "\nSyncing zarr store and cellpose segmentation from S3...\n"
mkdir -p {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}
mkdir -p {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ --no-progress
aws s3 sync s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/cellpose_cell/ {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell/ --no-progress

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
echo -e "\nRunning ProSeg segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    run-segmentation-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/proseg_cell \
    --voxel-layers=7 \
    --ncomponents=10 \
    --enforce-connectivity=True \
    --nuclear-reassignment-prob=0.05 \
    --cell-compactness=0.05 \
    --diffusion-probability=0.01 \
    --overwrite=True

echo -e "\nLoading ProSeg segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli --config {CONFIG_PATH} \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/proseg_cell \
    --type proseg \
    --prefix_name proseg_cell

echo -e "\nFiltering cells - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli --config {CONFIG_PATH} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_cell \
    --seg_fam proseg \
    --cutoffs_path {CUTOFFS_PATH}

echo -e "\nSetting up AnnData - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli --config {CONFIG_PATH} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_cell \
    --suffix _filt

echo -e "\nGenerating segmentation QC figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site --config {CONFIG_PATH} \
    generate-seg-qc-figs \
    {EXPERIMENT} \
    {REGION} \
    proseg_cell \
    --brain-region {BR} \
    --lab salk

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
aws s3 sync {SEGMENTATION_DIR}/{EXPERIMENT}/proseg_cell/ s3://{S3_BUCKET}/spida_outputs/data/segmentation/{EXPERIMENT}/proseg_cell/ --only-show-errors
aws s3 sync {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ s3://{S3_BUCKET}/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ --only-show-errors
aws s3 sync {ROOT_DIR}/data/anndata/ s3://{S3_BUCKET}/spida_outputs/data/anndata/ --only-show-errors
aws s3 sync {ROOT_DIR}/images/ s3://{S3_BUCKET}/spida_outputs/images/ --only-show-errors
