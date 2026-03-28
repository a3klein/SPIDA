#!/bin/bash
# FILENAME: 3b_cellpose_cpu.sh
# DESCRIPTION: Load cellpose segmentation into zarr store, filter cells, and set up AnnData.

#SBATCH -J cellpose_cpu_{EXP_N}_{REG_N}
#SBATCH --partition=cpu-ondemand
#SBATCH --constraint=c6id.12xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/{BR}/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out

# --- Sync from S3 ---
echo -e "\nSyncing zarr store, segmentation, and images from S3...\n"
mkdir -p {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images
mkdir -p {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}
mkdir -p {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell
rsync -av /s3-data/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/ {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/
rsync -av /s3-data/spatial_data/{EXPERIMENT}/out/{REGION}/images/ {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/images/
rsync -av /s3-data/spida_outputs/data/segmentation/{EXPERIMENT}/cellpose_cell/ {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell/

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
echo -e "\nLoading deconvoluted images - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli \
    load-decon-images \
    {EXPERIMENT} \
    {REGION} \
    {ROOT_DIR} \
    --config {CONFIG_PATH}

echo -e "\nLoading cellpose segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --type vpt \
    --prefix_name cellpose_cell \
    --transcript-qc \
    --config {CONFIG_PATH}

echo -e "\nFiltering cells - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --cutoffs_path {CUTOFFS_PATH} \
    --config {CONFIG_PATH}

echo -e "\nSetting up AnnData - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --suffix _filt \
    --config {CONFIG_PATH}

echo -e "\nGenerating segmentation QC figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site \
    generate-seg-qc-figs \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --brain-region {BR} \
    --lab salk \
    --config {CONFIG_PATH}

# --- Sync to S3 ---
echo -e "\nSyncing results to S3...\n"
rsync -av {ROOT_DIR}/data/zarr_store/{EXPERIMENT}/{REGION}/ /s3-data/spida_outputs/data/zarr_store/{EXPERIMENT}/{REGION}/
rsync -av {ROOT_DIR}/data/anndata/ /s3-data/spida_outputs/data/anndata/
rsync -av {ROOT_DIR}/images/ /s3-data/spida_outputs/images/
