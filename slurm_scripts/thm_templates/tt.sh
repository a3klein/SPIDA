### CELLPOSE CELL GPU

#!/bin/bash
# FILENAME: cellpose_cell.sh
# DESCRIPTION: This script applies cellpose segmentation on the deconvoluted images

#SBATCH -J cellpose_cell_gpu_MTC_UCI5224
#SBATCH --partition=gpu-l40s-ondemand
#SBATCH --constraint=g6e.4xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/BICAN/CTX/MTC/3_cellpose_cell_MTC_UCI5224.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/BICAN/CTX/MTC/3_cellpose_cell_MTC_UCI5224.out

nvidia-smi
nvcc --version

EXPERIMENT=202505231106_BICAN-4x1-MTC-E-05_VMSC31810
REGION=region_UCI-5224
cd /home/ubuntu/aklein/SPIDA

echo -e "\n\tRunning Cellpose SAM on Region ${REGION} of Experiment ${EXPERIMENT}\n"

# Running Cellpose
pixi run --frozen -e cellpose \
    python -m spida.S.cli \
    run-segmentation-region \
    cellpose \
    ${EXPERIMENT} \
    ${REGION} \
    --input_dir /home/ubuntu/aklein/raw/${EXPERIMENT}/out/${REGION}/images \
    --output_dir /home/ubuntu/aklein/data/segmentation/${EXPERIMENT}/cellpose_cell \
    --scale=4 \
    --image_ext=.decon.tif \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-2 \
    --tile_norm_blocksize=2960 \


#### CELLPOSE CELL CPU

#!/bin/bash
# FILENAME: cellpose_cell_cpu.sh
# DESCRIPTION: This script loads the cellpose segmentation results

#SBATCH -J cellpose_cell_cpu_MTC_UCI5224
#SBATCH --partition=cpu
#SBATCH --constraint=c6id.12xlarge
#SBATCH --output=/home/ubuntu/aklein/spida_logs/BICAN/CTX/MTC/3_cellpose_cell_cpu_MTC_UCI5224.out
#SBATCH --error=/home/ubuntu/aklein/spida_logs/BICAN/CTX/MTC/3_cellpose_cell_cpu_MTC_UCI5224.out
### ADD THE AWS SPECIFIC HEADER HERE ###%

EXPERIMENT=202505231106_BICAN-4x1-MTC-E-05_VMSC31810
REGION=region_UCI-5224
BR=CTX
cd /home/ubuntu/aklein/SPIDA

echo -e "\n\tLoading Deconvoluted Images - {REG_N} - {EXP_N}\n"
# Loading in deconvoluted images
pixi run --frozen -e preprocessing \
    python -m spida.S.cli \
    load-decon-images \
    ${EXPERIMENT} \
    ${REGION} \
    /home/ubuntu/aklein/raw

echo -e "\n\tLoading Segmentation - {REG_N} - {EXP_N}\n"
# Loading in the segmentation into the zarr store
pixi run --frozen -e preprocessing \
    python -m spida.S.cli \
    load-segmentation-region  \
    ${EXPERIMENT} \
    ${REGION} \
    /home/ubuntu/aklein/data/segmentation/${EXPERIMENT}/cellpose_cell \
    --type vpt \
    --prefix_name cellpose_cell \
    --transcript-qc

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
# FILTERING 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli \
    filter_cells_region \
    ${EXPERIMENT} \
    ${REGION} \
    cellpose_cell \
    --cutoffs_path /home/ubuntu/aklein/spida-config/filtering_cutoffs.json

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"
# SETUP ADATA 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli \
    setup_adata_region \
    ${EXPERIMENT} \
    ${REGION} \
    cellpose_cell \
    --suffix _filt \

echo -e "\n Generating QC Figures \n"
pixi run -e preprocessing \
    python -m spida.site \
    generate-seg-qc-figs \
    ${EXPERIMENT} \
    ${REGION} \
    cellpose_cell \
    --brain-region ${BR} \
    --lab salk


vpt --verbose --processes 32 \
    sum-signals \
    --input-images /home/ubuntu/aklein/raw/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/out/region_UCI-5224/images/ \
    --input-boundaries /home/ubuntu/aklein/data/segmentation/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/cellpose_cell/region_UCI-5224/cellpose_micron_space.parquet \
    --input-micron-to-mosaic /home/ubuntu/aklein/raw/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/out/region_UCI-5224/images/micron_to_mosaic_pixel_transform.csv \
    --output-csv /home/ubuntu/aklein/data/segmentation/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/cellpose_cell/region_UCI-5224/sum_signals.csv \
    --overwrite