#!/bin/bash
# FILENAME: cellpose_cell.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_cell_{EXP_N}_{REG_N}
#SBATCH -p gpu
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/cellpose_cell_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/cellpose_cell_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/BICAN/HIPP/hipp

echo "Running Cellpose SAM on Region {REG_N} of Experiment {EXP_N}"

# Loading in deconvoluted images
pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-decon-images \
    {EXPERIMENT} \
    {REGION} \
    {ROOT_PATH} \
    --plot

# Running Cellpose
pixi run -e cellpose \
    python -m spida.S.cli {CONFIG} \
    run-segmentation-region \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {ROOT_PATH}/{EXPERIMENT}/out/{REGION}/images \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --scale=4 \
    --image_ext=.decon.tif \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-4 \
    --tile_norm_blocksize=0 \


# Loading in the segmentation into the zarr store
pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region  \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --type vpt \
    --prefix_name cellpose_cell \
    --plot

# FILTERING 
pixi run -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --plot \
    --cutoffs_path {CUTOFFS_PATH}

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --suffix _filt \
    --plot