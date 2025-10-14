#!/bin/bash
# FILENAME: cellpose.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_{EXP_N}_{REG_N}
#SBATCH -p gpu
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/BG/logs/{EXP_N}/cellpose_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/BG/logs/{EXP_N}/cellpose_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Cellpose SAM on Region {REG_N} of Experiment {EXP_N}"

# running cellpose
pixi run -e cellpose \
    python -m spida.S run \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --input_dir /anvil/scratch/x-aklein2/BICAN/{EXPERIMENT}/out/{REGION}/images \
    --output_dir /anvil/projects/x-mcb130189/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose \
    --kwargs \
    scale=4 \
    image_ext=.decon.tif \
    nuc_stain_name=DAPI \
    cyto_stain_name=PolyT \
    flow_threshold=0 \
    cellprob_threshold=-4 \
    tile_norm_blocksize=0


# Loading in the segmentation into the zarr store
pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose \
    --type vpt \
    --prefix_name cellpose_SAM \
    --plot

# Loading in deconvoluted images
pixi run -e preprocessing python \
    -m spida.S load_decon_images \
    {EXPERIMENT} \
    {REGION} \
    --plot


# FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot