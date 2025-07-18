#!/bin/bash
# FILENAME: cellpose.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_{EXP_N}_{REG_N}
#SBATCH -p gpu
#SBATCH --time=1:15:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/cellpose_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/cellpose_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Cellpose SAM on Region {REG_N} of Experiment {EXP_N}"

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
    python spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose \
    --type vpt \
    --prefix_name cellpose_SAM \
    --plot

