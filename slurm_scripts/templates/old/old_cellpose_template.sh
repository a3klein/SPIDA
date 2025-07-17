#!/bin/bash
# FILENAME: cellpose_{ENAME}_{RNAME}.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_{ENAME}_{RNAME}
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH -o script_outputs/cellpose_{ENAME}_{RNAME}.o
#SBATCH -e script_outputs/cellpose_{ENAME}_{RNAME}.o
#SBATCH -p gpu
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
export PATH="/home/x-aklein2/.pixi/bin:$PATH"

cd /anvil/projects/x-mcb130189/aklein/SPIDA
echo "Running Cellpose SAM on Region"


pixi run -e cellpose \
    python src/spida/segmentation/main.py run_segmentation \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    -i /anvil/scratch/x-aklein2/BICAN/{EXPERIMENT}/out \
    -o /anvil/projects/x-mcb130189/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose \
    --scale=4 \
    --image_ext={IMAGE_EXT} \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=None \
    --flow_threshold=0 \
    --cellprob_threshold=-4 \
    --tile_norm_blocksize=0
