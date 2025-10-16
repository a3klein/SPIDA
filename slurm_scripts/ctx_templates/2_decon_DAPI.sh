#!/bin/bash
# FILENAME: whole_image_dw.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J decon_{EXP_N}_{REG_N}_D
#SBATCH -p gpu
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=4
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/CORTEX/logs/{EXP_N}/decon_{EXP_N}_{REG_N}_D.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/CORTEX/logs/{EXP_N}/decon_{EXP_N}_{REG_N}_D.out
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module load gsl
module load libtiff
module load libpng
module load fftw

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running whole image deconvolution - DAPI - {REG_N} - {EXP_N}"

pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    decon_image \
    -i {ROOT_PATH}/{EXPERIMENT}/out/{REGION}/images \
    --data_org_path {ROOT_PATH}/{EXPERIMENT}/raw/dataorganization.csv \
    -o {ROOT_PATH}/{EXPERIMENT}/analysis/{REGION}/tile_images \
    --channels DAPI \
    --tile_size 2960 \
    --overlap 400 \
    --z_step 1.5 \
    --filter deconwolf \
    --filter_args tilesize 1500 \
    --gpu \
    --plot_thr
