#!/bin/bash
# FILENAME: whole_image_dw.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J decon_{EXP}_{REGION}
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=8:00:00
#SBATCH -o script_outputs/decon_{EXP}_{REGION}.o
#SBATCH -e script_outputs/decon_{EXP}_{REGION}.o
#SBATCH -p gpu
#SBATCH --export=ALL

# module purge
module load modtree/gpu

module load gsl
module load libtiff
module load libpng
module load fftw

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib

echo "Running whole image deconvolution"


python spida_dev/ultra_scripts/cli.py decon_image \
-i /anvil/scratch/x-aklein2/BICAN/{EXPERIMENT}/out/{REGION}/images \
--data_org_path /anvil/scratch/x-aklein2/BICAN/{EXPERIMENT}/raw/dataorganization.csv \
-o /anvil/scratch/x-aklein2/BICAN/{EXPERIMENT}/analysis/{REGION}/tile_images \
--channels PolyT \
-ts 2960 \
--overlap 300 \
--z_step 1.5 \
--filter deconwolf \
--filter_args tilesize=900 \
--gpu True \

