#!/bin/bash
# FILENAME = proseg_nuc_v4.sh 

#SBATCH -A mcb130189
#SBATCH -J proseg_nuc_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/proseg_nuc_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/proseg_nuc_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL


module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/BICAN/HIPP/hipp

####
echo "Running ProSeg V3.8"
####

PREFIX=proseg_nuc

# SEGMENTATION
pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    run-segmentation-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_nuc/ \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/${{PREFIX}} \
    --voxel-layers=7 \
    --ncomponents=10 \
    --enforce-connectivity=True \
    --nuclear-reassignment-prob=0.05 \
    --cell-compactness=0.05 \
    --diffusion-probability=0.01

# LOAD PROSEG SEGMENTATION
pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/${{PREFIX}} \
    --type proseg \
    --prefix_name ${{PREFIX}} \
    --plot

# FILTERING 
pixi run -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --seg_fam proseg \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/HIPP/config/filtering_cutoffs_proseg.json

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --suffix _filt \
    --plot