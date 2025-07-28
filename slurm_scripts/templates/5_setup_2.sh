#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189-gpu
#SBATCH -J setup_{EXP_N}_{REG_N}_2
#SBATCH -p gpu
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=16gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}_P2.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}_P2.out
#SBATCH --export=ALL

module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Setup - {EXP_N} - {REG_N}"

# FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json

### NOT REMOVING DOUBLETS ON CELLPOSE BECUASE I SAW THAT THEY WERE NOT REAL DOUBLETS

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot

# RESOLVI
pixi run -e preprocessing-gpu \
    python -m spida.P resolvi_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --suffix _filt \
    --plot

# PROSEG ALIGNED

# FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json

# # REMOVE DOUBLETS? 
# pixi run -e preprocessing-gpu \
#     python -m spida.P remove-doublets-region \
#     {EXPERIMENT} \
#     {REGION} \
#     proseg_aligned \
#     --threshold 0.5 \
#     --plot

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --suffix _filt \
    --plot

# RESOLVI
pixi run -e preprocessing-gpu \
    python -m spida.P resolvi_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --suffix _filt \
    --plot \
    --model_kwargs \
    max_epochs=300
