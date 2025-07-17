#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189
#SBATCH -J setup_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

odule load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Setup - {REG_N} - {EXP_N}"

# FILTERING 
pixi run -e preprocessing \
    python spida.P.cli filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot=True \
    --cutoffs_path=/home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json # Set Cutoffs Path

# SETUP ADATA 
pixi run -e preprocessing \
    python spida.P.cli setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_SAM \
    --plot=True


# FILTERING 
pixi run -e preprocessing \
    python spida.P.cli filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --plot=True \
    --cutoffs_path=/home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json # Set Cutoffs Path

# SETUP ADATA 
pixi run -e preprocessing \
    python spida.P.cli setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --plot=True