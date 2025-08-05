#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189
#SBATCH -J setup_{EXP_N}_{REG_N}
#SBATCH -p shared
#SBATCH --time=00:25:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=16gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Setup - {REG_N} - {EXP_N}"

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


# FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_aligned \
    --plot