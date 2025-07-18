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


    # FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    202506171319_BICAN-4x1-GP-E-05_VMSC31910 \
    region_UWA7648 \
    cellpose_SAM \
    --plot \
    --cutoffs_path /ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    202506171319_BICAN-4x1-GP-E-05_VMSC31910 \
    region_UWA7648 \
    cellpose_SAM \
    --plot