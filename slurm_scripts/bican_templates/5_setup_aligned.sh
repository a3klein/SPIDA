#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189-gpu
#SBATCH -J {EXP_N}_{REG_N}_setup_aligned
#SBATCH -p gpu
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=64gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_aligned_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_aligned_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

module load modtree/gpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Setup - {EXP_N} - {REG_N}"

# This loading function is used after re-running transcript assignment on the proseg shapes
pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM \
    --plot \
    --type vpt \
    --prefix_name proseg_aligned \
    --load_kwargs \
    cell_metadata_fname=merged_cell_metadata.csv \
    cell_by_gene_fname=cell_by_gene.csv \
    cellpose_micron_space_fname=merged_converted_boundaries.parquet \
    detected_transcripts_fname=detected_transcripts.csv

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
    max_epochs=100
