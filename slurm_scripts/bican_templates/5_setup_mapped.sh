#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189-gpu
#SBATCH -J {EXP_N}_{REG_N}_setup_mapped
#SBATCH -p gpu
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=64gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_mapped_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/setup_mapped_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

module load modtree/gpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running Setup - {EXP_N} - {REG_N}"

# This loading function is used with the cellpose_SAM shapes but with the direct proseg output cell_by_gene
cp /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose/{REGION}/cell_metadata.csv \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/cellpose_metadata.csv
cp /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose/{REGION}/cellpose_micron_space.parquet \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/cellpose_polygons.parquet

# Fixing index mismatching + naming schemes 
pixi run -e preprocessing \
    python -m spida.S filt-to-ids \
    --meta_path /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/cellpose_metadata.csv \
    --tz_path /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/merged_transcript_metadata.csv \
    --cbg_path /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/merged_cell_by_gene.csv \
    --geom_path /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM/{REGION}/cellpose_polygons.parquet

pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM \
    --plot \
    --type vpt \
    --prefix_name proseg_mapped \
    --load_kwargs \
    cell_metadata_fname=cellpose_metadata.csv \
    cell_by_gene_fname=merged_cell_by_gene.csv \
    cellpose_micron_space_fname=cellpose_polygons.parquet \
    detected_transcripts_fname=merged_transcript_metadata.csv

# FILTERING 
pixi run -e preprocessing \
    python -m spida.P filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_mapped \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json

# SETUP ADATA 
pixi run -e preprocessing \
    python -m spida.P setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_mapped \
    --suffix _filt \
    --plot

# RESOLVI
pixi run -e preprocessing-gpu \
    python -m spida.P resolvi_region \
    {EXPERIMENT} \
    {REGION} \
    proseg_mapped \
    --suffix _filt \
    --plot \
    --model_kwargs \
    max_epochs=100
