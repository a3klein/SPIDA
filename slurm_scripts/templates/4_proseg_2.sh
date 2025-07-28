#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189
#SBATCH -J proseg_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/proseg_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/logs/{EXP_N}/proseg_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL

module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

echo "Running ProSeg - {REG_N} - {EXP_N}"


# SEGMENTATION
pixi run -e preprocessing \
    python -m spida.S run \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --input_dir /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/cellpose/ \
    --output_dir /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM \
    --kwargs \
    detect-layers=True \
    voxel-layers=3 \
    ncomponents=10 \
    nbglayers=10 \
    enforce-connectivity=True \
    max-transcript-nucleus-distance=10 \
    nuclear-reassignment-prob=0.01 \

# LOAD PROSEG SEGMENTATION
pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM \
    --type proseg \
    --prefix_name proseg_SAM \
    --plot

# Align PROSEG Transcripts # Remove loading from here?
pixi run -e preprocessing \
    python -m spida.S align-proseg \
    {EXPERIMENT} \
    {REGION} \
    --seed_prefix_name cellpose_SAM \
    --prefix_name proseg_SAM \
    --out_prefix_name proseg_aligned \
    --seg_dir /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM

# # Fixing Loading proseg_aligned of there was an error in the above script
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
    --prefix_name proseg_aligned \
    --load_kwargs \
    cell_metadata_fname=cellpose_metadata.csv \
    cell_by_gene_fname=merged_cell_by_gene.csv \
    cellpose_micron_space_fname=cellpose_polygons.parquet \
    detected_transcripts_fname=merged_transcript_metadata.csv
    
    # cell_metadata_fname=merged_cell_metadata.csv \
    # cell_by_gene_fname=cell_by_gene.csv \
    # cellpose_micron_space_fname=merged_converted_boundaries.parquet \
    # detected_transcripts_fname=detected_transcripts.csv