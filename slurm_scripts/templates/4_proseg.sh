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
    x-column=global_x \
    y-column=global_y \
    z-column=global_z \
    gene-column=barcode_id \
    cell-id-column=cell_id \
    cell-id-unassigned=0 \
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
pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    {EXPERIMENT} \
    {REGION} \
    /home/x-aklein2/projects/aklein/BICAN/data/segmented/{EXPERIMENT}/proseg_SAM \
    --plot=True \
    --type=vpt \
    --prefix_name=proseg_aligned \
    --cell_metadata_fname=merged_cell_metadata.csv \
    --cell_by_gene_fname=cell_by_gene.csv \
    --cellpose_micron_space_fname=merged_converted_boundaries.parquet \
    --detected_transcripts_fname=detected_transcripts.csv

pixi run -e preprocessing \
    python -m spida.S load_segmentation_region \
    202506171319_BICAN-4x1-GP-E-05_VMSC31910 region_UWA7648 \
    /ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/proseg_SAM \
    --plot \
    --type vpt \
    --prefix_name proseg_aligned \
    --cell_metadata_fname merged_cell_metadata.csv \
    --cell_by_gene_fname cell_by_gene.csv \
    --cellpose_micron_space_fname merged_converted_boundaries.parquet \
    --detected_transcripts_fname detected_transcripts.csv


pixi run -e preprocessing \
    python -m spida.S align-proseg \
    202506171319_BICAN-4x1-GP-E-05_VMSC31910 region_UWA7648 \
    --seed_prefix_name cellpose_SAM \
    --prefix_name proseg_SAM \
    --out_prefix_name proseg_aligned \
    --seg_dir /ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/proseg_SAM