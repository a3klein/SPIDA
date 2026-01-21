#!/bin/bash
# FILENAME = proseg_v4.sh 

#SBATCH -A mcb130189
#SBATCH -J proseg_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/proseg_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/proseg_{EXP_N}_{REG_N}.out
#SBATCH --export=ALL


module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
cd /anvil/projects/x-mcb130189/aklein/SPIDA

####
echo "Running ProSeg V3.8"
####

PREFIX=proseg_fv38

# SEGMENTATION
pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    run-segmentation-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose/ \
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

# # ALLCOOLS Integration 
# pixi run -e preprocessing python src/spida/I/main.py allcools-integration-region \
#     {EXPERIMENT} \
#     {REGION} \
#     ${{PREFIX}} \
#     /home/x-aklein2/projects/aklein/BICAN/data/reference/AIT/AIT_{EXP_N}.h5ad \
#     --gene_rename_dict /home/x-aklein2/projects/aklein/BICAN/data/reference/AIT/BG_gene_rename.json \
#     --max_cells_per_cluster 3000 \
#     --min_cells_per_cluster 20 \
#     --top_deg_genes 50 \
#     --rna_cell_type_column Subclass \
#     --qry_cluster_column base_leiden \
#     --run_joint_embeddings \
#     --plot
