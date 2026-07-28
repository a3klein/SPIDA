#!/bin/bash
# FILENAME = proseg_cell_v4.sh 

#SBATCH -A mcb130189
#SBATCH -J proseg_cell_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/4a_proseg_cell_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/4a_proseg_cell_{EXP_N}_{REG_N}.out

module purge
module load modtree/cpu
module list

export PATH="/home/x-aklein2/.pixi/bin:$PATH"
export SLURM_TMPDIR="${{SLURM_TMPDIR:-${{TMPDIR:-/tmp}}}}"
export PIXI_HOME="${{SLURM_TMPDIR}}/pixi_home_${{SLURM_JOB_ID:-$$}}"
export PIXI_CACHE_DIR="/anvil/projects/x-mcb130189/aklein/.cache/rattler/cache"
export RATTLER_CACHE_DIR="/anvil/projects/x-mcb130189/aklein/.cache/rattler/cache"
mkdir -p "$PIXI_HOME" "$PIXI_CACHE_DIR" "$RATTLER_CACHE_DIR"

PROJECT_ROOT="/anvil/projects/x-mcb130189/aklein/BICAN/HIPP/hipp"
WORKDIR="${{SLURM_TMPDIR}}/hipp_${{SLURM_JOB_ID:-$$}}"
rsync -a --delete --exclude '.pixi/envs' "$PROJECT_ROOT/" "$WORKDIR/"

echo -e ${{SLURM_TMPDIR}}
echo -e ${{PIXI_HOME}}
echo -e ${{PIXI_CACHE_DIR}}
echo -e ${{RATTLER_CACHE_DIR}}

cd "$WORKDIR"
sleep 10
pixi --version
pixi install --frozen -e preprocessing

PREFIX=proseg_cell

echo -e "\n\tRunning ProSeg Cell on Region {REG_N} of Experiment {EXP_N} - ${{PREFIX}}\n"

# SEGMENTATION
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    run-segmentation-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell/ \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/${{PREFIX}} \
    --voxel-layers=7 \
    --ncomponents=10 \
    --enforce-connectivity=True \
    --nuclear-reassignment-prob=0.05 \
    --cell-compactness=0.05 \
    --diffusion-probability=0.01 \
    --overwrite=True

echo -e "\n\tLoading ProSeg Segmentation - {REG_N} - {EXP_N}\n"
# LOAD PROSEG SEGMENTATION
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/${{PREFIX}} \
    --type proseg \
    --prefix_name ${{PREFIX}} \
    --plot

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
# FILTERING 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --seg_fam proseg \
    --plot \
    --cutoffs_path /home/x-aklein2/projects/aklein/BICAN/HIPP/config/filtering_cutoffs_proseg.json

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"
# SETUP ADATA 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --suffix _filt \
    --plot