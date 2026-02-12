#!/bin/bash
# FILENAME = cellpose_align.sh 

#SBATCH -A mcb130189
#SBATCH -J cellpose_align_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=0:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/3c_cellpose_align_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/3c_cellpose_align_{EXP_N}_{REG_N}.out

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

####
echo "Running Cellpose Alignment"
####

PREFIX=cellpose_align

sleep 10
pixi --version
echo -e "\n\tRunning Cellpose Alignment - {REG_N} - {EXP_N}\n"

pixi install --frozen -e preprocessing
pixi run --frozen -e preprocessing \
    spida-S {CONFIG} \
    align-segmentation \
    {EXPERIMENT} \
    {REGION} \
    --out_dir_name cellpose_align

echo -e "\n\tLoading Segmentation - {REG_N} - {EXP_N}\n"

pixi run --frozen -e preprocessing \
    spida-S {CONFIG} \
    load-segmentation-region  \
    {EXPERIMENT} \
    {REGION} \
    /anvil/projects/x-mcb130189/aklein/BICAN/HIPP/data/segmented/{EXPERIMENT}/cellpose_align \
    --type vpt \
    --prefix_name cellpose_align \
    --plot

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
# FILTERING 
pixi run --frozen -e preprocessing \
    spida-P {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_align \
    --plot \
    --cutoffs_path {CUTOFFS_PATH}

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"

# SETUP ADATA 
pixi run --frozen -e preprocessing \
    spida-P {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_align \
    --suffix _filt \
    --plot
