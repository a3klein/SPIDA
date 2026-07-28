#!/bin/bash
# FILENAME: 3b_cellpose_cpu.sh
# CPU step (after 3a): post-process the raw cellpose output into the segmentation schema
# (native), load it into the zarr, filter, and set up the AnnData.

#SBATCH -A mcb130189
#SBATCH -J cellpose_cpu_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3b_cellpose_cpu_{EXP_N}_{REG_N}.out

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

cd "$WORKDIR"
sleep 10
pixi install --frozen -e preprocessing

echo -e "\n\tLoading Deconvoluted Images - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-decon-images \
    {EXPERIMENT} \
    {REGION} \
    {ROOT_PATH} \
    --plot

echo -e "\n\tPost-processing cellpose output into the segmentation schema (native) - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    process-segmentation-region \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --prefix_name cellpose_cell \
    --root_path {ROOT_PATH} \
    --segmentation_store {SEGMENTATION_DIR}

echo -e "\n\tLoading Segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    --segmentation_store {SEGMENTATION_DIR} \
    --type cellpose \
    --prefix_name cellpose_cell \
    --transcript-qc \
    --plot

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --plot \
    --cutoffs_path {CUTOFFS_PATH}

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --suffix _filt \
    --plot
