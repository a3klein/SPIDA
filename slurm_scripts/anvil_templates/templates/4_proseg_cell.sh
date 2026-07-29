#!/bin/bash
# FILENAME: 4_proseg_cell.sh
# proseg segmentation (CPU; shells the rust binary) + native post-process -> load ->
# filter -> setup. proseg is CPU-only, so segment + process share one job.

#SBATCH -A mcb130189
#SBATCH -J proseg_cell_{EXP_N}_{REG_N}
#SBATCH -p shared
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
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

cd "$WORKDIR"
sleep 10
pixi install --frozen -e preprocessing

PREFIX=proseg_cell

echo -e "\n\tRunning ProSeg segmentation (backend) on Region {REG_N} of Experiment {EXP_N} - ${{PREFIX}}\n"
# Segment only -> raw proseg output at {SEGMENTATION_DIR}/{EXP}/{REG}/proseg_cell
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    segment-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --prefix_name ${{PREFIX}} \
    --root_path {ROOT_PATH} \
    --segmentation_store {SEGMENTATION_DIR} \
    --voxel-layers=7 \
    --ncomponents=10 \
    --enforce-connectivity=True \
    --nuclear-reassignment-prob=0.05 \
    --cell-compactness=0.05 \
    --diffusion-probability=0.01 \
    --overwrite=True

echo -e "\n\tPost-processing proseg output into the segmentation schema (native) - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    process-segmentation-region \
    proseg \
    {EXPERIMENT} \
    {REGION} \
    --prefix_name ${{PREFIX}} \
    --root_path {ROOT_PATH} \
    --segmentation_store {SEGMENTATION_DIR}

echo -e "\n\tLoading ProSeg Segmentation - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region \
    {EXPERIMENT} \
    {REGION} \
    --segmentation_store {SEGMENTATION_DIR} \
    --type proseg \
    --prefix_name ${{PREFIX}} \
    --plot

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --seg_fam proseg \
    --plot \
    --cutoffs_path {CUTOFFS_PATH}

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    --suffix _filt \
    --plot
