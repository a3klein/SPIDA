#!/bin/bash
# FILENAME: 3a_cellpose_gpu.sh
# GPU step: run ONLY the cellpose segmentation backend (raw output). Post-processing
# (CPU/IO-bound) is split into 3b so the GPU isn't held during it.

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_gpu_{EXP_N}_{REG_N}
#SBATCH -p gpu
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3a_cellpose_gpu_{EXP_N}_{REG_N}.out

module purge
module load modtree/gpu
module load ngc
module load mpc
module load cuda/12.0.1
module load pytorch/21.09-py3
module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib
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
pixi install --frozen -e cellpose

echo -e "\n\tRunning Cellpose segmentation (backend only) on Region {REG_N} of Experiment {EXP_N}\n"
# Segment only -> raw boundaries at {SEGMENTATION_DIR}/{EXP}/{REG}/cellpose_cell
pixi run --frozen -e cellpose \
    python -m spida.S.cli {CONFIG} \
    segment-region \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --prefix_name cellpose_cell \
    --root_path {ROOT_PATH} \
    --segmentation_store {SEGMENTATION_DIR} \
    --scale=6 \
    --image_ext=.decon.tif \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-2 \
    --tile_norm_blocksize=2960
