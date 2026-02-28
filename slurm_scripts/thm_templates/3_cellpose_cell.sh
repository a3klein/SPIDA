#!/bin/bash
# FILENAME: cellpose_cell.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_cell_{EXP_N}_{REG_N}
#SBATCH -p gpu
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3a_cellpose_cell_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/THM/logs/{EXP_N}/3a_cellpose_cell_{EXP_N}_{REG_N}.out

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

echo -e ${{SLURM_TMPDIR}}
echo -e ${{PIXI_HOME}}
echo -e ${{PIXI_CACHE_DIR}}
echo -e ${{RATTLER_CACHE_DIR}}

cd "$WORKDIR"
sleep 10

pixi install --frozen -e preprocessing
pixi install --frozen -e cellpose

echo -e "\n\tLoading Deconvoluted Images - {REG_N} - {EXP_N}\n"
# Loading in deconvoluted images
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-decon-images \
    {EXPERIMENT} \
    {REGION} \
    {ROOT_PATH} \
    --plot

echo -e "\n\tRunning Cellpose SAM on Region {REG_N} of Experiment {EXP_N}\n"
# Running Cellpose
pixi run --frozen -e cellpose \
    python -m spida.S.cli {CONFIG} \
    run-segmentation-region \
    cellpose \
    {EXPERIMENT} \
    {REGION} \
    --input_dir {ROOT_PATH}/{EXPERIMENT}/out/{REGION}/images \
    --output_dir {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --scale=4 \
    --image_ext=.decon.tif \
    --nuc_stain_name=DAPI \
    --cyto_stain_name=PolyT \
    --flow_threshold=0 \
    --cellprob_threshold=-2 \
    --tile_norm_blocksize=2960 \


echo -e "\n\tLoading Segmentation - {REG_N} - {EXP_N}\n"
# Loading in the segmentation into the zarr store
pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    load-segmentation-region  \
    {EXPERIMENT} \
    {REGION} \
    {SEGMENTATION_DIR}/{EXPERIMENT}/cellpose_cell \
    --type vpt \
    --prefix_name cellpose_cell \
    --transcript-qc \
    --plot

echo -e "\n\tFiltering Cells - {REG_N} - {EXP_N}\n"
# FILTERING 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    filter_cells_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --plot \
    --cutoffs_path {CUTOFFS_PATH}

echo -e "\n\tSetting up AnnData - {REG_N} - {EXP_N}\n"
# SETUP ADATA 
pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    setup_adata_region \
    {EXPERIMENT} \
    {REGION} \
    cellpose_cell \
    --suffix _filt \
    --plot