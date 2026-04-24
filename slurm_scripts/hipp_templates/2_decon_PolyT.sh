#!/bin/bash
# FILENAME: whole_image_dw.sh

#SBATCH -A mcb130189-gpu
#SBATCH -J decon_{EXP_N}_{REG_N}_P
#SBATCH -p gpu
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/2_decon_{EXP_N}_{REG_N}_P.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/2_decon_{EXP_N}_{REG_N}_P.out

module purge
module load modtree/gpu

module load gsl
module load libtiff
module load libpng
module load fftw

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
echo -e "\n\tRunning whole image deconvolution - PolyT - {REG_N} - {EXP_N}\n"

pixi install --frozen -e preprocessing-gpu
pixi run --frozen -e preprocessing-gpu \
    python -m spida.S.cli {CONFIG} \
    decon_image \
    -i {ROOT_PATH}/{EXPERIMENT}/out/{REGION}/images \
    --data_org_path {ROOT_PATH}/{EXPERIMENT}/raw/dataorganization.csv \
    -o {ROOT_PATH}/{EXPERIMENT}/analysis/{REGION}/tile_images \
    --channels PolyT \
    --tile_size 2960 \
    --overlap 400 \
    --z_step 1.5 \
    --filter deconwolf \
    --filter_args tilesize 1500 \
    --gpu \
    --plot_thr
