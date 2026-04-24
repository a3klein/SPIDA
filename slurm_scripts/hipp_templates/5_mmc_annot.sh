#!/bin/bash
# FILENAME = mmc_annot.sh 

#SBATCH -A mcb130189
#SBATCH -J mmc_annot_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/5_mmc_annot_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/5_mmc_annot_{EXP_N}_{REG_N}.out

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

pixi install --frozen -e mmc

PREFIX=cellpose_cell
echo -e "\n\tRunning MapMyCells Annotations - Region {REG_N} of Experiment {EXP_N} - ${{PREFIX}}\n"
pixi run --frozen -e mmc \
    spida-mmc {CONFIG} \
    mmc-annotation-region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    hippocampus-franjic-deep \
    --suffix _filt \
    --mmc_store_path /home/x-aklein2/projects/aklein/BICAN/HIPP/data/mmc \
    --n_cpu 8 \
    --bootstrap_factor 0.9 \
    --bootstrap_iterations 500 \
    --plot \
    --filter_annotations \
    --palette_path /home/x-aklein2/projects/aklein/BICAN/HIPP/data/reference/Franjic/color_palette.json \
    --l1_thr=0.75 \
    --l2_thr=0.75 \
    --l3_thr=0.75 \
    --l4_thr=0.5 \
    --l1_count_thr=20 \
    --l2_count_thr=20 \
    --l3_count_thr=20
sleep 2

PREFIX=proseg_cell
echo -e "\n\tRunning MapMyCells Annotations - Region {REG_N} of Experiment {EXP_N} - ${{PREFIX}}\n"
pixi run --frozen -e mmc \
    spida-mmc {CONFIG} \
    mmc-annotation-region \
    {EXPERIMENT} \
    {REGION} \
    ${{PREFIX}} \
    hippocampus-franjic-deep \
    --suffix _filt \
    --mmc_store_path /home/x-aklein2/projects/aklein/BICAN/HIPP/data/mmc \
    --n_cpu 8 \
    --bootstrap_factor 0.9 \
    --bootstrap_iterations 500 \
    --plot \
    --filter_annotations \
    --palette_path /home/x-aklein2/projects/aklein/BICAN/HIPP/data/reference/Franjic/color_palette.json \
    --l1_thr=0.75 \
    --l2_thr=0.75 \
    --l3_thr=0.75 \
    --l4_thr=0.5 \
    --l1_count_thr=20 \
    --l2_count_thr=20 \
    --l3_count_thr=20
    