#!/bin/bash
# FILENAME: start_script.sh
# DESCRIPTION: This script loads the spatialdata object for a given sample

#SBATCH -A mcb130189
#SBATCH -J start_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH -t 00:40:00
#SBATCH --nodes=1
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/1_start_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/1_start_{EXP_N}_{REG_N}.out

# Load the necessary modules
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
    
echo -e "\nSetting up .zarr file for region {REG_N} of experiment {EXP_N}\n"

pixi install --frozen -e preprocessing

pixi run --frozen -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    ingest-region \
    {EXPERIMENT} \
    {REGION} \
    --plot

echo -e "\n Calling White Matter regions using BCAS1 \n"

pixi run --frozen -e preprocessing \
    python -m spida.P.cli {CONFIG} \
    call-region-tz \
    {EXPERIMENT} \
    {REGION} \
    default \
    wm_regions \
    --plot

echo -e "\n Running Transcript QC \n"
pixi run --frozen -e preprocessing \
    spida-P {CONFIG} \
    transcript-qc \
    {EXPERIMENT} \
    {REGION} \
    --hex_size 30 \
    --min_transcripts 100 \
    --plot
    
echo -e "\n Clustering Hexagons \n"
pixi run --frozen -e preprocessing \
    spida-P {CONFIG} \
    cluster-hexes \
    {EXPERIMENT} \
    {REGION} \
    --hex_size 30 \
    --leiden_resolution 0.6 \
    --min_transcripts 100 \
    --plot

echo -e "\nGenerating Load QC Figures\n"
pixi run --frozen -e preprocessing \
    python -m spida.site {CONFIG} \
    generate-load-qc-figs \
    {EXPERIMENT} \
    {REGION} \
    --brain-region HIPP \
    --lab salk