#!/bin/bash
# FILENAME: start_script.sh
# DESCRIPTION: This script loads the spatialdata object for a given sample

#SBATCH -A mcb130189
#SBATCH -J start_{EXP_N}_{REG_N}
#SBATCH -p wholenode
#SBATCH -t 00:20:00
#SBATCH --nodes=1
#SBATCH -o /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/start_{EXP_N}_{REG_N}.out
#SBATCH -e /home/x-aklein2/projects/aklein/BICAN/HIPP/logs/{EXP_N}/start_{EXP_N}_{REG_N}.out

# Load the necessary modules
module load modtree/cpu
module list 

export PATH="/home/x-aklein2/.pixi/bin:$PATH"

cd /anvil/projects/x-mcb130189/aklein/BICAN/HIPP/hipp
echo "Setting up .zarr file for region {REG_N} of experiment {EXP_N}"


pixi run -e preprocessing \
    python -m spida.S.cli {CONFIG} \
    ingest-region \
    {EXPERIMENT} \
    {REGION} \
    --plot
