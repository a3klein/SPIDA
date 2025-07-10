#!/bin/bash
# FILENAME: start.sh

#SBATCH -A mcb130189
#SBATCH -J start-MGM1
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -p wholenode
#SBATCH --time=1:30:00
#SBATCH -o test_run.out
#SBATCH -e test_run.out


cd /home/x-aklein2/projects/aklein/SPIDA/
echo "Run Template 19"
./run_templates/template_19.sh

echo "Done"
