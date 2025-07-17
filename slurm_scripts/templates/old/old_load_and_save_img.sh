#!/bin/bash
# FILENAME = load_and_save_img.sh

#SBATCH -A mcb130189
#SBATCH -J load_and_save
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH -p shared
#SBATCH -o script_outputs/job-load.o
#SBATCH -e script_outputs/job-load.o
#SBATCH --export=ALL


# module purge
module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib

FOV_OFFSET={FOV_OFFSET}
fov=$((SLURM_ARRAY_TASK_ID + FOV_OFFSET))

python spida_dev/ultra_scripts/cli.py load_and_save \
    --fov $fov \
    --input {EXPERIMENT_RAW_FILEPATH} \
    --output {REGION_POLYGON_DIR} \
    --data_path {{input}}/data-compressed \
    --data_org_path {{input}}/dataorganization.csv \
    --channels DAPI,PolyT