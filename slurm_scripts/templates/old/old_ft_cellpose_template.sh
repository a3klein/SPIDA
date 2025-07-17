#!/bin/bash
# FILENAME = cellpose_slurm.sh 

#SBATCH -A mcb130189-gpu
#SBATCH -J cellpose_run
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH -p gpu
#SBATCH -o script_outputs/job-cellpose-{ENAME}-{RNAME}-%a.o
#SBATCH -e script_outputs/job-cellpose-{ENAME}-{RNAME}-%a.o
#SBATCH --export=ALL
#SBATCH --array={FOV_RANGE}


# module purge
module load modtree/gpu

module load gsl
module load libtiff
module load libpng
module load fftw

module list

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/anvil/projects/x-mcb130189/aklein/programs/gsl/lib

FOV_OFFSET={FOV_OFFSET}
fov=$((SLURM_ARRAY_TASK_ID + FOV_OFFSET))

echo "Running Cellpose $fov" 
python spida_dev/ultra_scripts/cli.py cellpose \
	--fov $fov \
	--input {EXPERIMENT_RAW_FILEPATH} \
	--output {REGION_POLYGON_DIR} \
	--data_path {{input}}/data-compressed \
    --data_org_path {{input}}/dataorganization.csv \
	--channels DAPI,PolyT \
	--downsample 4 \
	--min_size 1000 \
	--do_3D True \
 	--gpu {gpu} {filter}

echo "Cellpose $fov finished"
echo "SUCCESS"

