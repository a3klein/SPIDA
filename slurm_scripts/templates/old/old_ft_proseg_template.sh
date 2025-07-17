#!/bin/bash
# FILENAME = proseg.sh 

#SBATCH -A mcb130189
#SBATCH -J proseg
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=128gb
#SBATCH --time=1:15:00
#SBATCH -p wholenode
#SBATCH -o script_outputs/job-proseg-{ENAME}-{RNAME}.o
#SBATCH -e script_outputs/job-proseg-{ENAME}-{RNAME}.e
#SBATCH --export=ALL

module load modtree/cpu
module list

echo "Aggregating polygons"
fishtank aggregate-polygons \
    -i {REGION_POLYGON_DIR} \
    -o {REGION_POLYGON_FILE} \
    --min_size 200 \
    --z_column "global_z" \
    --save_union True \
    --scale_factor {SCALE_FACTOR} \
    --x_off {X_OFFSET} \
    --y_off {Y_OFFSET} \

echo "Assigning transcripts to polygons"
fishtank assign-spots \
    -i {TRANSCRIPT_INPUT_FILE} \
    -p {REGION_POLYGON_FILE} \
    -o {CELLPOSE_TZ_OUTPUT_FILE} \
    --max_dist 0 \
    --z_column global_z \
    --cell_fill 0 \
    --map_z True \

echo "Running ProSeg"
proseg \
    -x global_x \
    -y global_y \
    -z global_z \
    --gene-column barcode_id \
    --cell-id-column cell \
    --cell-id-unassigned 0 \
    --detect-layers \
    --voxel-layers 3 \
    --ncomponents 10 \
    --nbglayers 10 \
    --enforce-connectivity \
    --max-transcript-nucleus-distance 10 \
    --nuclear-reassignment-prob .2 \
    --output-path {PROSEG_OUTPUT_DIR} \
    {CELLPOSE_TZ_OUTPUT_FILE} \

# NEED TO DECIDE ON REMOVING UNMAPPED CELLPOSE CELLS
echo "Assigning additional transcripts to polygons with ProSeg"
python spida_dev/ultra_scripts/cli.py assign-proseg \
    -t {CELLPOSE_TZ_OUTPUT_FILE} \
    -p {PROSEG_OUTPUT_DIR}/transcript-metadata.csv.gz \
    -o {PROSEG_TZ_OUTPUT_FILE} \
    -gi {PROSEG_GEOM_OUTPUT_FILE} \
    -go {GEOM_OUTPUT_FILE} \
    --codebook {CODEBOOK_PATH} \
    --barcode_column barcode_id \
    --cell_missing 0 \
    --min_jaccard 0.4 \

echo "SUCCESS"