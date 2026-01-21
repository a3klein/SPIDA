#!/bin/bash
# FILENAME: chain.sh
# DESCRIPTION: This script chains together the processing steps for a given experiment and region.


cd /home/x-aklein2/projects/aklein/BICAN/HIPP/image_processing/{EXP_N}/{REG_N}
echo "Starting processing chain for Experiment {EXP_N} and Region {REG_N}"

jid1=$(sbatch 1_start.sh | awk '{{print $4}}')
echo "Submitted start script with job ID: $jid1"

jid2=$(sbatch 2_decon_DAPI.sh | awk '{{print $4}}')
jid3=$(sbatch 2_decon_PolyT.sh | awk '{{print $4}}')
echo "Submitted DAPI deconvolution script with job ID: $jid2"
echo "Submitted PolyT deconvolution script with job ID: $jid3"

jid4=$(sbatch --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_nuc.sh | awk '{{print $4}}')
echo "Submitted Cellpose Nuclei segmentation script with job ID: $jid4"
jid5=$(sbatch --dependency=afterok:$jid4 4_proseg_nuc.sh | awk '{{print $4}}')
echo "Submitted PROSEG Nuclei segmentation script with job ID: $jid5"

jid6=$(sbatch --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_cell.sh | awk '{{print $4}}')
echo "Submitted Cellpose Cell segmentation script with job ID: $jid6"
jid7=$(sbatch --dependency=afterok:$jid6 4_proseg_cell.sh | awk '{{print $4}}')
echo "Submitted PROSEG Cell segmentation script with job ID: $jid7"

echo "Processing chain for Experiment {EXP_N} and Region {REG_N} has been initiated."