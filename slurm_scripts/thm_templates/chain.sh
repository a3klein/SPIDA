#!/bin/bash
# FILENAME: chain.sh
# DESCRIPTION: This script chains together the processing steps for a given experiment and region.


cd /home/x-aklein2/projects/aklein/BICAN/THM/image_processing/{EXP_N}/{REG_N}
echo "Starting processing chain for Experiment {EXP_N} and Region {REG_N}"

DO_1="{STEP_1}"
DO_2="{STEP_2}"
DO_3="{STEP_3}"
DO_4="{STEP_4}"
DO_5="{STEP_5}"

if [[ $DO_1 == "true" ]]; then
    echo "Step 1: Starting initial processing..."
    jid1=$(sbatch 1_start.sh | awk '{{print $4}}')
    echo "Submitted start script with job ID: $jid1"
else
    echo "Step 1: Skipping initial processing."
fi

if [[ $DO_2 == "true" ]]; then
    echo "Step 2: Starting deconvolution..."
    jid2=$(sbatch 2_decon_DAPI.sh | awk '{{print $4}}')
    jid3=$(sbatch 2_decon_PolyT.sh | awk '{{print $4}}')
    echo "Submitted DAPI deconvolution script with job ID: $jid2"
    echo "Submitted PolyT deconvolution script with job ID: $jid3"
else
    echo "Step 2: Skipping deconvolution."
fi

if [[ $DO_3 == "true" ]]; then
    echo "Step 3: Starting Cellpose Segmentation post dependencies..."
    if [[ -v jid1 ]] && [[ -v jid2 ]] && [[ -v jid3 ]]; then
        jid4=$(sbatch --exclude=g004 --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_cell.sh | awk '{{print $4}}')
    elif [[ -v jid1 ]]; then
        jid4=$(sbatch --exclude=g004 --dependency=afterok:$jid1 3_cellpose_cell.sh | awk '{{print $4}}')
    elif [[ -v jid2 ]]; then
        jid4=$(sbatch --exclude=g004 --dependency=afterok:$jid2:$jid3 3_cellpose_cell.sh | awk '{{print $4}}')
    else 
        jid4=$(sbatch --exclude=g004 3_cellpose_cell.sh | awk '{{print $4}}')
    fi
    echo "Submitted Cellpose Cell segmentation script with job ID: $jid4"
else
    echo "Step 3: Skipping Cellpose Segmentation."
fi

if [[ $DO_4 == "true" ]]; then
    echo "Step 4: Starting PROSEG Segmentation post dependencies..."
    if [[ -v jid4 ]]; then
        jid5=$(sbatch --dependency=afterok:$jid4 4_proseg_cell.sh | awk '{{print $4}}')
    else
        jid5=$(sbatch 4_proseg_cell.sh | awk '{{print $4}}')
    fi
    echo "Submitted PROSEG Cell segmentation script with job ID: $jid5"
else
    echo "Step 4: Skipping PROSEG Segmentation."
fi

if [[ $DO_5 == "true" ]]; then
    echo "Step 5: Starting MMC Annotations post dependencies..."
    if [[ -v jid5 ]]; then
        jid10=$(sbatch --dependency=afterok:$jid5 5_mmc_annot.sh | awk '{{print $4}}')
    else
        jid10=$(sbatch 5_mmc_annot.sh | awk '{{print $4}}')
    fi
    echo "Submitted MMC Annotations script with job ID: $jid10"
else
    echo "Step 5: Skipping MMC Annotations."
fi

# jid4=$(sbatch --exclude=g004 --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_cell.sh | awk '{{print $4}}')
# echo "Submitted Cellpose Cell segmentation script with job ID: $jid4"
# jid5=$(sbatch --dependency=afterok:$jid4 4_proseg_cell.sh | awk '{{print $4}}')
# echo "Submitted PROSEG Cell segmentation script with job ID: $jid5"

# jid10=$(sbatch --dependency=afterok:$jid5 5_mmc_annot.sh | awk '{{print $4}}')
# echo "Submitted MMC Annotations script with job ID: $jid10"
# jid4=$(sbatch --exclude=g004 --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_nuc.sh | awk '{{print $4}}')
# echo "Submitted Cellpose Nuclei segmentation script with job ID: $jid4"
# jid5=$(sbatch --dependency=afterok:$jid4 4_proseg_nuc.sh | awk '{{print $4}}')
# echo "Submitted PROSEG Nuclei segmentation script with job ID: $jid5"

# jid6=$(sbatch --exclude=g004 --dependency=afterok:$jid1:$jid2:$jid3 3_cellpose_cell.sh | awk '{{print $4}}')
# echo "Submitted Cellpose Cell segmentation script with job ID: $jid6"
# jid7=$(sbatch --dependency=afterok:$jid6 4_proseg_cell.sh | awk '{{print $4}}')
# echo "Submitted PROSEG Cell segmentation script with job ID: $jid7"

# jid8=$(sbatch --dependency=afterok:$jid4:$jid6 3_cellpose_align.sh | awk '{{print $4}}')
# echo "Submitted Cellpose Align segmentations script with job ID: $jid8"
# jid9=$(sbatch --dependency=afterok:$jid8 4_proseg_align.sh | awk '{{print $4}}')
# echo "Submitted PROSEG Align segmentation script with job ID: $jid9"

# jid10=$(sbatch --dependency=afterok:$jid5:$jid7:$jid9 5_mmc_annot.sh | awk '{{print $4}}')
# echo "Submitted MMC Annotations script with job ID: $jid10"

echo "Processing chain for Experiment {EXP_N} and Region {REG_N} has been initiated."