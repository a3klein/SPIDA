#!/bin/bash
# FILENAME: chain.sh
# DESCRIPTION: Submit the processing chain for experiment {EXPERIMENT}, region {REGION}.
# Steps 1, 2a, and 2b are submitted immediately (no inter-dependencies).
# Step 3a waits on 2a and 2b; step 3b waits on 1 and 3a; step 4 waits on 3b.

cd /home/ubuntu/aklein/slurm_jobs/{BR}/{EXP_N}/{REG_N}
echo "Starting processing chain for {EXP_N} / {REG_N}"

DO_1={STEP_1}
DO_2={STEP_2}
DO_3={STEP_3}
DO_4={STEP_4}

if [[ $DO_1 == "true" ]]; then
    echo "Step 1: Submitting ingest / transcript QC..."
    jid1=$(sbatch 1_start.sh | awk '{{print $4}}')
    echo "  job ID: $jid1"
else
    echo "Step 1: Skipping."
fi

if [[ $DO_2 == "true" ]]; then
    echo "Step 2: Submitting deconvolution (DAPI + PolyT)..."
    jid2a=$(sbatch 2_decon_DAPI.sh | awk '{{print $4}}')
    jid2b=$(sbatch 2_decon_PolyT.sh | awk '{{print $4}}')
    echo "  DAPI job ID: $jid2a  |  PolyT job ID: $jid2b"
else
    echo "Step 2: Skipping."
fi

if [[ $DO_3 == "true" ]]; then
    echo "Step 3: Submitting Cellpose (GPU + CPU)..."

    # 3a: cellpose GPU — depends on decon steps if they ran
    if [[ -v jid2a ]] && [[ -v jid2b ]]; then
        jid3a=$(sbatch --dependency=afterok:$jid2a:$jid2b 3a_cellpose_gpu.sh | awk '{{print $4}}')
    else
        jid3a=$(sbatch 3a_cellpose_gpu.sh | awk '{{print $4}}')
    fi
    echo "  cellpose GPU job ID: $jid3a"

    # 3b: cellpose CPU — depends on zarr store (step 1) AND cellpose GPU (step 3a)
    if [[ -v jid1 ]] && [[ -v jid3a ]]; then
        jid3b=$(sbatch --dependency=afterok:$jid1:$jid3a 3b_cellpose_cpu.sh | awk '{{print $4}}')
    elif [[ -v jid1 ]]; then
        jid3b=$(sbatch --dependency=afterok:$jid1 3b_cellpose_cpu.sh | awk '{{print $4}}')
    elif [[ -v jid3a ]]; then
        jid3b=$(sbatch --dependency=afterok:$jid3a 3b_cellpose_cpu.sh | awk '{{print $4}}')
    else
        jid3b=$(sbatch 3b_cellpose_cpu.sh | awk '{{print $4}}')
    fi
    echo "  cellpose CPU job ID: $jid3b"
else
    echo "Step 3: Skipping."
fi

if [[ $DO_4 == "true" ]]; then
    echo "Step 4: Submitting ProSeg..."
    if [[ -v jid3b ]]; then
        jid4=$(sbatch --dependency=afterok:$jid3b 4_proseg.sh | awk '{{print $4}}')
    else
        jid4=$(sbatch 4_proseg.sh | awk '{{print $4}}')
    fi
    echo "  proseg job ID: $jid4"
else
    echo "Step 4: Skipping."
fi

echo "Processing chain for {EXP_N} / {REG_N} has been submitted."
