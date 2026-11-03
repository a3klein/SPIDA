#!/bin/bash
# FILENAME: chain.sh
# DESCRIPTION: Submit the processing chain for experiment {EXPERIMENT}, region {REGION}.
#
# Dependency graph (edges only apply when both ends are enabled):
#   1  ingest / tissue mask ─┬─> 3a cellpose GPU ─> 3b ─> 4 proseg ─> 5 MMC (optional)
#   2a decon DAPI  ──────────┤
#   2b decon PolyT ──────────┘
#
# 3a waits on step 1 because that is where tissue_mask.tif is built, which
# tissue-gated segmentation and --gpu_tile_size=auto both need. 3b waits on
# step 1 as well, for the zarr store. 4 waits on 3b: proseg starts from
# cellpose_cell's detected_transcripts.csv, which 3b writes.
#
# PREP_JOBS (optional, colon-separated job IDs) adds dependencies to 3a, for a
# chain that reuses another chain's decon images / tissue_mask.tif (run with
# steps 3-5 only):  PREP_JOBS=123:124:125 bash chain.sh
# The submitted job IDs are written to jobs.env next to this script.

# Fail fast: if any sbatch returns non-zero, abort before submitting downstream steps.
set -euo pipefail

cd /home/ubuntu/aklein/slurm_jobs/{BR}/{EXP_N}/{REG_N}
echo "Starting processing chain for {EXP_N} / {REG_N}"

DO_1={STEP_1}
DO_2={STEP_2}
DO_3={STEP_3}
DO_4={STEP_4}
DO_5={STEP_5}

# Join the job IDs that actually got submitted into an --dependency=afterok arg.
# Skipped steps contribute nothing, so re-running a subset never waits on a job
# that was never queued.
dep_arg() {{
    local ids=""
    for jid in "$@"; do
        [[ -n "$jid" ]] && ids="${{ids}}:${{jid}}"
    done
    [[ -n "$ids" ]] && echo "--dependency=afterok${{ids}}"
    return 0
}}

jid1="" jid2a="" jid2b="" jid3a="" jid3b="" jid4="" jid5=""
prep="${{PREP_JOBS:-}}"

if [[ $DO_1 == "true" ]]; then
    echo "Step 1: Submitting ingest / tissue mask..."
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
    # 3a needs the deconvolved images (step 2) and the tissue mask (step 1).
    jid3a=$(sbatch $(dep_arg "$jid2a" "$jid2b" "$jid1" ${{prep//:/ }}) 3a_cellpose_gpu.sh | awk '{{print $4}}')
    echo "  cellpose GPU job ID: $jid3a"
    # 3b needs the zarr store (step 1) and the raw boundaries (step 3a).
    jid3b=$(sbatch $(dep_arg "$jid1" "$jid3a") 3b_cellpose_cpu.sh | awk '{{print $4}}')
    echo "  cellpose CPU job ID: $jid3b"
else
    echo "Step 3: Skipping."
fi

if [[ $DO_4 == "true" ]]; then
    echo "Step 4: Submitting ProSeg..."
    jid4=$(sbatch $(dep_arg "$jid3b") 4_proseg.sh | awk '{{print $4}}')
    echo "  proseg job ID: $jid4"
else
    echo "Step 4: Skipping."
fi

if [[ $DO_5 == "true" ]]; then
    echo "Step 5: Submitting MapMyCells..."
    jid5=$(sbatch $(dep_arg "$jid4") 5_mmc_annot.sh | awk '{{print $4}}')
    echo "  MMC job ID: $jid5"
else
    echo "Step 5: Skipping."
fi

cat > jobs.env <<JOBS
jid1=$jid1
jid2a=$jid2a
jid2b=$jid2b
jid3a=$jid3a
jid3b=$jid3b
jid4=$jid4
jid5=$jid5
JOBS

echo "Processing chain for {EXP_N} / {REG_N} has been submitted (job IDs in jobs.env)."
