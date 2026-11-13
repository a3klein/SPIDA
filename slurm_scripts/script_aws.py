"""
script_aws.py — Generate AWS ParallelCluster Slurm job scripts for the SPIDA pipeline.

Each (EXPERIMENT, REGION, BRAIN_REGION) triple gets its own directory under:
    {output_dir}/{BR}/{EXP_N}/{REG_N}/

containing populated versions of every *.sh template in --template_dir, plus a
chain.sh that chains them with sbatch --dependency.

The brain-region config file ({config_dir}/{BR}_config.json) is required and
provides ROOT_DIR, SEGMENTATION_DIR, and CUTOFFS_PATH.  It is also passed
verbatim as CONFIG_PATH so SPIDA commands receive --config at runtime.

Templates: aws_templates/aws_s3_template/ (data moves via 'aws s3 sync' using the EC2 instance IAM role).
The rsync / FUSE-mount variant was archived in 2026-10 (aws_templates/archive/aws_rsync_template/).

Deconvolved planes go to {ROOT_DIR}/{EXPERIMENT}/out/{REGION}/decon-images/ (raw images/ untouched apart
from tissue_mask.tif). Step 5 (MapMyCells) is written only with --mmc_reference, since only some brain
regions have a reference. --overwrite adds --overwrite to every writing command (re-runs).

Output guards: each job starts on an empty /scratch, so steps that do not first sync their
outputs down (1, 2, 3a, 4's proseg dir) check S3 before any compute and exit unless
--overwrite. With --overwrite, steps 1, 3a and 4 replace their S3 prefix (scoped --delete).
Every step is --requeue'd; a requeued attempt (SLURM_RESTART_COUNT > 0) runs as --overwrite.

Step 3a segmentation is configured here rather than left to the CellposeConfig
defaults, so a generated script records the run it performed. By default it tiles
the region ('auto' grid) and gates inference to the tissue mask that step 1
builds; --gpu_tile_size none --no_tissue_mask restores the older whole-region,
ungated behaviour.

Usage examples
--------------
# Single region (aws s3 sync, default):
python script_aws.py 202505231106_BICAN-4x1-MTC-E-05_VMSC31810 \\
    --brain_region CTX --region_name region_UCI-5224 --lab salk


# All regions in an experiment (discovers region_* dirs from a local path):
python script_aws.py 202505231106_BICAN-4x1-MTC-E-05_VMSC31810 \\
    --brain_region CTX \\
    --data_path /home/ubuntu/aklein/raw/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/out \\
    --lab salk

# Re-run only steps 3 and 4:
python script_aws.py ... --steps 3,4

# Untiled whole-region segmentation (the pre-tiling behaviour):
python script_aws.py ... --gpu_tile_size none --no_tissue_mask

# True-3D segmentation:
python script_aws.py ... --seg_mode 3d_true
"""

import os
import sys
import json
import stat
from pathlib import Path

import click


# ---------------------------------------------------------------------------
# Default paths (relative to this script so they work in any SPIDA clone)
# ---------------------------------------------------------------------------

_SCRIPT_DIR = Path(__file__).parent
_TEMPLATE_DIR = _SCRIPT_DIR / "aws_templates" / "aws_s3_template"
_DEFAULT_OUTPUT_DIR = Path("/home/ubuntu/aklein/slurm_jobs")
_DEFAULT_CONFIG_DIR = Path("/home/ubuntu/aklein/spida_config")
_DEFAULT_LOG_DIR = Path("/home/ubuntu/aklein/spida_logs")
_DEFAULT_S3_BUCKET = "salk-workstation-data-dev-020125249408"


# ---------------------------------------------------------------------------
# Name-derivation helpers (same logic as script.py / spida._constants)
# ---------------------------------------------------------------------------

def rename_exp_salk(x: str) -> str:
    """Derive a short experiment name from a Salk-style full experiment string."""
    exp_n = x.split("_")[1]
    temp = exp_n.split("-")[-2]
    if len(temp) == 1:
        temp = exp_n.split("-")[-3]
    return temp


def rename_reg_salk(x: str) -> str:
    """Derive a short region name from a Salk-style full region string."""
    reg_n = x.split("_")[1]
    temp = reg_n.split("-")[0]
    if len(temp) == 3:
        temp = "".join(reg_n.split("-")[0:2])
    return temp


def rename_exp_ucsd(x: str) -> str:
    """Derive a short experiment name from a UCSD-style full experiment string."""
    return x.split("_")[1].split("BICAN")[-1]


def rename_reg_ucsd(x: str) -> str:
    """Derive a short region name from a UCSD-style full region string."""
    return x.split("_")[1].split("Q0")[0]


def _parse_list(arg: str) -> list:
    """Parse a comma-separated list CLI argument."""
    if arg is None:
        return None
    return arg.split(",")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command("config-aws")
@click.argument("experiment_name")
@click.option(
    "--brain_region", "-b",
    required=True,
    help="Brain region abbreviation (e.g. CTX, THM, HIPP). "
         "Used for directory structure and to locate {config_dir}/{BR}_config.json.",
)
@click.option(
    "--region_name", "-r",
    default=None,
    help="Full region name (e.g. region_UCI-5224). "
         "If omitted, all region_* directories found under --data_path are processed.",
)
@click.option(
    "--lab",
    default=None,
    type=click.Choice([None, "salk", "ucsd"], case_sensitive=False),
    help="Lab naming convention for deriving EXP_N / REG_N. "
         "If omitted, the full experiment/region name is used as-is.",
)
@click.option(
    "data_path",
    "--data_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default=None,
    help="Path to the experiment's 'out/' directory on the local filesystem. "
         "Required when --region_name is omitted (used to discover region_* dirs).",
)
@click.option(
    "--output_dir",
    type=click.Path(dir_okay=True, path_type=Path),
    default=_DEFAULT_OUTPUT_DIR,
    show_default=True,
    help="Root directory for generated job scripts.",
)
@click.option(
    "--template_dir",
    type=click.Path(exists=False, dir_okay=True, path_type=Path),
    default=None,
    show_default=False,
    help="Override the template directory (default: aws_templates/aws_s3_template/ next to this script).",
)
@click.option(
    "--s3_bucket",
    default=_DEFAULT_S3_BUCKET,
    show_default=True,
    help="S3 bucket name (without s3:// prefix). Fills {S3_BUCKET} in s3 templates.",
)
@click.option(
    "--config_dir",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default=_DEFAULT_CONFIG_DIR,
    show_default=True,
    help="Directory containing per-brain-region config JSON files "
         "(e.g. CTX_config.json, THM_config.json).",
)
@click.option(
    "--exp_n",
    default=None,
    help="Override the derived short experiment name (EXP_N).",
)
@click.option(
    "--reg_n",
    default=None,
    help="Override the derived short region name (REG_N). "
         "Only applies when a single --region_name is provided.",
)
@click.option(
    "--steps",
    type=_parse_list,
    default="all",
    show_default=True,
    help="Comma-separated list of steps to enable in chain.sh (1-5). "
         "Default 'all' enables every step (5 only with --mmc_reference). Step 1 builds the "
         "tissue mask that step 3a's tiling and gating need.",
)
@click.option("--proseg_cutoffs_path", default=None, help="QC cutoffs JSON for proseg_cell filtering (default: CUTOFFS_PATH in the config).")
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Pass --overwrite to every writing command (re-running a region over existing outputs). "
         "Off by default: a re-run without it stops before any compute.",
)
@click.option(
    "--tissue_raster/--no_tissue_raster",
    "tissue_raster",
    default=True,
    show_default=True,
    help="Step 1 writes images/tissue_mask.tif. --no_tissue_raster writes only the store's mask, "
         "for a second chain of the same region that reuses another chain's tissue_mask.tif.",
)
@click.option(
    "--mmc_reference",
    default=None,
    help="MapMyCells reference identifier (under --mmc_store_path). Enables step 5; "
         "omit for brain regions without a reference.",
)
@click.option(
    "--mmc_store_path",
    default=None,
    help="MapMyCells reference store (default: MMC_DIR in the brain-region config).",
)
@click.option(
    "--mmc_args",
    default="--n_cpu 8 --bootstrap_iterations 100",
    show_default=True,
    help="Further mmc-annotation-region options, verbatim (e.g. --filter_annotations --l1_thr=0.75).",
)
@click.option(
    "--force_2d_segmentation",
    is_flag=True,
    default=False,
    help="Use the middle z-slice (_z3.decon.tif) for cellpose input instead of the full "
         "deconvolved image (.decon.tif). Use for 3D experiments (PSC-01, MTC_REPEAT) "
         "where segmentation should stay consistent with 2D brain regions.",
)
# --- deconvolution (step 2) -------------------------------------------------
# These were hardcoded in the templates as `--transform_args tilesize 1000`, with the
# iteration count not reachable at all. Passing them explicitly follows the same rule
# the segmentation options already do: a run should be reproducible from the generated
# script rather than from whatever a default happened to be at submission time.
@click.option(
    "--decon_iter",
    type=int,
    default=100,
    show_default=True,
    help="deconwolf Richardson-Lucy iterations. NOTE: 100 is the long-standing "
         "production value and is deep in the RL 'night sky' regime -- measured "
         "fragmentation at iter=100 is ~24x worse than not deconvolving at all. "
         "Lower values are under evaluation; see dark_staining/docs/PRODUCTION_PLAN.md.",
)
@click.option(
    "--decon_tilesize",
    type=int,
    default=1000,
    show_default=True,
    help="deconwolf's INTERNAL tile size, applied inside each SPIDA --tile_size tile. "
         "0 disables it (untiled). Untiled is ~3x faster at a 2960px tile (measured: "
         "47.5s fixed overhead tiled vs 14.3s untiled) but raises peak VRAM from "
         "~2.5 GB to ~12 GB, so it needs a GPU with headroom.",
)
@click.option(
    "--decon_psigma",
    type=float,
    default=None,
    help="deconwolf --psigma, its Anscombe-domain regulariser for the RL sparsity "
         "artifact. REQUIRES --decon_tilesize 0: deconwolf silently ignores psigma "
         "whenever a tilesize is set (verified -- tiled psigma 0.5 and 1.5 give "
         "bit-identical output).",
)
@click.option(
    "--cellpose_prefix",
    default="cellpose_cell",
    show_default=True,
    help="Segmentation name for the cellpose run (3a/3b; proseg's prior). Distinct names let "
         "two chains of one region coexist, e.g. cellpose_cell_2d and cellpose_cell_3d.",
)
@click.option(
    "--proseg_prefix",
    default="proseg_cell",
    show_default=True,
    help="Segmentation name for the proseg run (step 4).",
)
@click.option(
    "--seg_mode",
    type=click.Choice(["2d", "3d_stitched", "3d_true"], case_sensitive=False),
    default="3d_stitched",
    show_default=True,
    help="Cellpose segmentation mode for step 3a. '2d' segments one z-slice; "
         "'3d_stitched' segments each slice and links them by IoU; "
         "'3d_true' uses joint 3D flows and anisotropy.",
)
@click.option(
    "--cellpose_model",
    type=click.Choice(["cpsam", "cpsam_v2", "cpdino", "cpdino-vitb"], case_sensitive=False),
    default="cpsam",
    show_default=True,
    help="Cellpose model for step 3a. Passed explicitly so the run is reproducible "
         "from the generated script rather than from the config default.",
)
@click.option(
    "--gpu_tile_size",
    default="auto",
    show_default=True,
    help="Outer GPU tile edge for step 3a, in full-resolution pixels: 'auto' searches "
         "for the grid sending the least area to the GPU (requires the step-1 tissue "
         "mask), an int gives a square tile, '<y>x<x>' a non-square one, and 'none' "
         "disables tiling so the whole region goes to the GPU in one model.eval call.",
)
@click.option(
    "--gpu_tile_overlap",
    type=int,
    default=500,
    show_default=True,
    help="Outer-tile overlap in full-resolution pixels. Only used when tiling is on. "
         "Too small an overlap loses cells that straddle a tile edge.",
)
@click.option(
    "--tissue_mask/--no_tissue_mask",
    "use_tissue_mask",
    default=True,
    show_default=True,
    help="Gate decon tile selection (step 2) and step 3a to the step-1 tissue mask, passed as an "
         "explicit path so a missing mask fails the job instead of silently using the full extent. "
         "--no_tissue_mask: decon uses otsu, 3a segments ungated (incompatible with --gpu_tile_size=auto).",
)
@click.option(
    "--tissue_rho_min",
    type=float,
    default=5.0,
    show_default=True,
    help="Min transcript density (tx/100 um^2) for step 1 to call a pixel tissue.",
)
@click.option(
    "--gpu_partition",
    default="gpu-l40s-ondemand",
    show_default=True,
    help="SLURM partition(s) for the GPU steps (2 decon, 3a cellpose), comma-separated. "
         "On-demand by default: spot reclaims killed jobs without requeueing them.",
)
@click.option(
    "--gpu_constraint",
    default="g6e.4xlarge|g6e.16xlarge",
    show_default=True,
    help="SLURM --constraint for the GPU steps. '|' means OR, so the job takes "
         "whichever instance type has capacity.",
)
@click.option(
    "--cpu_partition",
    default="cpu-ondemand",
    show_default=True,
    help="SLURM partition(s) for the CPU steps (1, 3b, 4, 5), comma-separated. "
         "On-demand by default: spot reclaims killed jobs without requeueing them.",
)
@click.option(
    "--cpu_constraint",
    default="cpu-32vcpu|cpu-48vcpu|cpu-64vcpu",
    show_default=True,
    help="SLURM --constraint for the CPU steps. Regions above ~150k cells need "
         "'cpu-64vcpu|cpu-96vcpu' (3b consensus clustering OOMs on 62 GB).",
)
def config_templates_aws(
    experiment_name: str,
    brain_region: str,
    region_name: str | None,
    lab: str | None,
    data_path: Path | None,
    output_dir: Path,
    template_dir: Path | None,
    s3_bucket: str,
    config_dir: Path,
    exp_n: str | None,
    reg_n: str | None,
    steps: list,
    proseg_cutoffs_path: str | None,
    overwrite: bool,
    tissue_raster: bool,
    mmc_reference: str | None,
    mmc_store_path: str | None,
    mmc_args: str,
    force_2d_segmentation: bool,
    decon_iter: int,
    decon_tilesize: int,
    decon_psigma: float | None,
    cellpose_prefix: str,
    proseg_prefix: str,
    seg_mode: str,
    cellpose_model: str,
    gpu_tile_size: str,
    gpu_tile_overlap: int,
    use_tissue_mask: bool,
    tissue_rho_min: float,
    gpu_partition: str,
    gpu_constraint: str,
    cpu_partition: str,
    cpu_constraint: str,
):
    """Generate AWS Slurm job scripts for EXPERIMENT_NAME."""

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if region_name is None and data_path is None:
        click.echo(
            "ERROR: provide --region_name for a single region, "
            "or --data_path to discover all region_* directories.",
            err=True,
        )
        sys.exit(1)

    # Caught at generation time, not just at runtime: deconwolf silently ignores
    # --psigma whenever a tilesize is set, so this combination would produce a script
    # that runs for hours and quietly does not do what was asked. Warning rather than
    # error, matching decon_image: dropping tilesize is the more dangerous fix, since
    # tilesize is what bounds VRAM.
    if decon_psigma is not None and decon_tilesize and decon_tilesize > 0:
        click.echo(
            f"  WARNING: --decon_psigma {decon_psigma} will be IGNORED because "
            f"--decon_tilesize is {decon_tilesize}. deconwolf drops psigma whenever "
            f"tiling is on. Pass --decon_tilesize 0 to use psigma (peak VRAM rises "
            f"from ~2.5 GB to ~12 GB at a 2960px tile).",
            err=True,
        )

    if template_dir is None:
        template_dir = _TEMPLATE_DIR

    if not template_dir.exists():
        click.echo(f"ERROR: template directory not found: {template_dir}", err=True)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Load brain-region config
    # ------------------------------------------------------------------
    config_file = config_dir / f"{brain_region}_config.json"
    if not config_file.exists():
        click.echo(
            f"ERROR: brain-region config not found: {config_file}\n"
            f"Create it from aws_templates/config.json and adjust paths for {brain_region}.",
            err=True,
        )
        sys.exit(1)

    with open(config_file) as f:
        config = json.load(f)

    ROOT_DIR = config.get("PROCESSED_ROOT_PATH")
    SEGMENTATION_DIR = config.get("SEGMENTATION_OUT_PATH")
    CUTOFFS_PATH = config.get("CUTOFFS_PATH")
    CONFIG_PATH = str(config_file)

    missing = [k for k, v in [
        ("PROCESSED_ROOT_PATH", ROOT_DIR),
        ("SEGMENTATION_OUT_PATH", SEGMENTATION_DIR),
        ("CUTOFFS_PATH", CUTOFFS_PATH),
    ] if v is None]
    if missing:
        click.echo(
            f"ERROR: {config_file} is missing required keys: {missing}",
            err=True,
        )
        sys.exit(1)
    if config.get("REGION_LABELS") is None:
        click.echo(
            f"  WARNING: {config_file} has no REGION_LABELS: donors become the region name without "
            "'region_' (e.g. region_UCI-5224 -> UCI-5224).",
        )
    mmc_store_path = mmc_store_path or config.get("MMC_DIR")
    if mmc_reference is not None and mmc_store_path is None:
        click.echo("ERROR: --mmc_reference needs --mmc_store_path or MMC_DIR in the config.", err=True)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Derive EXP_N
    # ------------------------------------------------------------------
    if exp_n is not None:
        pass  # user override
    elif lab is None:
        exp_n = experiment_name
    else:
        lab = lab.lower()
        if lab == "salk":
            exp_n = rename_exp_salk(experiment_name)
        elif lab == "ucsd":
            exp_n = rename_exp_ucsd(experiment_name)
        else:
            click.echo(f"ERROR: unknown lab '{lab}'. Use 'salk' or 'ucsd'.", err=True)
            sys.exit(1)

    # ------------------------------------------------------------------
    # Collect regions to process
    # ------------------------------------------------------------------
    if region_name is not None:
        region_names = [region_name]
    else:
        region_names = sorted(p.name for p in data_path.glob("region_*") if p.is_dir())
        if not region_names:
            click.echo(f"ERROR: no region_* directories found in {data_path}", err=True)
            sys.exit(1)
        click.echo(f"Discovered {len(region_names)} regions: {region_names}")

    # ------------------------------------------------------------------
    # Step flags
    # ------------------------------------------------------------------
    if len(steps) == 1 and steps[0] == "all":
        step1 = step2 = step3 = step4 = "true"
        step5 = "true" if mmc_reference else "false"
    else:
        steps = [x.strip().lower() for x in steps]
        step1 = "true" if "1" in steps else "false"
        step2 = "true" if "2" in steps else "false"
        step3 = "true" if "3" in steps else "false"
        step4 = "true" if "4" in steps else "false"
        step5 = "true" if "5" in steps else "false"
    if step5 == "true" and not mmc_reference:
        click.echo("ERROR: step 5 (MapMyCells) needs --mmc_reference.", err=True)
        sys.exit(1)
    overwrite_arg = "--overwrite" if overwrite else ""
    tissue_raster_arg = "" if tissue_raster else "--no_raster "

    # ------------------------------------------------------------------
    # Cellpose image extension
    # ------------------------------------------------------------------
    image_ext = "_z3.decon.tif" if force_2d_segmentation else ".decon.tif"
    if force_2d_segmentation:
        click.echo("  force_2d_segmentation=True → cellpose will use _z3.decon.tif")
        if seg_mode == "3d_true":
            click.echo("ERROR: --force_2d_segmentation gives cellpose one plane; use --seg_mode 2d.", err=True)
            sys.exit(1)

    # ------------------------------------------------------------------
    # Segmentation wiring for step 3a
    # ------------------------------------------------------------------
    # 'auto' optimises the tile grid against the tissue mask, so it cannot run
    # ungated. Catch that here rather than letting the GPU job start and die.
    tiling_off = str(gpu_tile_size).lower() in ("", "none", "off", "false")
    if not tiling_off and str(gpu_tile_size).lower() == "auto" and not use_tissue_mask:
        click.echo(
            "ERROR: --gpu_tile_size=auto needs a tissue mask to optimise the grid "
            "against. Either drop --no_tissue_mask, or give an explicit tile size.",
            err=True,
        )
        sys.exit(1)

    seg_extra = []
    if not tiling_off:
        seg_extra.append(f"--gpu_tile_size={gpu_tile_size}")
        seg_extra.append(f"--gpu_tile_overlap={gpu_tile_overlap}")
    seg_extra_args = " ".join(seg_extra)

    if use_tissue_mask and step1 == "false":
        click.echo(
            "  WARNING: step 1 builds tissue_mask.tif but is disabled — 3a will fail "
            "unless the mask already exists in S3 for this region.",
        )

    click.echo(
        f"  segmentation: seg_mode={seg_mode} model={cellpose_model} "
        f"tiling={'off' if tiling_off else gpu_tile_size} "
        f"tissue_mask={'on' if use_tissue_mask else 'off'}"
    )

    # ------------------------------------------------------------------
    # Process each region
    # ------------------------------------------------------------------
    for current_region in region_names:

        # Derive REG_N
        if reg_n is not None and len(region_names) == 1:
            current_reg_n = reg_n
        elif lab is None:
            current_reg_n = current_region
        elif lab == "salk":
            current_reg_n = rename_reg_salk(current_region)
        elif lab == "ucsd":
            current_reg_n = rename_reg_ucsd(current_region)
            click.echo(f"  Region {current_region} → {current_reg_n}")

        click.echo(f"\nConfiguring {exp_n} / {current_reg_n}  ({current_region})")

        # Explicit path rather than the backend's auto-discovery, so a region whose
        # mask was never built fails the job instead of quietly segmenting the full
        # extent -- a silent difference that is hard to spot downstream. Trailing
        # space: the template appends a line continuation right after it.
        tissue_mask_arg = (
            f"--tissue_mask_path {ROOT_DIR}/{experiment_name}/out/{current_region}"
            "/images/tissue_mask.tif "
            if use_tissue_mask else ""
        )

        # Deconvolution args, one `--transform_args KEY VALUE` per pair (the CLI option is
        # multiple=True). Trailing space: the template appends a line continuation.
        _decon_pairs = [("iter", decon_iter)]
        if decon_tilesize and decon_tilesize > 0:
            _decon_pairs.append(("tilesize", decon_tilesize))
        if decon_psigma is not None:
            _decon_pairs.append(("psigma", decon_psigma))
        decon_filter_args = "".join(
            f"--transform_args {k} {v} " for k, v in _decon_pairs)

        # Create job output directory
        job_dir = output_dir / brain_region / exp_n / current_reg_n
        job_dir.mkdir(parents=True, exist_ok=True)
        click.echo(f"  Job directory: {job_dir}")

        # Create log directory so SBATCH --output paths exist
        log_dir = _DEFAULT_LOG_DIR / brain_region / exp_n
        log_dir.mkdir(parents=True, exist_ok=True)

        decon_dir = f"{ROOT_DIR}/{experiment_name}/out/{current_region}/decon-images"

        # Fill and write each template
        for template_file in sorted(template_dir.glob("*.sh")):
            if template_file.name.startswith("5_") and not mmc_reference:
                continue  # no reference for this brain region
            out_path = job_dir / template_file.name
            with open(template_file) as fin:
                content = fin.read()

            try:
                content = content.format(
                    EXPERIMENT=experiment_name,
                    REGION=current_region,
                    BR=brain_region,
                    EXP_N=exp_n,
                    REG_N=current_reg_n,
                    ROOT_DIR=ROOT_DIR,
                    SEGMENTATION_DIR=SEGMENTATION_DIR,
                    CUTOFFS_PATH=CUTOFFS_PATH,
                    PROSEG_CUTOFFS_PATH=proseg_cutoffs_path or CUTOFFS_PATH,
                    CONFIG_PATH=CONFIG_PATH,
                    S3_BUCKET=s3_bucket,
                    STEP_1=step1,
                    STEP_2=step2,
                    STEP_3=step3,
                    STEP_4=step4,
                    STEP_5=step5,
                    OVERWRITE_ARG=overwrite_arg,
                    TISSUE_RASTER_ARG=tissue_raster_arg,
                    DECON_DIR=decon_dir,
                    MMC_REFERENCE=mmc_reference or "",
                    MMC_STORE_PATH=mmc_store_path or "",
                    MMC_ARGS=mmc_args,
                    NAMING_MAP_ARG="\\\n    --naming-map /home/ubuntu/aklein/site-images/naming_map.csv",
                    IMAGE_EXT=image_ext,
                    CELLPOSE_PREFIX=cellpose_prefix,
                    PROSEG_PREFIX=proseg_prefix,
                    SEG_MODE=seg_mode,
                    CELLPOSE_MODEL=cellpose_model,
                    SEG_EXTRA_ARGS=seg_extra_args,
                    TISSUE_MASK_ARG=tissue_mask_arg,
                    DECON_FILTER_ARGS=decon_filter_args,
                    DECON_MASK_ARG=(
                        f"--select_tiles tissue_mask --mask_path {ROOT_DIR}/{experiment_name}/out/"
                        f"{current_region}/images/tissue_mask.tif "
                        if use_tissue_mask else "--select_tiles otsu "),
                    TISSUE_RHO_MIN=tissue_rho_min,
                    GPU_PARTITION=gpu_partition,
                    GPU_CONSTRAINT=gpu_constraint,
                    CPU_PARTITION=cpu_partition,
                    CPU_CONSTRAINT=cpu_constraint,
                )
            except KeyError as e:
                click.echo(
                    f"  WARNING: unknown placeholder {e} in {template_file.name} — skipping.",
                    err=True,
                )
                continue

            with open(out_path, "w") as fout:
                fout.write(content)

            # Make all generated scripts executable
            out_path.chmod(out_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
            click.echo(f"  Written: {out_path.name}")

        click.echo(f"  Done — submit with: bash {job_dir}/chain.sh")

    click.echo(f"\nAll done. {len(region_names)} region(s) configured.")


if __name__ == "__main__":
    config_templates_aws()
