"""
script_aws.py — Generate AWS ParallelCluster Slurm job scripts for the SPIDA pipeline.

Each (EXPERIMENT, REGION, BRAIN_REGION) triple gets its own directory under:
    {output_dir}/{BR}/{EXP_N}/{REG_N}/

containing populated versions of every *.sh template in --template_dir, plus a
chain.sh that chains them with sbatch --dependency.

The brain-region config file ({config_dir}/{BR}_config.json) is required and
provides ROOT_DIR, SEGMENTATION_DIR, and CUTOFFS_PATH.  It is also passed
verbatim as CONFIG_PATH so SPIDA commands receive --config at runtime.

Usage examples
--------------
# Single region:
python script_aws.py 202505231106_BICAN-4x1-MTC-E-05_VMSC31810 \\
    --brain_region CTX --region_name region_UCI-5224 --lab salk

# All regions in an experiment (discovers region_* dirs from a local path):
python script_aws.py 202505231106_BICAN-4x1-MTC-E-05_VMSC31810 \\
    --brain_region CTX \\
    --data_path /home/ubuntu/aklein/raw/202505231106_BICAN-4x1-MTC-E-05_VMSC31810/out \\
    --lab salk

# Re-run only steps 3 and 4:
python script_aws.py ... --steps 3,4
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
_DEFAULT_TEMPLATE_DIR = _SCRIPT_DIR / "aws_templates"
_DEFAULT_OUTPUT_DIR = Path("/home/ubuntu/aklein/slurm_jobs")
_DEFAULT_CONFIG_DIR = Path("/home/ubuntu/aklein/spida_config")
_DEFAULT_S3_BUCKET = "s3://salk-workstation-data-dev-020125249408"
_DEFAULT_LOG_DIR = Path("/home/ubuntu/aklein/spida_logs")


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
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default=_DEFAULT_TEMPLATE_DIR,
    show_default=True,
    help="Directory containing the *.sh template files. "
         "Defaults to slurm_scripts/aws_templates/ next to this script.",
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
    "--s3_bucket",
    default=_DEFAULT_S3_BUCKET,
    show_default=True,
    help="S3 bucket base URL (no trailing slash).",
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
    help="Comma-separated list of steps to enable in chain.sh (1-4). "
         "Default 'all' enables every step.",
)
def config_templates_aws(
    experiment_name: str,
    brain_region: str,
    region_name: str | None,
    lab: str | None,
    data_path: Path | None,
    output_dir: Path,
    template_dir: Path,
    config_dir: Path,
    s3_bucket: str,
    exp_n: str | None,
    reg_n: str | None,
    steps: list,
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

    if not template_dir.exists():
        click.echo(
            f"ERROR: template directory not found: {template_dir}\n"
            "Copy the aws/template/ files to that location or pass --template_dir.",
            err=True,
        )
        sys.exit(1)

    # ------------------------------------------------------------------
    # Load brain-region config
    # ------------------------------------------------------------------
    config_file = config_dir / f"{brain_region}_config.json"
    if not config_file.exists():
        click.echo(
            f"ERROR: brain-region config not found: {config_file}\n"
            f"Create it from aws/template/config.json and adjust paths for {brain_region}.",
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
        step1, step2, step3, step4 = "true", "true", "true", "true"
    else:
        step1 = "true" if "1" in steps else "false"
        step2 = "true" if "2" in steps else "false"
        step3 = "true" if "3" in steps else "false"
        step4 = "true" if "4" in steps else "false"

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

        # Create job output directory
        job_dir = output_dir / brain_region / exp_n / current_reg_n
        job_dir.mkdir(parents=True, exist_ok=True)
        click.echo(f"  Job directory: {job_dir}")

        # Create log directory so SBATCH --output paths exist
        log_dir = _DEFAULT_LOG_DIR / brain_region / exp_n
        log_dir.mkdir(parents=True, exist_ok=True)

        # Fill and write each template
        for template_file in sorted(template_dir.glob("*.sh")):
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
                    S3_BUCKET=s3_bucket,
                    ROOT_DIR=ROOT_DIR,
                    SEGMENTATION_DIR=SEGMENTATION_DIR,
                    CUTOFFS_PATH=CUTOFFS_PATH,
                    CONFIG_PATH=CONFIG_PATH,
                    STEP_1=step1,
                    STEP_2=step2,
                    STEP_3=step3,
                    STEP_4=step4,
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
