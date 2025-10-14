import os
import sys
from dotenv import load_dotenv
from pathlib import Path
import click 
import json

load_dotenv()


def rename_exp_salk(x):
    exp_n = x.split("_")[1]
    temp = exp_n.split("-")[-2]
    if len(temp) == 1:
        temp = exp_n.split("-")[-3]
    return temp


def rename_reg_salk(x):
    reg_n = x.split("_")[1]
    temp = reg_n.split("-")[0]
    if len(temp) == 3:
        temp = "".join(reg_n.split("-")[0:2])
    return temp


def rename_exp_ucsd(x):
    exp_n = x.split("_")[1].split("BICAN")[-1]
    return exp_n


def rename_reg_ucsd(x):
    reg_n = x.split("_")[1].split("Q0")[0]
    return reg_n



@click.command("config-template")
@click.argument("experiment_name")
@click.option("--region_name", default=None, help="Region name (if None, generate all regions in an experiment)")
@click.option(
    "--lab",
    default=None,
    type=click.Choice([None, "salk", "ucsd"], case_sensitive=False),
    help="Lab name (salk/ucsd) (default: None, do no renaming)")
@click.option(
    "data_path",
    "--data_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default=None,
    help="If region_name is None, provide the path to the experiment directory to generate all regions")
@click.option(
    "--output_dir",
    type=click.Path(exists=False, dir_okay=True, path_type=Path),
    default="/anvil/projects/x-mcb130189/aklein/BICAN/image_processing",
    help="Output directory for the scripts")
@click.option(
    "--template_dir",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default="/home/x-aklein2/projects/aklein/SPIDA/slurm_scripts/templates",
    help="Template directory for the scripts",
)
@click.option("--reg_n", default=None, help="Shortened Region name to be used in the scripts (derived from region_name if None)")
@click.option("--exp_n", default=None, help="Shortened Experiment name to be used in the scripts (derived from experiment_name if None)")
@click.option(
    "--config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Path to the config file for the project (if None, use environment variables)",)
def config_templates(
    experiment_name : str,
    region_name : str | None = None,
    lab : str = None, 
    data_path : str | Path = None,
    output_dir : str = None,
    template_dir : str = None,
    reg_n : str = None,
    exp_n : str = None,
    config_path : str | Path = None,
):
    """ Configures the templates in the current directory for the provided EXPERIMENT_NAME and REGION_NAME, ran at the given LAB """
    click.echo("Configuring templates...")
    click.echo(f"Experiment Name: {experiment_name}")
    click.echo(f"Region Name: {region_name if region_name is not None else 'All regions in the experiment'}")
    click.echo(f"Lab: {lab if lab is not None else 'No renaming'}")
    click.echo(f"Data Path: {data_path if data_path is not None else 'Not provided'}")
    click.echo(f"Output Directory: {output_dir}")
    click.echo(f"Template Directory: {template_dir}")
    click.echo(f"Config Path: {config_path if config_path is not None else 'Not provided, using environment variables'}")
    
    if region_name is None and data_path is None: 
        print("If region_name is None, provide the path to the experiment directory to generate all regions")
        sys.exit(1)
    if lab is None: 
        exp_n = exp_n if exp_n is not None else experiment_name
    else: 
        lab = lab.lower()
        if lab == "salk":
            exp_n = rename_exp_salk(experiment_name)
        elif lab == "ucsd":
            exp_n = rename_exp_ucsd(experiment_name)
        else:
            print("Invalid lab name. Use 'salk' or 'ucsd'.")
            sys.exit(1)
    
    if region_name is None:
        data_path = Path(data_path)
        if not data_path.exists():
            print(f"Data path {data_path} does not exist")
            sys.exit(1)
        region_names = sorted([f.name for f in data_path.glob("region_*") if f.is_dir()])
        print(region_names)
    else: 
        region_names = [region_name]

    for region_name in region_names:
        if lab is None:
            reg_n = reg_n if reg_n is not None else region_name
        else:
            if lab == "salk":
                reg_n = rename_reg_salk(region_name)
            elif lab == "ucsd":
                reg_n = rename_reg_ucsd(region_name)
                click.echo(f"Using region {region_name} with shortened name {reg_n}")
            else:
                print("Invalid lab name. Use 'salk' or 'ucsd'.")
                sys.exit(1)

        click.echo(f"Experiment Name: {exp_n}")
        click.echo(f"Region Name: {reg_n}")

        if config_path is not None: 
            with open(config_path, "r") as f:
                data = json.load(f)
                ROOT_PATH = data.get("PROCESSED_ROOT_PATH", None)
                SEGMENTATION_DIR = data.get("SEGMENTATION_OUT_PATH", None)
                CUTOFFS_PATH = data.get("CUTOFFS_PATH", None)
                CONFIG=f"\\\n\t--config {config_path}" 
            if ROOT_PATH is None or SEGMENTATION_DIR is None or CUTOFFS_PATH is None:
                raise ValueError("Config file must contain ROOT_PATH, SEGMENTATION_DIR, and CUTOFFS_PATH")
        else: # Defaults
            ROOT_PATH = os.getenv("ROOT_PATH", "/anvil/projects/x-mcb130189/aklein/BICAN")
            SEGMENTATION_DIR = os.getenv("SEGMENTATION_DIR", "/anvil/projects/x-mcb130189/aklein/BICAN/data/segmentated")
            CUTOFFS_PATH = os.getenv("CUTOFFS_PATH", "/home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json")

        output_dir_2 = output_dir / exp_n / reg_n
        output_dir_2.mkdir(parents=True, exist_ok=True)
        click.echo(f"Output directory created at: {output_dir_2}")

        # configure_chain
        # in_path = template_dir / "chain.sh"
        # out_path = output_dir / "chain.sh"
        # with open(out_path, "w") as fout:
        #     with open(in_path, "r") as fin:
        #         cont = fin.read()
        #         cont = cont.format(OUTPUT_DIR=output_dir, EXP_N=exp_n, REG_N=reg_n)
        #     fout.write(cont)

        # configure scripts:
        for proc_file in template_dir.glob("*.sh"):
            # if proc_file.name == "chain.sh":
            #     continue
            print(f"Configuring {proc_file.name}")
            out_path = output_dir_2 / proc_file.name
            with open(out_path, "w") as fout:
                with open(proc_file, "r") as fin:
                    cont = fin.read().format(
                        EXP_N=exp_n,
                        REG_N=reg_n,
                        EXPERIMENT=experiment_name,
                        REGION=region_name,
                        ROOT_PATH=ROOT_PATH,
                        SEGMENTATION_DIR=SEGMENTATION_DIR,
                        CUTOFFS_PATH=CUTOFFS_PATH,
                        CONFIG=CONFIG if config_path is not None else "",
                        OUTPUT_DIR=output_dir_2,
                    )
                fout.write(cont)

        click.echo("configured scripts")

if __name__ == "__main__": 
    config_templates()

# for proc_file in [
#     "1_start.sh",
#     "2_decon_DAPI.sh",
#     "2_decon_PolyT.sh",
#     "3_cellpose.sh",
#     "4_proseg.sh",
#     "4_proseg_2.sh",
#     "5_setup.sh",
#     "5_setup_2.sh",
#     "5_setup_mapped.sh",
#     "5_setup_aligned.sh",
# ]:
# in_path = template_dir / proc_file
# out_path = output_dir / proc_file