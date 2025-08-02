import os
import sys
from pathlib import Path
import click 


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
@click.argument("region_name")
@click.option("--lab", default="salk", type=click.Choice(["salk", "ucsd"], case_sensitive=False), help="Lab name (salk/ucsd)")
@click.option(
    "--output_dir",
    type=click.Path(exists=False, dir_okay=True, path_type=Path),
    default="/anvil/projects/x-mcb130189/aklein/BICAN/image_processing",
    help="Output directory for the scripts")
@click.option(
    "--template_dir",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    default="/anvil/projects/x-mcb130189/aklein/BICAN/image_processing",
    help="Template directory for the scripts",
)
@click.option("--reg_n", default=None, help="Shortened Region name to be used in the scripts (derived from region_name if None)")
@click.option("--exp_n", default=None, help="Shortened Experiment name to be used in the scripts (derived from experiment_name if None)")
def config_templates(
    experiment_name : str,
    region_name : str,
    lab : str, 
    output_dir : str = None,
    template_dir : str = None,
    reg_n : str = None,
    exp_n : str = None
):
    """ Configures the templates in the current directory for the provided EXPERIMENT_NAME and REGION_NAME, ran at the given LAB """

    lab = lab.lower()

    if exp_n is None: 
        if lab == "salk":
            exp_n = rename_exp_salk(experiment_name)
        elif lab == "ucsd":
            exp_n = rename_exp_ucsd(experiment_name)
        else:
            print("Invalid lab name. Use 'salk' or 'ucsd'.")
            sys.exit(1)
    else:
        exp_n = exp_n.replace(" ", "_").replace("-", "_").replace(".", "_")
    
    if reg_n is None:
        if lab == "salk":
            reg_n = rename_reg_salk(region_name)
        elif lab == "ucsd":
            reg_n = rename_reg_ucsd(region_name)
        else:
            print("Invalid lab name. Use 'salk' or 'ucsd'.")
            sys.exit(1)
    else:
        reg_n = reg_n.replace(" ", "_").replace("-", "_").replace(".", "_")

    click.echo(f"Experiment Name: {exp_n}")
    click.echo(f"Region Name: {reg_n}")

    # SET OUTPUT DIRECTORY
    # template_path = Path("/anvil/projects/x-mcb130189/aklein/SPIDA/slurm_scripts/templates")
    # output_dir = Path(
    #     f"/anvil/projects/x-mcb130189/aklein/BICAN/image_processing/{exp_n}/{reg_n}"
    # )
    output_dir = output_dir / exp_n / reg_n
    output_dir.mkdir(parents=True, exist_ok=True)
    click.echo(f"Output directory created at: {output_dir}")

    # configure_chain
    in_path = template_dir / "chain.sh"
    out_path = output_dir / "chain.sh"
    with open(out_path, "w") as fout:
        with open(in_path, "r") as fin:
            cont = fin.read()
            cont = cont.format(EXP_N=exp_n, REG_N=reg_n)
        fout.write(cont)

    # configure scripts:
    for proc_file in [
        "1_start.sh",
        "2_decon_DAPI.sh",
        "2_decon_PolyT.sh",
        "3_cellpose.sh",
        "4_proseg.sh",
        "4_proseg_2.sh",
        "5_setup.sh",
        "5_setup_2.sh",
    ]:
        in_path = template_dir / proc_file
        out_path = output_dir / proc_file
        with open(out_path, "w") as fout:
            with open(in_path, "r") as fin:
                cont = fin.read().format(
                    EXP_N=exp_n, REG_N=reg_n, EXPERIMENT=experiment_name, REGION=region_name
                )
            fout.write(cont)

    click.echo("configured scripts")

if __name__ == "__main__": 
    config_templates()