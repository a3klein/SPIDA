import os
import sys
from pathlib import Path


if len(sys.argv) != 4: 
    print("Usage: python script.py <experiment_name> <region_name> <lab(ucsd/salk)>")
    sys.exit(1)

experiment_name = sys.argv[1]
region_name = sys.argv[2]
lab = sys.argv[3].lower()

def rename_exp_salk(x): 
    exp_n = x.split("_")[1]
    temp = exp_n.split("-")[-2]
    if len(temp) == 1: temp = exp_n.split("-")[-3]
    return temp

def rename_reg_salk(x): 
    reg_n = x.split("_")[1]
    temp = reg_n.split("-")[0]
    if len(temp) == 3: temp = "".join(reg_n.split("-")[0:2])
    return temp

def rename_exp_ucsd(x): 
    exp_n = x.split("_")[1].split("BICAN")[-1]
    return exp_n

def rename_reg_ucsd(x): 
    reg_n = x.split("_")[1].split("Q0")[0]
    return reg_n


if lab == "salk": 
    exp_n = rename_exp_salk(experiment_name)
    reg_n = rename_reg_salk(region_name)
elif lab == "ucsd":
    exp_n = rename_exp_ucsd(experiment_name)
    reg_n = rename_reg_ucsd(region_name)
else:
    print("Invalid lab name. Use 'salk' or 'ucsd'.")
    sys.exit(1)
    

print(f"Experiment Name: {exp_n}")
print(f"Region Name: {reg_n}")

#SET OUTPUT DIRECTORY 
template_path = Path("/anvil/projects/x-mcb130189/aklein/SPIDA/slurm_scripts/templates")
output_dir = Path(f"/anvil/projects/x-mcb130189/aklein/BICAN/image_processing/{exp_n}/{reg_n}")
output_dir.mkdir(parents=True, exist_ok=True)
print("Output directory created at:", output_dir)

# configure_chain 
in_path = template_path / "chain.sh"
out_path = output_dir / "chain.sh"
with open(out_path, "w") as fout: 
    with open(in_path, 'r') as fin:
        cont=fin.read()
        cont = cont.format(EXP_N=exp_n, REG_N=reg_n)
    fout.write(cont)

# configure scripts: 
for proc_file in ["1_start.sh", "2_decon_DAPI.sh", "2_decon_PolyT.sh", "3_cellpose.sh", "4_proseg.sh", "4_proseg_2.sh", "5_setup.sh", "5_setup_2.sh"]: 
    in_path = template_path / proc_file
    out_path = output_dir / proc_file
    with open(out_path, "w") as fout: 
        with open(in_path, 'r') as fin:
            cont=fin.read().format(EXP_N=exp_n, REG_N=reg_n, EXPERIMENT=experiment_name, REGION=region_name)
        fout.write(cont)

print("configured scripts")