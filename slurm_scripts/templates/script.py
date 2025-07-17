import os
import sys
from pathlib import Path


if len(sys.argv) != 3: 
    print("Usage: python script.py <experiment_name> <region_name>")
    sys.exit(1)

experiment_name = sys.argv[1]
region_name = sys.argv[2]

exp_n = experiment_name.split("_")[1]
temp = exp_n.split("-")[-2]
if len(temp) == 1: temp = exp_n.split("-")[-3]
exp_n = temp

reg_n = region_name.split("_")[1]
temp = reg_n.split("-")[0]
if len(temp) == 3: temp = "".join(reg_n.split("-")[0:2])
reg_n = temp

print(f"Experiment Name: {exp_n}")
print(f"Region Name: {reg_n}")

#SET OUTPUT DIRECTORY 
template_path = Path("/anvil/projects/x-mcb130189/aklein/BICAN/image_processing/templates")
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
for proc_file in ["1_start.sh", "2_decon_DAPI.sh", "2_decon_PolyT.sh", "3_cellpose.sh", "4_proseg.sh", "5_setup.sh"]: 
    in_path = template_path / proc_file
    out_path = output_dir / proc_file
    with open(out_path, "w") as fout: 
        with open(in_path, 'r') as fin:
            cont=fin.read().format(EXP_N=exp_n, REG_N=reg_n, EXPERIMENT=experiment_name, REGION=region_name)
        fout.write(cont)

print("configures scripts")