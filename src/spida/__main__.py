from pathlib import Path
import fire

# spida_path = "/ceph/cephatlas/aklein/spida/src"
# sys.path.append(spida_path)
# from spida._templates import *
from _templates import *

TEMPLATE_PREFIX = "template"

def _get_outfile(dir:str): 
    path = Path(dir)
    if not path.exists(): 
        path.mkdir(parents=True)
        return f"{dir}/{TEMPLATE_PREFIX}_1.sh"

    existing_files = list(path.glob(f"{TEMPLATE_PREFIX}*.sh"))
    numbers = []
    for file in existing_files: 
        try: 
            num = int(file.stem.split("_")[-1])
            numbers.append(num)
        except ValueError: 
            continue
    
    return f"{dir}/{TEMPLATE_PREFIX}_{max(numbers, default=0) + 1}.sh"


def _generate_run_command(type:str,
                          exp_name:str,
                          prefix_name:str,
                          **kwargs):
    """
    Generate the run command for a specific type of processing step in SPIDA.

    """
    
    seg_input_dir = kwargs.get("seg_input_dir", None)
    seg_output_dir = kwargs.get("seg_output_dir", None)
    seg_type = kwargs.get("seg_type", "")
    seg_fam = kwargs.get("seg_fam", "")
    plot_load = kwargs.get("plot_load", False)
    plot_filter = kwargs.get("plot_filtering", False)
    plot_setup = kwargs.get("plot_setup", False)
    cutoffs_path = kwargs.get("cutoffs_path", None)
    brain_region = kwargs.get("brain_region", None)
    codebook = kwargs.get("codebook", None)
    adata_path = kwargs.get("adata_path", None)
    adata_ref_path = kwargs.get("adata_ref_path", None)
    rna_cluster_col = kwargs.get("rna_cluster_col", "supercluster_name")


    if type == "start": 
        template = load_def_template()
        fout = template.format(
            EXPERIMENT=exp_name,
            PREFIX=prefix_name,
            LOAD_PLOT=plot_load,
            PLOT_FILTERING=plot_filter,
            PLOT_ADATA=plot_setup,
            CUTOFF_JSON_PATH=cutoffs_path,
        )
    elif type == "load_seg": 
        template = load_seg_template()
        fout = template.format(
            EXPERIMENT=exp_name,
            SEG_OUT_DIR=seg_output_dir,
            SEG_LOAD_PLOT=plot_load,
            SEG_FAM=seg_fam,
            PREFIX=prefix_name,
            PLOT_FILTERING=plot_filter,
            CUTOFF_JSON_PATH=cutoffs_path,
            PLOT_ADATA=plot_setup,
        )
    elif type == "seg":
        template = run_seg_template(seg_type)
        fout = template.format(
            SEG_TYPE=seg_type,
            EXPERIMENT=exp_name,
            SEG_IN_DIR=seg_input_dir,
            SEG_OUT_DIR=seg_output_dir,
            SEG_LOAD_PLOT=plot_load,
            SEG_FAM=seg_fam,
            PREFIX=prefix_name,
            PLOT_FILTERING=plot_filter,
            CUTOFF_JSON_PATH=cutoffs_path,
            PLOT_ADATA=plot_setup,
        )
    elif type == "mmc_annotation": 
        template = annot_mmc_template()
        fout = template.format(
            EXPERIMENT=exp_name,
            PREFIX=prefix_name,
            BRAIN_REGION=brain_region,
            CODEBOOK=codebook,
            ADATA_PATH=adata_path,
        )
    elif type == "allc_integration": 
        template = annot_allc_template()
        fout = template.format(
            EXPERIMENT=exp_name,
            PREFIX=prefix_name,
            ADATA_REF_PATH=adata_ref_path,
            RNA_CLUSTER_COLUMN=rna_cluster_col,
            ADATA_PATH=adata_path,
        )
    return fout

def gen_run(out_loc:str,
            type:str | list,
            exp_name:str | list, 
            prefix_name:str | list, 
            **kwargs
            ):
    """
    Generate the run commands for the main processing steps in SPIDA.
    """
    
    # Get the output file location 
    outfile = _get_outfile(out_loc)
    print(outfile)

    if isinstance(type, str): 
        type = [type]
    if isinstance(prefix_name, str):
        prefix_name = [prefix_name]
    if isinstance(exp_name, str):
        exp_name = [exp_name]
    
    
    output = "#!/bin/bash\n\n"
    for _e in exp_name: 
        for _p in prefix_name:
            for _t in type: 
                fout = _generate_run_command(_t, _e, _p, **kwargs)
                output += fout + "\n"
    
    with open(outfile, "w") as f:
        f.write(output)
    
    
if __name__ == "__main__": 
    fire.Fire({"run_template" : gen_run})

