def load_def_template():
    template = """    
    # DEFAULT 
    # pixi run ingest_all {EXPERIMENT} {LOAD_PLOT} 
    pixi run -e preprocessing python src/spida/cli.py io ingest_all {EXPERIMENT} --plot={LOAD_PLOT}

    # FILTERING 
    # pixi run filter_cells_all {EXPERIMENT} {PREFIX} {PLOT_FILTERING} {CUTOFF_JSON_PATH}
    pixi run -e preprocessing python src/spida/cli.py P filter_cells_all {EXPERIMENT} {PREFIX} --plot={PLOT_FILTERING} --cutoffs_path={CUTOFF_JSON_PATH}

    # SETUP ADATA 
    # pixi run setup_adata_all {EXPERIMENT} {PREFIX} {PLOT_ADATA}
    pixi run -e preprocessing python src/spida/cli.py P setup_adata_all {EXPERIMENT} {PREFIX} --plot={PLOT_ADATA}
    """
    return template


def load_seg_template():
    template = """
    # SEGMENTATION
    # pixi run load_seg_all {EXPERIMENT} {SEG_OUT_DIR} {SEG_LOAD_PLOT} {SEG_FAM} {PREFIX}
    pixi run -e preprocessing python src/spida/cli.py io load_segmentation_all {EXPERIMENT} {SEG_OUT_DIR} --plot={SEG_LOAD_PLOT} --type={SEG_FAM} --prefix_name={PREFIX}

    # FILTERING 
    # pixi run filter_cells_all {EXPERIMENT} {PREFIX} {PLOT_FILTERING} {CUTOFF_JSON_PATH}
    pixi run -e preprocessing python src/spida/cli.py P filter_cells_all {EXPERIMENT} {PREFIX} --plot={PLOT_FILTERING} --cutoffs_path={CUTOFF_JSON_PATH}

    # SETUP ADATA 
    # pixi run setup_adata_all {EXPERIMENT} {PREFIX} {PLOT_ADATA}
    pixi run -e preprocessing python src/spida/cli.py P setup_adata_all {EXPERIMENT} {PREFIX} --plot={PLOT_ADATA}
    """
    return template


def run_seg_template(seg_type: str):
    if seg_type == "cellpose":
        template = """
        # SEGMENTATION
        pixi run -e cellpose python src/spida/segmentation/main.py segment_experiment cellpose {EXPERIMENT} --input_dir={SEG_IN_DIR} --output_dir={SEG_OUT_DIR}
        pixi run -e preprocessing python src/spida/cli.py io load_segmentation_all {EXPERIMENT} {SEG_OUT_DIR} --plot={SEG_LOAD_PLOT} --type={SEG_FAM} --prefix_name={PREFIX}
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/cli.py P filter_cells_all {EXPERIMENT} {PREFIX} --plot={PLOT_FILTERING} --cutoffs_path={CUTOFF_JSON_PATH}

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/cli.py P setup_adata_all {EXPERIMENT} {PREFIX} --plot={PLOT_ADATA}
        """
    elif seg_type == "mesmer":
        template = """
        # SEGMENTATION
        pixi run -e cellpose python src/spida/segmentation/main.py segment_experiment mesmer {EXPERIMENT} --input_dir={SEG_IN_DIR} --output_dir={SEG_OUT_DIR}
        pixi run -e preprocessing python src/spida/cli.py io load_segmentation_all {EXPERIMENT} {SEG_OUT_DIR} --plot={SEG_LOAD_PLOT} --type={SEG_FAM} --prefix_name={PREFIX}

        # FILTERING
        pixi run -e preprocessing python src/spida/cli.py P filter_cells_all {EXPERIMENT} {PREFIX} --plot={PLOT_FILTERING} --cutoffs_path={CUTOFF_JSON_PATH}

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/cli.py P setup_adata_all {EXPERIMENT} {PREFIX} --plot={PLOT_ADATA}
        """
    else:
        template = """
        # SEGMENTATION
        pixi run -e preprocessing python src/spida/segmentation/main.py segment_experiment {SEG_TYPE} {EXPERIMENT} --input_dir={SEG_IN_DIR} --output_dir={SEG_OUT_DIR}
        pixi run -e preprocessing python src/spida/cli.py io load_segmentation_all {EXPERIMENT} {SEG_OUT_DIR} --plot={SEG_LOAD_PLOT} --type={SEG_FAM} --prefix_name={PREFIX}

        # FILTERING
        pixi run -e preprocessing python src/spida/cli.py P filter_cells_all {EXPERIMENT} {PREFIX} --plot={PLOT_FILTERING} --cutoffs_path={CUTOFF_JSON_PATH}

        # SETUP ADATA
        pixi run -e preprocessing python src/spida/cli.py P setup_adata_all {EXPERIMENT} {PREFIX} --plot={PLOT_ADATA}
        """
    return template


def seg_template_P():
    template = """
    # SEGMENTATION
    pixi run segment_all {SEG_TYPE} {EXPERIMENT}
    pixi run load_seg_all {EXPERIMENT} {SEG_DIR} {SEG_LOAD_PLOT} {SEG_FAM} {PREFIX}

    # FILTERING 
    pixi run filter_cells_all {EXPERIMENT} {PREFIX} {PLOT_FILTERING} {CUTOFF_JSON_PATH}

    # SETUP ADATA 
    pixi run setup_adata_all {EXPERIMENT} {PREFIX} {PLOT_ADATA}
    """
    return template


def annot_mmc_template():
    template = """
    # ANNOTATION
    pixi run -e mapmycells python src/spida/cli.py I mmc_annotation_experiment {EXPERIMENT} {PREFIX} {BRAIN_REGION} {CODEBOOK}
    pixi run -e preprocessing python src/spida/cli.py I backup_adata_experiment {EXPERIMENT} {PREFIX} {ADATA_PATH}
    """
    return template


def annot_allc_template():
    template = """
    # ANNOTATION
    pixi run -e preprocessing python src/spida/cli.py I allcools_integration_experiment {EXPERIMENT} {PREFIX} {ADATA_REF_PATH} --rna_cell_type_column={RNA_CLUSTER_COLUMN}
    pixi run -e preprocessing python src/spida/cli.py I backup_adata_experiment {EXPERIMENT} {PREFIX} {ADATA_PATH}
    """
    return template


### template sandbox:

setup_mmc_template = """
pixi run -e mapmycells python src/spida/I/main.py setup_mmc 
/ceph/cephatlas/aklein/bican/reference/cortex/RNA_ren/BICAN_all.h5ad 
[supercluster_name,cluster_name,subcluster_name] 
CORTEX 
D2M2006
"""

vpt_template = """
"""

proseg_template = """
"""


seg_template = """
# SEGMENTATION
pixi run segment_{SEG_TYPE}_all {EXPERIMENT}
pixi run load_seg_all {EXPERIMENT} {SEG_DIR} {SEG_LOAD_PLOT} {SEG_FAM} {PREFIX}

# FILTERING 
pixi run filter_cells_all {EXPERIMENT} {PREFIX} {PLOT_FILTERING} {CUTOFF_JSON_PATH}

# SETUP ADATA 
pixi run setup_adata_all {EXPERIMENT} {PREFIX} {PLOT_ADATA}
"""


seg_template_p = """
# SEGMENTATION
pixi run segment_all {SEG_TYPE} {EXPERIMENT}
pixi run load_seg_all {EXPERIMENT} {SEG_DIR} {SEG_LOAD_PLOT} {SEG_FAM} {PREFIX}

# FILTERING 
pixi run filter_cells_all {EXPERIMENT} {PREFIX} {PLOT_FILTERING} {CUTOFF_JSON_PATH}

# SETUP ADATA 
pixi run setup_adata_all {EXPERIMENT} {PREFIX} {PLOT_ADATA}
"""
