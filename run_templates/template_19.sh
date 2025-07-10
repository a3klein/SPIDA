#!/bin/bash

    
    # DEFAULT 
    pixi run -e preprocessing python src/spida/cli.py io ingest_all 202506291134_BICAN-4x1-CAT-E-03_VMSC31910 --plot=False  # Plotting is raising errors? 

    # FILTERING 
    pixi run -e preprocessing python src/spida/cli.py P filter_cells_all 202506291134_BICAN-4x1-CAT-E-03_VMSC31910 default --plot=True --cutoffs_path=/home/x-aklein2/projects/aklein/BICAN/data/filtering_cutoffs.json # Set Cutoffs Path

    # SETUP ADATA 
    pixi run -e preprocessing python src/spida/cli.py P setup_adata_all 202506291134_BICAN-4x1-CAT-E-03_VMSC31910 default --plot=True
    
