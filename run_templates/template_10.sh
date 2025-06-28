#!/bin/bash

    
    # DEFAULT 
    # pixi run ingest_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 True 
    pixi run -e preprocessing python src/spida/io/main.py ingest_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 --plot=True

    # FILTERING 
    # pixi run filter_cells_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 default True None
    pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 default --plot=True --cutoffs_path=None
    
    # SETUP ADATA 
    # pixi run setup_adata_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 default True
    pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506221112_BICAN-4x1-MGM1-E-03_VMSC31910 default --plot=True
    
