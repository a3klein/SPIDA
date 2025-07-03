#!/bin/bash


    # SEGMENTATION
    #  pixi run -e preprocessing python src/spida/cli.py io load_segmentation_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 /ceph/cephatlas/aklein/bican/data/segmented/202506151211_BICAN4x1-CAH-E03_VMSC31810/cellpose --plot=False --type=vpt --prefix_name=cellposeSAM

    # FILTERING 
    pixi run -e preprocessing python src/spida/cli.py P filter_cells_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 cellposeSAM --plot=True --cutoffs_path=None

    # SETUP
    pixi run -e preprocessing python src/spida/cli.py P setup_adata_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 cellposeSAM --plot=True
    
