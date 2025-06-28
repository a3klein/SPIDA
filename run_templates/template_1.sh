#!/bin/bash


        # SEGMENTATION
        pixi run -e cellpose python src/spida/segmentation/main.py segment_experiment cellpose 202506151211_BICAN4x1-CAH-E03_VMSC31810 --input_dir=None --output_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506151211_BICAN4x1-CAH-E03_VMSC31810/cellpose/
        pixi run -e preprocessing python src/spida/io/main.py load_segmentation_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 /ceph/cephatlas/aklein/bican/data/segmented/202506151211_BICAN4x1-CAH-E03_VMSC31810/cellpose/ --plot=False --type=vpt --prefix_name=cellposeSAM
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 cellposeSAM --plot=True --cutoffs_path=None

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 cellposeSAM --plot=True
        
