#!/bin/bash


        # SEGMENTATION
        pixi run -e preprocessing python src/spida/segmentation/main.py segment_experiment proseg 202506151211_BICAN4x1-CAH-E03_VMSC31810 --input_dir=None --output_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506151211_BICAN4x1-CAH-E03_VMSC31810/proseg/
        pixi run -e preprocessing python src/spida/io/main.py load_segmentation_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 /ceph/cephatlas/aklein/bican/data/segmented/202506151211_BICAN4x1-CAH-E03_VMSC31810/proseg/ --plot=False --type=proseg --prefix_name=proseg
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 proseg --plot=True --cutoffs_path=None

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506151211_BICAN4x1-CAH-E03_VMSC31810 proseg --plot=True
        
