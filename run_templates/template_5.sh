#!/bin/bash


        # SEGMENTATION
        pixi run -e preprocessing python src/spida/segmentation/main.py segment_experiment proseg 202506171319_BICAN-4x1-GP-E-05_VMSC31910 --input_dir=None --output_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/proseg/
        pixi run -e preprocessing python src/spida/io/main.py load_segmentation_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 /ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/proseg/ --plot=False --type=proseg --prefix_name=proseg
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 proseg --plot=True --cutoffs_path=None

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 proseg --plot=True
        
