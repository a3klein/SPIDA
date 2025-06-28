#!/bin/bash


        # SEGMENTATION
        pixi run -e cellpose python src/spida/segmentation/main.py segment_experiment cellpose 202506171319_BICAN-4x1-GP-E-05_VMSC31910 --input_dir=None --output_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/cellpose/
        pixi run -e preprocessing python src/spida/io/main.py load_segmentation_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 /ceph/cephatlas/aklein/bican/data/segmented/202506171319_BICAN-4x1-GP-E-05_VMSC31910/cellpose/ --plot=False --type=vpt --prefix_name=cellposeSAM
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 cellposeSAM --plot=True --cutoffs_path=None

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506171319_BICAN-4x1-GP-E-05_VMSC31910 cellposeSAM --plot=True
        
