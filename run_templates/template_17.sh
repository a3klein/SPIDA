#!/bin/bash


        # SEGMENTATION
        pixi run -e preprocessing python src/spida/segmentation/main.py segment_experiment proseg 202506221112_BICAN4x1-NAC-E05_VMSC31810 --input_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506221112_BICAN4x1-NAC-E05_VMSC31810/cellpose/ --output_dir=/ceph/cephatlas/aklein/bican/data/segmented/202506221112_BICAN4x1-NAC-E05_VMSC31810/proseg_SAM/
        pixi run -e preprocessing python src/spida/io/main.py load_segmentation_all 202506221112_BICAN4x1-NAC-E05_VMSC31810 /ceph/cephatlas/aklein/bican/data/segmented/202506221112_BICAN4x1-NAC-E05_VMSC31810/proseg_SAM/ --plot=False --type=proseg --prefix_name=proseg_SAM
        
        # FILTERING 
        pixi run -e preprocessing python src/spida/P/main.py filter_cells_all 202506221112_BICAN4x1-NAC-E05_VMSC31810 proseg_SAM --plot=False --cutoffs_path=None

        # SETUP ADATA 
        pixi run -e preprocessing python src/spida/P/main.py setup_adata_all 202506221112_BICAN4x1-NAC-E05_VMSC31810 proseg_SAM --plot=False
        
