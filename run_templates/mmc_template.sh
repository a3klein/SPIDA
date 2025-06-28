# Running mmc annotations on all the prefixes for the A8 experiment
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 default CORTEX D2M2006
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 proseg CORTEX D2M2006
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 cellpose_nuclei CORTEX D2M2006
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 proseg_nuclei CORTEX D2M2006
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 cellposeSAM CORTEX D2M2006
pixi run -e mapmycells python src/spida/I/main.py mmc_annotation_experiment 202505161056_BICAN-4x1-A8-E-03_VMSC31910 mesmer CORTEX D2M2006
