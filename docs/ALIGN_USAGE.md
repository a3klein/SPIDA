# Segmentation Alignment Module Usage Guide

## Overview

The `align_segmentations()` function in `align.py` aligns two sets of cell segmentations by their spatial overlap. This is useful when you have segmentations from different sources (e.g., nuclear segmentation and cellular segmentation) and want to create a unified dataset containing only cells that exist in both datasets.

## Key Features

- **One-to-one mapping**: Each cell in the first dataset is matched to at most one cell in the second dataset, and vice versa
- **Flexible geometry selection**: Choose which geometry to keep (larger cell, first prefix, second prefix, or intersection)
- **Parquet output**: Results are saved as parquet files compatible with the VPT segmentation pipeline
- **Spatial join based**: Uses geopandas spatial joins for efficient alignment

## Basic Usage

```python
from spida.S.segmentation import align_segmentations

result = align_segmentations(
    zarr_path="/path/to/zarr/store",
    exp_name="202511071313_BICAN-4x1-HIPU-E-01_VMSC31910",
    reg_name="region_UCI-5224",
    prefix1="cellpose_nuc",
    prefix2="cellpose_cell",
    geometry_mode="larger"  # Default mode
)

print(f"Aligned {result['n_aligned']} cell pairs")
print(f"Output: {result['output_path']}")
```

## Parameters

### Required
- **zarr_path** (str): Path to the spatialdata zarr store containing the segmentations
- **exp_name** (str): Experiment name (e.g., "202511071313_BICAN-4x1-HIPU-E-01_VMSC31910")
- **reg_name** (str): Region name (e.g., "region_UCI-5224")
- **prefix1** (str): Prefix for first shapes layer (e.g., "cellpose_nuc")
- **prefix2** (str): Prefix for second shapes layer (e.g., "cellpose_cell")

### Optional
- **output_dir** (str | Path): Custom output directory. If None, uses {SEGMENTATION_OUT_PATH}/{exp_name}/align/
- **segmentation_store** (str | Path): Custom segmentation store path. Defaults to SEGMENTATION_OUT_PATH env var
- **geometry_mode** (str, default="larger"): Which geometry to keep:
  - "larger": Keep the larger of the two overlapping geometries
  - "prefix1": Keep geometry from prefix1
  - "prefix2": Keep geometry from prefix2
  - "intersection": Keep the intersection polygon
- **min_intersection_area** (float, default=0.0): Minimum intersection area to consider a match

## Return Value

Returns a dictionary containing:
- **joint**: GeoDataFrame with all matched cell pairs and their intersection areas
- **aligned_gdf**: GeoDataFrame with the final aligned geometries
- **output_path**: Path where the parquet file was saved
- **n_aligned**: Number of aligned cell pairs
- **geometry_mode**: The geometry mode used

## Geometry Mode Examples

### 1. "larger" (default)
Keeps whichever overlapping cell has the larger area. Good for when one segmentation is more reliable than the other.

### 2. "prefix1" 
Always keeps the geometry from prefix1 (e.g., nuclear). Useful when prefix1 is your reference.

### 3. "prefix2"
Always keeps the geometry from prefix2 (e.g., cellular). Useful when prefix2 is your reference.

### 4. "intersection"
Keeps only the overlapping region. This creates the smallest possible geometries - useful for strict cell definitions.

## Integration with VPT

The aligned parquet file can be directly used with VPT to create new cell_by_gene matrices:

```python
from spida.S.segmentation import seg_to_vpt

# After alignment, use the output with VPT
seg_to_vpt(
    root_dir="/path/to/root",
    seg_out_dir="/path/to/segmentation/output/202511071313_BICAN-4x1-HIPU-E-01_VMSC31910/align",
    region="region_UCI-5224",
    vpt_filepaths={
        "input_boundaries": "aligned_cellpose_nuc_cellpose_cell_202511071313_BICAN-4x1-HIPU-E-01_VMSC31910_region_UCI-5224_polygons.parquet"
    }
)
```

## Output File Location

The aligned geometries are saved as:
```
{output_dir}/aligned_{prefix1}_{prefix2}_{exp_name}_{reg_name}_polygons.parquet
```

Example:
```
/segmentation/output/202511071313_BICAN-4x1-HIPU-E-01_VMSC31910/align/aligned_cellpose_nuc_cellpose_cell_202511071313_BICAN-4x1-HIPU-E-01_VMSC31910_region_UCI-5224_polygons.parquet
```

## GeoDataFrame Columns

The output GeoDataFrame contains:
- **geometry**: The aligned polygon geometry (selected based on geometry_mode)
- **{prefix1}_id**: Original cell ID from prefix1 segmentation
- **{prefix2}_id**: Original cell ID from prefix2 segmentation
- Any other columns from the original shape layers

## Filtering Weak Overlaps

Use `min_intersection_area` to filter out weak overlaps:

```python
result = align_segmentations(
    zarr_path="/path/to/zarr",
    exp_name="experiment",
    reg_name="region",
    prefix1="cellpose_nuc",
    prefix2="cellpose_cell",
    min_intersection_area=100.0  # Only keep matches with intersection >= 100 pixels
)
```

## Logging

The function logs detailed information about the alignment process:
- Number of initial intersections found
- Number of one-to-one mappings after deduplication
- Final number of aligned cells
- Output file location

Check the logs to monitor progress and verify the alignment quality.

## Example Workflow

```python
from spida.S.segmentation import align_segmentations, seg_to_vpt

# Step 1: Align segmentations
result = align_segmentations(
    zarr_path="/home/x-aklein2/projects/aklein/BICAN/HIPP/data/zarr_store/202511071313_BICAN-4x1-HIPU-E-01_VMSC31910/region_UCI-5224",
    exp_name="202511071313_BICAN-4x1-HIPU-E-01_VMSC31910",
    reg_name="region_UCI-5224",
    prefix1="cellpose_nuc",
    prefix2="cellpose_cell",
    geometry_mode="larger"
)

print(f"✓ Aligned {result['n_aligned']} cells")
print(f"✓ Saved to {result['output_path']}")

# Step 2: Generate new cell_by_gene matrix using VPT
seg_to_vpt(
    root_dir="/path/to/root",
    seg_out_dir=str(result['output_path']).rsplit('/', 1)[0],
    region="region_UCI-5224",
    vpt_filepaths={
        "input_boundaries": Path(result['output_path']).name
    }
)

print("✓ New cell_by_gene matrix generated!")
```
