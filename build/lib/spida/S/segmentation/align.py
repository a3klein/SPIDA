import os
from pathlib import Path
import logging
from dotenv import load_dotenv  # type: ignore

import spatialdata as sd
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

load_dotenv()
logger = logging.getLogger(__package__)


def align_segmentations(
    zarr_path: str | Path,
    exp_name: str,
    reg_name: str,
    prefix1: str,
    prefix2: str,
    output_dir: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    geometry_mode: str = "larger",
    cell_id: str = "EntityID",
    coordinate_system: str = "global",
    min_intersection_area: float = 0.0,
    
):
    """
    Align two sets of segmentations (e.g., nuclear and cellular) by their spatial overlap,
    keeping only cells present in both datasets with one-to-one mappings.

    This function identifies matching cells between two shape layers based on their intersection area,
    keeps only one-to-one mappings, and saves the aligned geometries as a parquet file for downstream
    processing with VPT functions to create new cell_by_gene matrices.

    The aligned parquet file can be directly used with the VPT segmentation pipeline to generate
    new cell_by_gene matrices based on the aligned geometries.

    Parameters:
    -----------
    zarr_path: str | Path
        Path to the spatialdata zarr store containing the segmentations
    exp_name : str
        Experiment name used to construct keys (e.g., '202511071313_BICAN-4x1-HIPU-E-01_VMSC31910')
    reg_name : str
        Region name used to construct keys (e.g., 'region_UCI-5224')
    prefix1 : str
        Prefix for the first shapes layer (e.g., 'cellpose_nuc')
    prefix2 : str
        Prefix for the second shapes layer (e.g., 'cellpose_cell')
    output_dir : str | Path, optional
        Directory to save the output parquet file. If None, uses segmentation_store environment variable.
        Default structure: {segmentation_store}/{exp_name}/align/{reg_name}/
    segmentation_store : str | Path, optional
        Default path to store the segmentation results. If None, uses SEGMENTATION_OUT_PATH environment variable.
    geometry_mode : str, default='larger'
        How to select geometry for aligned cells:
        - 'larger': Keep the larger of the two overlapping geometries
        - 'prefix1': Keep geometry from prefix1
        - 'prefix2': Keep geometry from prefix2
        - 'intersection': Keep the intersection polygon of the two geometries
    min_intersection_area : float, default=0.0
        Minimum intersection area (in pixels) required to consider a match. Set > 0 to filter
        out weak overlaps.
    cell_id : str, default='EntityID'
        Column name in the shape layers that contains the unique cell identifier.
    coordinate_system : str, default='global'
        Coordinate system to which shapes should be transformed before alignment.


    Returns:
    --------
    dict
        Dictionary containing:
        - 'joint': GeoDataFrame with matched cell pairs and their intersection areas
        - 'aligned_gdf': GeoDataFrame with the final aligned geometries (ready for VPT)
        - 'output_path': Path where the parquet file was saved
        - 'n_aligned': Number of aligned cell pairs
        - 'geometry_mode': The geometry mode used for alignment

    Raises:
    -------
    FileNotFoundError: If zarr_path or shapes don't exist
    ValueError: If geometry_mode is not valid or if no intersecting cells are found

    Example:
    --------
    >>> from spida.S.segmentation import align_segmentations
    >>> result = align_segmentations(
    ...     zarr_path="/path/to/zarr/store",
    ...     exp_name="202511071313_BICAN-4x1-HIPU-E-01_VMSC31910",
    ...     reg_name="region_UCI-5224",
    ...     prefix1="cellpose_nuc",
    ...     prefix2="cellpose_cell",
    ...     geometry_mode="larger"
    ... )
    >>> print(f"Aligned {result['n_aligned']} cell pairs")
    >>> print(f"Output saved to {result['output_path']}")
    """

    # Validate geometry_mode
    valid_modes = {"larger", "prefix1", "prefix2", "intersection"}
    if geometry_mode not in valid_modes:
        raise ValueError(
            f"Invalid geometry_mode '{geometry_mode}'. Must be one of {valid_modes}"
        )

    # Setup output directory
    if output_dir is None:
        if segmentation_store is None:
            seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        else:
            seg_out_path = segmentation_store
        output_dir = f"{seg_out_path}/{reg_name}"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading spatialdata from {zarr_path}")
    sdata = sd.read_zarr(zarr_path)

    # Generate shape keys
    shape_key1 = f"{prefix1}_{exp_name}_{reg_name}_polygons"
    shape_key2 = f"{prefix2}_{exp_name}_{reg_name}_polygons"

    logger.info(f"Loading shapes from {shape_key1} and {shape_key2}")

    if shape_key1 not in sdata.shapes:
        raise FileNotFoundError(f"Shape layer '{shape_key1}' not found in spatialdata")
    if shape_key2 not in sdata.shapes:
        raise FileNotFoundError(f"Shape layer '{shape_key2}' not found in spatialdata")

    shapes1 = sd.transform(sdata[shape_key1], to_coordinate_system=coordinate_system)
    shapes2 = sd.transform(sdata[shape_key2], to_coordinate_system=coordinate_system)
    # shapes1 = sdata[shape_key1]
    # shapes2 = sdata[shape_key2]

    # Set index names for clarity
    shapes1.index.name = f"{prefix1}_index"
    shapes2.index.name = f"{prefix2}_index"

    logger.info(
        f"Shape counts - {prefix1}: {len(shapes1)}, {prefix2}: {len(shapes2)}"
    )

    # Perform spatial join to find intersecting shapes
    logger.info("Finding intersecting shapes...")
    joint = gpd.sjoin(
        shapes1, shapes2, how="inner", predicate="intersects",
        lsuffix=prefix1, rsuffix=prefix2
    )

    if len(joint) == 0:
        raise ValueError(
            f"No intersecting cells found between {prefix1} and {prefix2}"
        )

    logger.info(f"Found {len(joint)} initial intersections")
    joint = joint.reset_index()

    # Calculate intersection areas
    logger.info("Calculating intersection areas...")
    joint['intersection_area'] = joint.apply(
        lambda row: shapes1.loc[str(row[f'{prefix1}_index'])].geometry.intersection(
            shapes2.loc[str(row[f'{prefix2}_index'])].geometry
        ).area,
        axis=1
    )

    # Filter by minimum intersection area
    if min_intersection_area > 0:
        joint = joint.query('intersection_area >= @min_intersection_area')
        logger.info(f"After filtering by min area: {len(joint)} intersections")

    # Keep only the best match for each cell in both directions (one-to-one mapping)
    logger.info("Enforcing one-to-one cell mappings...")

    # First, for each cell in shapes1, keep the shapes2 cell with largest intersection
    joint = joint.sort_values(
        by=[f"{prefix1}_index", 'intersection_area'],
        ascending=[True, False]
    )
    joint = joint.drop_duplicates(subset=[f'{prefix1}_index'], keep='first')

    # Then, for each cell in shapes2, keep the shapes1 cell with largest intersection
    joint = joint.sort_values(
        by=[f"{prefix2}_index", 'intersection_area'],
        ascending=[True, False]
    )
    joint = joint.drop_duplicates(subset=[f'{prefix2}_index'], keep='first')

    logger.info(f"After one-to-one mapping: {len(joint)} aligned cell pairs")

    # Extract the IDs of aligned cells
    aligned_ids1 = joint[f'{prefix1}_index'].astype(str).tolist()
    aligned_ids2 = joint[f'{prefix2}_index'].astype(str).tolist()

    # Create aligned shapes based on geometry_mode
    logger.info(f"Creating aligned geometries using mode: {geometry_mode}")
    aligned_shapes_list = []
    aligned_ids_list = []

    for idx_row in joint.itertuples():
        id1 = str(idx_row[idx_row._fields.index(f'{prefix1}_index')])
        id2 = str(idx_row[idx_row._fields.index(f'{prefix2}_index')])
        shape1 = shapes1.loc[[id1]].iloc[0]
        shape2 = shapes2.loc[[id2]].iloc[0]

        if geometry_mode == "larger":
            # Keep the larger shape
            if shape1.geometry.area >= shape2.geometry.area:
                aligned_shapes_list.append(shape1.copy())
            else:
                aligned_shapes_list.append(shape2.copy())
        elif geometry_mode == "prefix1":
            # Keep geometry from prefix1
            aligned_shapes_list.append(shape1.copy())
        elif geometry_mode == "prefix2":
            # Keep geometry from prefix2
            aligned_shapes_list.append(shape2.copy())
        elif geometry_mode == "intersection":
            # Keep intersection polygon
            intersection_geom = shape1.geometry.intersection(shape2.geometry)
            aligned_shape = shape1.copy()
            aligned_shape['geometry'] = intersection_geom
            aligned_shapes_list.append(aligned_shape)

        aligned_ids_list.append(id1)

    # Create GeoDataFrame with aligned geometries
    aligned_gdf = gpd.GeoDataFrame(
        aligned_shapes_list, crs=shapes1.crs, geometry="geometry"
    )

    # Reset index to use sequential IDs
    aligned_gdf = aligned_gdf.reset_index(drop=True)

    # Add aligned IDs as a column for reference
    aligned_gdf[f'{prefix1}_id'] = aligned_ids_list
    aligned_gdf[f'{prefix2}_id'] = joint[f'{prefix2}_index'].astype(str).values

    logger.info(
        f"Created {len(aligned_gdf)} aligned cell geometries using mode '{geometry_mode}'"
    )

    # For compatibility with VPT, rename geometry column to "Geometry"
    aligned_gdf.rename_geometry("Geometry", inplace=True)
    logger.info("Renamed geometry column to 'Geometry' for VPT compatibility")

    # Save to parquet
    output_filename = f"aligned_{prefix1}_{prefix2}_polygons.parquet"
    output_path = output_dir / output_filename

    logger.info(f"Saving aligned geometries to {output_path}")
    aligned_gdf.to_parquet(str(output_path), index=False)

    logger.info(f"Successfully saved {len(aligned_gdf)} aligned cells to {output_path}")

    return {
        'joint': joint,
        'aligned_gdf': aligned_gdf,
        'output_path': str(output_path),
        'n_aligned': len(aligned_gdf),
        'geometry_mode': geometry_mode,
    }
