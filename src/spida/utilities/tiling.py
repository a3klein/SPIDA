from pathlib import Path
from collections.abc import Callable

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely as shp
import skimage as ski
import imageio.v3 as iio

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm


def tile_image_with_overlap(
    image: np.ndarray,
    tile_size: int | tuple,
    overlap: int | tuple = 0,
    color: str = None,
):
    """
    Tile a large image into smaller tiles with specified overlap.

    Parameters:
    -----------
    image : numpy.ndarray
        Input image array (2D or 3D)
    tile_size : int or tuple
        Size of each tile. If int, creates square tiles.
        If tuple, should be (height, width) or (depth, height, width) for 3D
    overlap : int or tuple, optional
        Overlap between tiles in pixels. If int, same overlap for all dimensions.
        If tuple, should match tile_size dimensions. Default is 0.
    color : str, optional
        Color of channel to append to tile info

    Returns:
    --------
    tiles : list of numpy.ndarray
        List of tile arrays
    tile_info : list of dict
        List of dictionaries containing tile information:
        - 'position': starting position (row, col) or (depth, row, col)
        - 'end_position': ending position
        - 'tile_id': unique identifier for the tile
        - 'overlap_region': dict with overlap information

    Examples:
    ---------
    # For 2D image
    tiles, info = tile_image_with_overlap(img, tile_size=512, overlap=64)

    # For 3D image
    tiles, info = tile_image_with_overlap(img_3d, tile_size=(7, 512, 512), overlap=(0, 64, 64))
    """

    # Handle input validation
    if not isinstance(image, np.ndarray):
        raise ValueError("Image must be a numpy array")

    # Convert tile_size to tuple if it's an integer
    if isinstance(tile_size, int):
        if image.ndim == 2:
            tile_size = (tile_size, tile_size)
        elif image.ndim == 3:
            tile_size = (image.shape[0], tile_size, tile_size)
        else:
            raise ValueError("Unsupported image dimensions")

    # Convert overlap to tuple if it's an integer
    if isinstance(overlap, int):
        overlap = tuple([overlap] * len(tile_size))

    # Validate dimensions
    if len(tile_size) != image.ndim:
        raise ValueError(
            f"tile_size dimensions ({len(tile_size)}) must match image dimensions ({image.ndim})"
        )

    if len(overlap) != image.ndim:
        raise ValueError(
            f"overlap dimensions ({len(overlap)}) must match image dimensions ({image.ndim})"
        )

    tiles = []
    tile_info = []
    tile_id = 0

    # Calculate step size (tile_size - overlap)
    step_size = tuple(ts - ov for ts, ov in zip(tile_size, overlap))

    # Generate tile positions
    if image.ndim == 2:
        height, width = image.shape
        tile_h, tile_w = tile_size
        step_h, step_w = step_size

        for row in range(0, height - tile_h + 1, step_h):
            for col in range(0, width - tile_w + 1, step_w):
                # Handle edge cases where we might go beyond image boundaries
                end_row = min(row + tile_h, height)
                end_col = min(col + tile_w, width)

                # Extract tile
                tile = image[row:end_row, col:end_col]

                # Store tile information
                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (end_row, end_col),
                    "actual_size": tile.shape,
                    "overlap_region": {
                        "top": row > 0,
                        "bottom": end_row < height,
                        "left": col > 0,
                        "right": end_col < width,
                    },
                    "color": color,
                }

                tiles.append(tile)
                tile_info.append(info)
                tile_id += 1

        # Handle remaining edge tiles
        # Right edge
        if width % step_w != 0:
            col = width - tile_w
            for row in range(0, height - tile_h + 1, step_h):
                end_row = min(row + tile_h, height)
                tile = image[row:end_row, col:width]

                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (end_row, width),
                    "actual_size": tile.shape,
                    "overlap_region": {
                        "top": row > 0,
                        "bottom": end_row < height,
                        "left": True,
                        "right": False,
                    },
                    "color": color,
                }

                tiles.append(tile)
                tile_info.append(info)
                tile_id += 1

        # Bottom edge
        if height % step_h != 0:
            row = height - tile_h
            for col in range(0, width - tile_w + 1, step_w):
                end_col = min(col + tile_w, width)
                tile = image[row:height, col:end_col]

                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (height, end_col),
                    "actual_size": tile.shape,
                    "overlap_region": {
                        "top": True,
                        "bottom": False,
                        "left": col > 0,
                        "right": end_col < width,
                    },
                    "color": color,
                }

                tiles.append(tile)
                tile_info.append(info)
                tile_id += 1

        # Bottom-right corner
        if height % step_h != 0 and width % step_w != 0:
            row = height - tile_h
            col = width - tile_w
            tile = image[row:height, col:width]

            info = {
                "tile_id": tile_id,
                "position": (row, col),
                "end_position": (height, width),
                "actual_size": tile.shape,
                "overlap_region": {
                    "top": True,
                    "bottom": False,
                    "left": True,
                    "right": False,
                },
                "color": color,
            }

            tiles.append(tile)
            tile_info.append(info)
            tile_id += 1

    # TODO: Implement for 3D images if needed
    else:
        raise ValueError("Only 2D images are supported")

    return tiles, tile_info


def colormap_mapper(values, colormap="viridis"):
    """
    Map values to a colormap.

    Parameters:
    - values: array-like, values to map
    - colormap: str, name of the colormap to use

    Returns:
    - mapped_colors: array-like, colors corresponding to the input values
    """
    norm = Normalize(vmin=min(values), vmax=max(values))
    cmap = cm.get_cmap(colormap)
    return lambda x: cmap(norm(x))


def visualize_tiling_grid(tile_info: list, image_shape: tuple, color_prop: str = None):
    """
    Visualize the tiling grid on the image to preview how tiles will be arranged.

    Parameters:
    -----------
    image_shape : tuple
        Shape of the original image
    tile_size : int or tuple
        Size of each tile
    overlap : int or tuple, optional
        Overlap between tiles
    max_tiles_to_show : int, optional
        Maximum number of tiles to show in visualization

    Returns:
    --------
    fig : matplotlib figure
        Figure showing the tiling grid
    """
    if len(image_shape) == 2:
        height, width = image_shape

        fig, ax = plt.subplots(figsize=(12, 8))

        # Draw image boundary
        ax.add_patch(
            plt.Rectangle(
                (0, 0), width, height, fill=False, edgecolor="black", linewidth=2
            )
        )

        if color_prop is not None:
            # Create a colormap based on the color property
            values = [tile[color_prop] for tile in tile_info]
            mapper = colormap_mapper(values, colormap="Grays")

        # Draw tile boundaries
        for tile in tile_info:
            id = tile["tile_id"]
            row, col = tile["position"]
            end_row, end_col = tile["end_position"]
            tile_h, tile_w = tile["actual_size"]
            cprop = tile.get(color_prop, None)

            if color_prop is not None:
                # Draw tile rectangle
                rect = plt.Rectangle(
                    (col, height - end_row),
                    end_col - col,
                    end_row - row,
                    facecolor=mapper(cprop),
                    edgecolor="red",
                    linewidth=1,
                    alpha=1,
                )
                ax.add_patch(rect)
            else:
                # Draw tile rectangle
                rect = plt.Rectangle(
                    (col, height - end_row),
                    end_col - col,
                    end_row - row,
                    fill=False,
                    edgecolor="red",
                    linewidth=1,
                    alpha=0.7,
                )
                ax.add_patch(rect)

                # Add tile number
                ax.text(
                    col + (end_col - col) / 2,
                    height - (row + (end_row - row) / 2),
                    str(id),
                    ha="center",
                    va="center",
                    fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                )

        overlap = tile_info[0]["end_position"][1] - tile_info[1]["position"][1]

        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect("equal")
        ax.set_title(
            f"Tiling Grid Preview\nImage: {width}x{height}, Tile: {tile_w}x{tile_h}, Overlap: {overlap}"
        )
        ax.set_xlabel("Width")
        ax.set_ylabel("Height")

        plt.tight_layout()
        return fig

    else:
        print("Visualization only supports 2D images")
        return None


def save_tiles(
    tiles,
    tile_info,
    output_dir,
    prefix="tile",
    format="tif",
    func: Callable[[np.ndarray], np.ndarray] = None,
):
    """
    Save tiles to individual files.

    Parameters:
    -----------
    tiles : list of numpy.ndarray
        List of tile arrays
    tile_info : list of dict
        List of tile information dictionaries
    output_dir : str or Path
        Directory to save tiles
    prefix : str, optional
        Prefix for tile filenames
    format : str, optional
        Image format ('tif', 'png', 'jpg')
    func : callable, optional
        Function to apply to each tile before saving (e.g., z_stack generation)

    Returns:
    --------
    saved_files : list
        List of saved file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    saved_files = []

    for tile, info in zip(tiles, tile_info):
        # Create filename
        tile_id = info["tile_id"]
        pos = info["position"]
        col = info["color"]

        if len(pos) == 2:  # 2D
            filename = f"{prefix}_{tile_id:04d}_r{pos[0]}_c{pos[1]}_col{col}.{format}"
        else:  # 3D
            filename = f"{prefix}_{tile_id:04d}_d{pos[0]}_r{pos[1]}_c{pos[2]}_col{col}.{format}"

        filepath = output_dir / filename

        # applying func to the tile
        if func is not None:
            tile = func(tile)

        # Save tile
        if format.lower() == "tif":
            iio.imwrite(filepath, tile)
        else:
            iio.imwrite(filepath, tile)

        saved_files.append(filepath)

    return saved_files


def reconstruct_image_from_tile_files(
    output_dir,
    tile_info,
    original_shape,
    overlap_strategy="average",
    prefix="tile",
    format="tif",
    suffix=None,
    match_pre: bool = False,
):
    """
    Reconstruct the original image from tiles stored in a directory.

    Parameters:
    -----------
    output_dir : str or Path
        Directory containing the saved tile files
    tile_info : list of dict
        List of tile information dictionaries
    original_shape : tuple
        Shape of the original image
    overlap_strategy : str, optional
        How to handle overlapping regions ('average', 'max')
    prefix : str, optional
        Prefix used for tile filenames (default: "tile")
    format : str, optional
        Image format ('tif', 'png', 'jpg') (default: "tif")
    suffix : str, optional
        Additional suffix for tile filenames (e.g., ".decon" for deconvolved tiles)

    Returns:
    --------
    reconstructed : numpy.ndarray
        Reconstructed image
    """
    output_dir = Path(output_dir)

    # Initialize output array and weight map
    reconstructed = None
    weight_map = None

    # Load and place tiles back into the reconstructed image one at a time
    for info in tile_info:
        # Build filename based on tile info
        tile_id = info["tile_id"]
        pos = info["position"]
        col = info["color"]

        if len(pos) == 2:  # 2D
            filename = f"{prefix}_{tile_id:04d}_r{pos[0]}_c{pos[1]}_col{col}"
        else:  # 3D
            filename = f"{prefix}_{tile_id:04d}_d{pos[0]}_r{pos[1]}_c{pos[2]}_col{col}"

        # Add suffix if provided
        input_filename = filename
        if suffix:
            filename += suffix

        input_filename += f".{format}"
        input_filepath = output_dir / input_filename
        filename += f".{format}"
        filepath = output_dir / filename

        # Skip if file doesn't exist
        if not filepath.exists():
            print(f"Warning: Tile file {filepath} not found, skipping...")
            continue

        # Load tile
        tile = iio.imread(filepath)

        # If adjusting the tile color histogram
        if match_pre:
            input_tile = iio.imread(input_filepath)
            input_tile = input_tile.mean(axis=0).astype("int16")
            input_tile[input_tile < 0] = 0  # Ensure no negative values
            tile = ski.exposure.match_histograms(tile, input_tile, channel_axis=None)
            del input_tile  # clear memory

        # Initialize arrays on first tile to get dtype
        if reconstructed is None:
            reconstructed = np.zeros(original_shape, dtype=tile.dtype)
            weight_map = np.zeros(original_shape, dtype=np.float32)

        # Place tile in reconstructed image
        pos = info["position"]
        actual_size = info["actual_size"]

        if len(pos) == 2:  # 2D
            r, c = pos
            reconstructed[r : r + actual_size[0], c : c + actual_size[1]] += tile
            weight_map[r : r + actual_size[0], c : c + actual_size[1]] += 1
        else:  # 3D
            d, r, c = pos
            reconstructed[
                d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
            ] += tile
            weight_map[
                d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
            ] += 1

    # Handle case where no tiles were found
    if reconstructed is None:
        raise ValueError("No valid tile files found in the output directory")

    # Handle overlapping regions
    if overlap_strategy == "average":
        # Avoid division by zero
        weight_map[weight_map == 0] = 1
        reconstructed = reconstructed / weight_map
    elif overlap_strategy == "max":
        # TODO: For max strategy, we need to reconstruct differently
        pass  # Implementation would be more complex

    reconstructed[reconstructed < 0] = 0  # Making sure no negative values!
    return reconstructed.astype(reconstructed.dtype)


def reconstruct_image_from_tiles(
    tiles, tile_info, original_shape, overlap_strategy="average"
):
    """
    Reconstruct the original image from tiles.

    Parameters:
    -----------
    tiles : list of numpy.ndarray
        List of tile arrays
    tile_info : list of dict
        List of tile information dictionaries
    original_shape : tuple
        Shape of the original image
    overlap_strategy : str, optional
        How to handle overlapping regions ('average', 'max')

    Returns:
    --------
    reconstructed : numpy.ndarray
        Reconstructed image
    """
    # Initialize output array
    reconstructed = np.zeros(original_shape, dtype=tiles[0].dtype)
    weight_map = np.zeros(original_shape, dtype=np.float32)

    # Place tiles back into the reconstructed image
    for tile, info in zip(tiles, tile_info):
        pos = info["position"]
        actual_size = info["actual_size"]

        if len(pos) == 2:  # 2D
            r, c = pos
            reconstructed[r : r + actual_size[0], c : c + actual_size[1]] += tile
            weight_map[r : r + actual_size[0], c : c + actual_size[1]] += 1
        else:  # 3D
            d, r, c = pos
            reconstructed[
                d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
            ] += tile
            weight_map[
                d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
            ] += 1

    # Handle overlapping regions
    if overlap_strategy == "average":
        # Avoid division by zero
        weight_map[weight_map == 0] = 1
        reconstructed = reconstructed / weight_map
    elif overlap_strategy == "max":
        # TODO: For max strategy, we need to reconstruct differently
        pass  # Implementation would be more complex

    return reconstructed.astype(tiles[0].dtype)


def save_projected_tiles_from_files(
    input_dir,
    tile_info,
    output_dir=None,
    prefix="tile",
    input_format="tif",
    output_format="tif",
    input_suffix=".decon",
    output_suffix=".decon.2d",
    projection_method="max",
    z_slice_range=None,
):
    """
    Load 3D tiles from files, project them to 2D, and save the projected tiles.

    Parameters:
    -----------
    input_dir : str or Path
        Directory containing the 3D tile files
    tile_info : list of dict
        List of tile information dictionaries
    output_dir : str or Path, optional
        Directory to save projected tiles (default: same as input_dir)
    prefix : str, optional
        Prefix for tile filenames (default: "tile")
    input_format : str, optional
        Input image format (default: "tif")
    output_format : str, optional
        Output image format (default: "tif")
    input_suffix : str, optional
        Additional suffix for input tile filenames (default: ".decon")
    output_suffix : str, optional
        Additional suffix for output tile filenames (default: ".decon.2d")
    projection_method : str, optional
        Projection method ('max', 'mean', 'sum') (default: "max")
    z_slice_range : tuple, optional
        Range of z-slices to use for projection (start, end). If None, uses all slices.
        If (1, -1), removes first and last slices.

    Returns:
    --------
    saved_files : list
        List of saved projected file paths
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir) if output_dir else input_dir
    output_dir.mkdir(exist_ok=True, parents=True)

    saved_files = []

    for info in tile_info:
        # Build input filename
        tile_id = info["tile_id"]
        pos = info["position"]
        col = info["color"]

        if len(pos) == 2:  # 2D
            base_filename = f"{prefix}_{tile_id:04d}_r{pos[0]}_c{pos[1]}_col{col}"
        else:  # 3D
            base_filename = (
                f"{prefix}_{tile_id:04d}_d{pos[0]}_r{pos[1]}_c{pos[2]}_col{col}"
            )

        input_filename = base_filename + input_suffix + f".{input_format}"
        output_filename = base_filename + output_suffix + f".{output_format}"

        input_filepath = input_dir / input_filename
        output_filepath = output_dir / output_filename

        # Skip if output already exists
        if output_filepath.exists():
            saved_files.append(output_filepath)
            continue

        # Skip if input doesn't exist
        if not input_filepath.exists():
            print(f"Warning: Input file {input_filepath} not found, skipping...")
            continue

        # Load 3D tile
        tile_3d = ski.io.imread(input_filepath)

        # Apply z-slice range if specified
        if z_slice_range is not None:
            start, end = z_slice_range
            if end == -1:
                tile_3d = tile_3d[start:]
            else:
                tile_3d = tile_3d[start:end]

        # Project to 2D
        if tile_3d.ndim == 3 and tile_3d.shape[0] > 1:
            if projection_method == "max":
                tile_2d = tile_3d.max(axis=0)
            elif projection_method == "mean":
                tile_2d = tile_3d.mean(axis=0)
            elif projection_method == "sum":
                tile_2d = tile_3d.sum(axis=0)
            else:
                raise ValueError(f"Unknown projection method: {projection_method}")
        else:
            tile_2d = tile_3d.squeeze()

        # Save projected tile
        iio.imwrite(output_filepath, tile_2d)
        saved_files.append(output_filepath)

        # Clear memory
        del tile_3d, tile_2d

    return saved_files


def merge_tile_polygons(
    tile_info_list, polygon_gdfs, coordinate_offset=(0, 0), geom_col="geometry"
):
    """
    Merge polygons from different tiles into a single GeoDataFrame.

    This function takes polygon GeoDataFrames from individual tiles and merges them
    into a single coordinate system by transforming the polygon coordinates based
    on each tile's position information.

    Parameters:
    -----------
    tile_info_list : list of dict
        List of tile information dictionaries from tile_image_with_overlap().
        Each dict should contain 'position' key with (row, col) or (z, row, col) coordinates.
    polygon_gdfs : list of geopandas.GeoDataFrame
        List of GeoDataFrames containing polygons for each tile.
        Should be in the same order as tile_info_list.
    coordinate_offset : tuple, optional
        Global coordinate offset to apply to all polygons. Default is (0, 0).
        Should be (row_offset, col_offset) for 2D or (z_offset, row_offset, col_offset) for 3D.

    Returns:
    --------
    merged_gdf : geopandas.GeoDataFrame
        Single GeoDataFrame containing all polygons with transformed coordinates.
        Includes additional columns:
        - 'tile_id': Original tile identifier
        - 'original_tile_position': Original position within the tile

    Examples:
    ---------
    # Merge polygons from multiple tiles
    merged_polygons = merge_tile_polygons(tile_info, polygon_list)

    # With coordinate offset
    merged_polygons = merge_tile_polygons(tile_info, polygon_list,
                                        coordinate_offset=(1000, 1000))
    """

    if len(tile_info_list) != len(polygon_gdfs):
        raise ValueError("tile_info_list and polygon_gdfs must have the same length")

    merged_polygons = []

    for tile_info, gdf in zip(tile_info_list, polygon_gdfs):
        if gdf.empty:
            continue

        # Get tile position - handle both 2D and 3D cases
        tile_position = tile_info["position"]
        if len(tile_position) == 2:
            tile_row, tile_col = tile_position
            tile_z = 0
        elif len(tile_position) == 3:
            tile_z, tile_row, tile_col = tile_position
        else:
            raise ValueError(
                f"Unsupported tile position dimensions: {len(tile_position)}"
            )

        # Apply coordinate offset
        if len(coordinate_offset) == 2:
            offset_row, offset_col = coordinate_offset
            offset_z = 0
        elif len(coordinate_offset) == 3:
            offset_z, offset_row, offset_col = coordinate_offset
        else:
            raise ValueError(
                f"Unsupported coordinate_offset dimensions: {len(coordinate_offset)}"
            )

        # Calculate total offset
        total_row_offset = tile_row + offset_row
        total_col_offset = tile_col + offset_col
        total_z_offset = tile_z + offset_z

        # Transform each polygon in the GeoDataFrame
        transformed_polygons = []
        for idx, row in gdf.iterrows():
            polygon = row[geom_col]

            if polygon is None or polygon.is_empty:
                continue

            # Transform polygon coordinates
            if hasattr(polygon, "exterior"):  # Single polygon
                # Get exterior coordinates
                exterior_coords = list(polygon.exterior.coords)
                # Transform coordinates: (x, y) -> (x + col_offset, y + row_offset)
                transformed_coords = [
                    (x + total_col_offset, y + total_row_offset)
                    for x, y in exterior_coords
                ]
                transformed_polygon = shp.geometry.Polygon(transformed_coords)

                # Handle holes if they exist
                if polygon.interiors:
                    holes = []
                    for interior in polygon.interiors:
                        hole_coords = list(interior.coords)
                        transformed_hole = [
                            (x + total_col_offset, y + total_row_offset)
                            for x, y in hole_coords
                        ]
                        holes.append(transformed_hole)
                    transformed_polygon = shp.geometry.Polygon(
                        transformed_coords, holes
                    )

            else:
                # Handle other geometry types (MultiPolygon, etc.)
                transformed_polygon = shp.affinity.translate(
                    polygon, xoff=total_col_offset, yoff=total_row_offset
                )

            # Create new row with transformed geometry
            new_row = row.copy()
            new_row[geom_col] = transformed_polygon
            new_row["tile_id"] = tile_info["tile_id"]
            new_row["original_tile_position"] = tile_position
            new_row["tile_row_offset"] = total_row_offset
            new_row["tile_col_offset"] = total_col_offset

            # Add z information if relevant
            if total_z_offset != 0:
                new_row["tile_z_offset"] = total_z_offset
                new_row["z"] = getattr(row, "z", 0) + total_z_offset

            transformed_polygons.append(new_row)

        # Add transformed polygons to the list
        if transformed_polygons:
            tile_gdf = gpd.GeoDataFrame(transformed_polygons)
            merged_polygons.append(tile_gdf)

    # Concatenate all GeoDataFrames
    if merged_polygons:
        merged_gdf = gpd.GeoDataFrame(pd.concat(merged_polygons, ignore_index=True))

        # Ensure CRS is set (use the CRS from the first non-empty GDF)
        for gdf in polygon_gdfs:
            if not gdf.empty and gdf.crs is not None:
                merged_gdf.crs = gdf.crs
                break

        return merged_gdf.set_geometry(geom_col)
    else:
        # Return empty GeoDataFrame with expected structure
        return gpd.GeoDataFrame(
            columns=[
                geom_col,
                "tile_id",
                "original_tile_position",
                "tile_row_offset",
                "tile_col_offset",
            ]
        ).set_geometry(geom_col)


def merge_overlapping_polygons(merged_gdf, buffer_distance=0.1, area_threshold=None):
    """
    Merge overlapping polygons that might result from tile boundaries.

    This function identifies and merges polygons that overlap significantly,
    which commonly occurs when the same object spans multiple tiles.

    Parameters:
    -----------
    merged_gdf : geopandas.GeoDataFrame
        GeoDataFrame with polygons from merge_tile_polygons()
    buffer_distance : float, optional
        Buffer distance for polygon intersection detection. Default is 0.1.
    area_threshold : float, optional
        Minimum overlap area threshold for merging. If None, uses 50% of smaller polygon.

    Returns:
    --------
    merged_gdf : geopandas.GeoDataFrame
        GeoDataFrame with overlapping polygons merged.
    """

    if merged_gdf.empty:
        return merged_gdf

    # Create spatial index for efficient overlap detection
    spatial_index = merged_gdf.sindex
    merged_gdf = merged_gdf.copy()
    merged_gdf["merged"] = False

    # Track polygons to merge
    to_merge = []
    processed = set()

    for idx, row in merged_gdf.iterrows():
        if idx in processed:
            continue

        polygon = row.geometry
        if polygon is None or polygon.is_empty:
            continue

        # Buffer polygon slightly for intersection detection
        buffered = polygon.buffer(buffer_distance)

        # Find potential intersecting polygons
        possible_matches_idx = list(spatial_index.intersection(buffered.bounds))

        overlapping_indices = [idx]

        for match_idx in possible_matches_idx:
            if match_idx == idx or match_idx in processed:
                continue

            other_polygon = merged_gdf.loc[match_idx, "geometry"]
            if other_polygon is None or other_polygon.is_empty:
                continue

            # Check for significant overlap
            intersection = polygon.intersection(other_polygon)
            if not intersection.is_empty:
                intersection_area = intersection.area

                # Determine threshold
                if area_threshold is None:
                    min_area = min(polygon.area, other_polygon.area)
                    threshold = min_area * 0.5  # 50% overlap
                else:
                    threshold = area_threshold

                if intersection_area >= threshold:
                    overlapping_indices.append(match_idx)

        if len(overlapping_indices) > 1:
            to_merge.append(overlapping_indices)
            processed.update(overlapping_indices)

    # Merge overlapping polygons
    for merge_group in to_merge:
        polygons_to_merge = [merged_gdf.loc[idx, "geometry"] for idx in merge_group]

        # Union all polygons in the group
        merged_polygon = polygons_to_merge[0]
        for poly in polygons_to_merge[1:]:
            merged_polygon = merged_polygon.union(poly)

        # Update the first polygon in the group with merged geometry
        first_idx = merge_group[0]
        merged_gdf.loc[first_idx, "geometry"] = merged_polygon
        merged_gdf.loc[first_idx, "merged"] = True

        # Mark other polygons for removal
        for idx in merge_group[1:]:
            merged_gdf.loc[idx, "merged"] = True
            merged_gdf.loc[idx, "geometry"] = None

    # Remove merged polygons (keep only the first one from each group)
    result_gdf = merged_gdf[merged_gdf.geometry.notna()].copy()

    return result_gdf
