from pathlib import Path
from collections.abc import Callable
from dataclasses import dataclass
import logging

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely as shp
import skimage as ski
import imageio.v3 as iio
import math

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

logger = logging.getLogger(__name__)


def tile_grid_with_overlap(
    shape: tuple,
    tile_size: int | tuple,
    overlap: int | tuple = 0,
    color: str = None,
):
    """
    Compute tile positions covering an array shape, WITHOUT needing the array.

    This is the geometry half of :func:`tile_image_with_overlap`, split out so
    callers that cannot hold the data in memory can still lay out a grid. That is
    the whole reason it exists: a full-resolution MERSCOPE z-stack is ~64 GB, so
    GPU-scale segmentation tiling has to plan from coordinates and then read each
    window lazily.

    Positions, ordering and edge handling are identical to
    :func:`tile_image_with_overlap` — that function now calls this one and simply
    slices the result, so the two can never drift apart.

    Parameters:
    -----------
    shape : tuple
        Shape of the array to be tiled (2D or 3D).
    tile_size : int or tuple
        Size of each tile. If int, square/cubic tiles (for 3D the depth is taken
        as the full extent). If tuple, must match ``len(shape)``.
    overlap : int or tuple, optional
        Overlap in pixels. If int, applied to every axis (for 3D, depth gets 0).
    color : str, optional
        Passed through onto each tile_info dict.

    Returns:
    --------
    tile_info : list of dict
        One dict per tile: ``position``, ``end_position``, ``actual_size``,
        ``tile_id``, ``overlap_region``, ``color``.

    Examples:
    ---------
    info = tile_grid_with_overlap((82428, 55905), tile_size=(12250, 9500), overlap=500)
    """

    # Convert tile_size to tuple if it's an integer
    if isinstance(tile_size, int):
        if len(shape) == 2:
            tile_size = (tile_size, tile_size)
        elif len(shape) == 3:
            tile_size = (shape[0], tile_size, tile_size)
        else:
            raise ValueError("Unsupported image dimensions")

    # Convert overlap to tuple if it's an integer
    if isinstance(overlap, int):
        if len(shape) == 2:
            overlap = tuple([overlap] * len(tile_size))
        elif len(shape) == 3:
            overlap = (0,) + tuple([overlap] * (len(tile_size) - 1))
        else:
            raise ValueError("Unsupported image dimensions")

    # Validate dimensions
    if len(tile_size) != len(shape):
        raise ValueError(
            f"tile_size dimensions ({len(tile_size)}) must match image dimensions ({len(shape)})"
        )

    if len(overlap) != len(shape):
        raise ValueError(
            f"overlap dimensions ({len(overlap)}) must match image dimensions ({len(shape)})"
        )

    # Clamp a tile that is LARGER than the extent, before step_size is derived.
    tile_size = tuple(min(ts, dim) for ts, dim in zip(tile_size, shape))
    # Get overlap while settin goverlap to 0 for tile that spans an axis entirely.
    overlap = tuple(0 if ts >= dim else ov
                    for ts, ov, dim in zip(tile_size, overlap, shape))

    tile_info = []
    tile_id = 0

    # Calculate step size (tile_size - overlap)
    step_size = tuple(ts - ov for ts, ov in zip(tile_size, overlap))

    # Generate tile positions
    if len(shape) == 2:
        height, width = shape
        tile_h, tile_w = tile_size
        step_h, step_w = step_size

        for row in range(0, height - tile_h + 1, step_h):
            for col in range(0, width - tile_w + 1, step_w):
                # Handle edge cases where we might go beyond image boundaries
                end_row = min(row + tile_h, height)
                end_col = min(col + tile_w, width)

                # Extract tile

                # Store tile information
                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (end_row, end_col),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "top": row > 0,
                        "bottom": end_row < height,
                        "left": col > 0,
                        "right": end_col < width,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Handle remaining edge tiles
        # Right edge
        if (width - tile_w) % step_w != 0:
            col = width - tile_w
            for row in range(0, height - tile_h + 1, step_h):
                end_row = min(row + tile_h, height)

                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (end_row, width),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "top": row > 0,
                        "bottom": end_row < height,
                        "left": True,
                        "right": False,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Bottom edge
        if (height - tile_h) % step_h != 0:
            row = height - tile_h
            for col in range(0, width - tile_w + 1, step_w):
                end_col = min(col + tile_w, width)

                info = {
                    "tile_id": tile_id,
                    "position": (row, col),
                    "end_position": (height, end_col),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "top": True,
                        "bottom": False,
                        "left": col > 0,
                        "right": end_col < width,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Bottom-right corner
        if (height - tile_h) % step_h != 0 and (width - tile_w) % step_w != 0:
            row = height - tile_h
            col = width - tile_w

            info = {
                "tile_id": tile_id,
                "position": (row, col),
                "end_position": (height, width),
                "actual_size": None,  # filled below
                "overlap_region": {
                    "top": True,
                    "bottom": False,
                    "left": True,
                    "right": False,
                },
                "color": color,
            }

            tile_info.append(info)
            tile_id += 1

    elif len(shape) == 3:
        depth, height, width = shape
        tile_d, tile_h, tile_w = tile_size
        step_d, step_h, step_w = step_size

        # Main grid of tiles
        for d in range(0, depth - tile_d + 1, step_d):
            for row in range(0, height - tile_h + 1, step_h):
                for col in range(0, width - tile_w + 1, step_w):
                    # Handle edge cases where we might go beyond image boundaries
                    end_d = min(d + tile_d, depth)
                    end_row = min(row + tile_h, height)
                    end_col = min(col + tile_w, width)

                    # Extract tile

                    # Store tile information
                    info = {
                        "tile_id": tile_id,
                        "position": (d, row, col),
                        "end_position": (end_d, end_row, end_col),
                        "actual_size": None,  # filled below
                        "overlap_region": {
                            "front": d > 0,
                            "back": end_d < depth,
                            "top": row > 0,
                            "bottom": end_row < height,
                            "left": col > 0,
                            "right": end_col < width,
                        },
                        "color": color,
                    }

                    tile_info.append(info)
                    tile_id += 1

        # Handle edge tiles for 3D
        # Right edge tiles
        if (width - tile_w) % step_w != 0:
            col = width - tile_w
            for d in range(0, depth - tile_d + 1, step_d):
                for row in range(0, height - tile_h + 1, step_h):
                    end_d = min(d + tile_d, depth)
                    end_row = min(row + tile_h, height)

                    info = {
                        "tile_id": tile_id,
                        "position": (d, row, col),
                        "end_position": (end_d, end_row, width),
                        "actual_size": None,  # filled below
                        "overlap_region": {
                            "front": d > 0,
                            "back": end_d < depth,
                            "top": row > 0,
                            "bottom": end_row < height,
                            "left": True,
                            "right": False,
                        },
                        "color": color,
                    }

                    tile_info.append(info)
                    tile_id += 1

        # Bottom edge tiles
        if (height - tile_h) % step_h != 0:
            row = height - tile_h
            for d in range(0, depth - tile_d + 1, step_d):
                for col in range(0, width - tile_w + 1, step_w):
                    end_d = min(d + tile_d, depth)
                    end_col = min(col + tile_w, width)

                    info = {
                        "tile_id": tile_id,
                        "position": (d, row, col),
                        "end_position": (end_d, height, end_col),
                        "actual_size": None,  # filled below
                        "overlap_region": {
                            "front": d > 0,
                            "back": end_d < depth,
                            "top": True,
                            "bottom": False,
                            "left": col > 0,
                            "right": end_col < width,
                        },
                        "color": color,
                    }

                    tile_info.append(info)
                    tile_id += 1

        # Back edge tiles
        if (depth - tile_d) % step_d != 0:
            d = depth - tile_d
            for row in range(0, height - tile_h + 1, step_h):
                for col in range(0, width - tile_w + 1, step_w):
                    end_row = min(row + tile_h, height)
                    end_col = min(col + tile_w, width)

                    info = {
                        "tile_id": tile_id,
                        "position": (d, row, col),
                        "end_position": (depth, end_row, end_col),
                        "actual_size": None,  # filled below
                        "overlap_region": {
                            "front": True,
                            "back": False,
                            "top": row > 0,
                            "bottom": end_row < height,
                            "left": col > 0,
                            "right": end_col < width,
                        },
                        "color": color,
                    }

                    tile_info.append(info)
                    tile_id += 1

        # Corner edge tiles
        # Bottom-right edge
        if (height - tile_h) % step_h != 0 and (width - tile_w) % step_w != 0:
            row = height - tile_h
            col = width - tile_w
            for d in range(0, depth - tile_d + 1, step_d):
                end_d = min(d + tile_d, depth)

                info = {
                    "tile_id": tile_id,
                    "position": (d, row, col),
                    "end_position": (end_d, height, width),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "front": d > 0,
                        "back": end_d < depth,
                        "top": True,
                        "bottom": False,
                        "left": True,
                        "right": False,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Back-right edge
        if (depth - tile_d) % step_d != 0 and (width - tile_w) % step_w != 0:
            d = depth - tile_d
            col = width - tile_w
            for row in range(0, height - tile_h + 1, step_h):
                end_row = min(row + tile_h, height)

                info = {
                    "tile_id": tile_id,
                    "position": (d, row, col),
                    "end_position": (depth, end_row, width),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "front": True,
                        "back": False,
                        "top": row > 0,
                        "bottom": end_row < height,
                        "left": True,
                        "right": False,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Back-bottom edge
        if (depth - tile_d) % step_d != 0 and (height - tile_h) % step_h != 0:
            d = depth - tile_d
            row = height - tile_h
            for col in range(0, width - tile_w + 1, step_w):
                end_col = min(col + tile_w, width)

                info = {
                    "tile_id": tile_id,
                    "position": (d, row, col),
                    "end_position": (depth, height, end_col),
                    "actual_size": None,  # filled below
                    "overlap_region": {
                        "front": True,
                        "back": False,
                        "top": True,
                        "bottom": False,
                        "left": col > 0,
                        "right": end_col < width,
                    },
                    "color": color,
                }

                tile_info.append(info)
                tile_id += 1

        # Back-bottom-right corner
        if (
            (depth - tile_d) % step_d != 0
            and (height - tile_h) % step_h != 0
            and (width - tile_w) % step_w != 0
        ):
            d = depth - tile_d
            row = height - tile_h
            col = width - tile_w

            info = {
                "tile_id": tile_id,
                "position": (d, row, col),
                "end_position": (depth, height, width),
                "actual_size": None,  # filled below
                "overlap_region": {
                    "front": True,
                    "back": False,
                    "top": True,
                    "bottom": False,
                    "left": True,
                    "right": False,
                },
                "color": color,
            }

            tile_info.append(info)
            tile_id += 1

    else:
        raise ValueError(f"Unsupported image dimensions: {len(shape)}D images are not supported")

    # `actual_size` must match what numpy would actually return. Rederiving from the slice.
    for info in tile_info:
        info["actual_size"] = tuple(
            len(range(*slice(s, e).indices(dim)))
            for s, e, dim in zip(info["position"], info["end_position"], shape)
        )

    return tile_info


def tile_image_with_overlap(
    image: np.ndarray,
    tile_size: int | tuple,
    overlap: int | tuple = 0,
    color: str = None,
):
    """
    Tile a large image into smaller tiles with specified overlap.

    Thin wrapper: :func:`tile_grid_with_overlap` decides the geometry, this slices
    the array accordingly. Behaviour is unchanged — positions, ordering, edge
    tiles and tile_info contents are all produced by the shared function.

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
        See :func:`tile_grid_with_overlap`.

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

    tile_info = tile_grid_with_overlap(image.shape, tile_size, overlap, color=color)
    tiles = [
        image[tuple(slice(s, e) for s, e in zip(i["position"], i["end_position"]))]
        for i in tile_info
    ]
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
    try:
        from matplotlib import colormaps

        cmap = colormaps.get_cmap(colormap)
    except Exception:
        cmap = cm.get_cmap(colormap)
    return lambda x: cmap(norm(x))


def visualize_tiling_grid(tile_info: list, image_shape: tuple, color_prop: str = None):
    """
    Visualize the tiling grid on the image to preview how tiles will be arranged.
    For 3D images, visualizes the middle z-slice.

    Parameters:
    -----------
    tile_info : list
        List of tile information dictionaries
    image_shape : tuple
        Shape of the original image (2D or 3D)
    color_prop : str, optional
        Property name to use for coloring tiles

    Returns:
    --------
    fig : matplotlib figure
        Figure showing the tiling grid
    """
    if len(image_shape) == 3:
        # Handle 3D case - extract middle z-slice and create 2D tile info
        depth, height, width = image_shape
        middle_z = depth // 2
        
        # Filter tiles that include the middle z-slice and convert to 2D
        middle_z_tiles = []
        for tile in tile_info:
            pos = tile["position"]
            end_pos = tile["end_position"]
            
            # Check if this tile includes the middle z-slice
            if len(pos) == 3 and pos[0] <= middle_z < end_pos[0]:
                # Convert 3D tile info to 2D for visualization
                tile_2d = tile.copy()
                tile_2d["position"] = (pos[1], pos[2])  # (row, col)
                tile_2d["end_position"] = (end_pos[1], end_pos[2])
                tile_2d["actual_size"] = (tile["actual_size"][1], tile["actual_size"][2])
                
                # Update overlap region to 2D
                overlap_3d = tile["overlap_region"]
                tile_2d["overlap_region"] = {
                    "top": overlap_3d["top"],
                    "bottom": overlap_3d["bottom"],
                    "left": overlap_3d["left"],
                    "right": overlap_3d["right"],
                }
                
                middle_z_tiles.append(tile_2d)
        
        # Recursively call with 2D image shape and filtered tiles
        return visualize_tiling_grid(middle_z_tiles, (height, width), color_prop)
    
    elif len(image_shape) == 2:
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

        overlap = tile_info[0]["end_position"][1] - tile_info[1]["position"][1] if len(tile_info) > 1 else 0

        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect("equal")
        
        # Update title based on whether this is 3D visualization
        if len(image_shape) == 3:
            depth = image_shape[0]
            middle_z = depth // 2
            ax.set_title(
                f"Tiling Grid Preview (Z-slice {middle_z}/{depth-1})\nImage: {width}x{height}x{depth}, Tile: {tile_w}x{tile_h}, Overlap: {overlap}"
            )
        else:
            ax.set_title(
                f"Tiling Grid Preview\nImage: {width}x{height}, Tile: {tile_w}x{tile_h}, Overlap: {overlap}"
            )
        ax.set_xlabel("Width")
        ax.set_ylabel("Height")

        plt.tight_layout()
        return fig

    else:
        print(f"Visualization only supports 2D and 3D images, got {len(image_shape)}D")
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
    z_index: int | None = None,
):
    """
    Reconstruct the original image from tiles stored in a directory.

    For 3D images, you may pass z_index to reconstruct a single z-slice (2D)
    without allocating the full 3D volume.

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
        Reconstructed image (2D, or 3D if original_shape is 3D and z_index is None)
    """
    output_dir = Path(output_dir)

    if z_index is not None and len(original_shape) != 3:
        raise ValueError("z_index is only valid when original_shape is 3D")

    def _read_tiff_plane(path: Path, plane: int) -> np.ndarray:
        """Read a single z-plane from a TIFF stack if possible.

        Falls back to reading the whole file and slicing if needed.
        """
        try:
            import tifffile  # type: ignore

            return tifffile.imread(path, key=plane)
        except Exception:
            logger.warning(
                "_read_tiff_plane: fast page read failed for %s plane %d — falling back to full load",
                path,
                plane,
            )
            arr = iio.imread(path)
            if arr.ndim == 3:
                return arr[plane]
            return arr

    # Initialize output array and weight map
    reconstructed = None
    weight_map = None
    output_dtype = None  # Captured from the first tile before any histogram matching.

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

        # Load tile (full tile if z_index is None, else only a single plane)
        if z_index is None:
            tile = iio.imread(filepath)
        else:
            pos = info["position"]
            if len(pos) != 3:
                raise ValueError("z_index reconstruction requires 3D tile_info")
            d, _, _ = pos
            actual_size = info["actual_size"]
            if not (d <= z_index < d + actual_size[0]):
                continue
            z_in_tile = int(z_index - d)
            tile = _read_tiff_plane(filepath, z_in_tile)

        # Capture the raw tile dtype before histogram matching may change it.
        if output_dtype is None:
            output_dtype = tile.dtype

        # If adjusting the tile color histogram
        if match_pre:
            if z_index is None:
                input_tile = iio.imread(input_filepath)
                if input_tile.ndim == 3:
                    input_tile = input_tile.mean(axis=0)
                input_tile = input_tile.astype("int16")
                input_tile[input_tile < 0] = 0  # Ensure no negative values
                tile = ski.exposure.match_histograms(tile, input_tile, channel_axis=None)
                del input_tile  # clear memory
            else:
                input_tile = _read_tiff_plane(input_filepath, z_in_tile)
                input_tile = input_tile.astype("int16")
                input_tile[input_tile < 0] = 0
                tile = ski.exposure.match_histograms(tile, input_tile, channel_axis=None)
                del input_tile

        # float32 accumulator avoids uint16 overflow in overlap regions.
        # uint8 weight map is sufficient (max overlap count per pixel is small).
        if reconstructed is None:
            if z_index is None:
                reconstructed = np.zeros(original_shape, dtype=np.float32)
                weight_map = np.zeros(original_shape, dtype=np.uint8)
            else:
                _, height, width = original_shape
                reconstructed = np.zeros((height, width), dtype=np.float32)
                weight_map = np.zeros((height, width), dtype=np.uint8)

        # Place tile in reconstructed image
        pos = info["position"]
        actual_size = info["actual_size"]

        if len(pos) == 2:  # 2D
            r, c = pos
            reconstructed[r : r + actual_size[0], c : c + actual_size[1]] += tile
            weight_map[r : r + actual_size[0], c : c + actual_size[1]] += 1
        else:  # 3D
            d, r, c = pos
            if z_index is None:
                reconstructed[
                    d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
                ] += tile
                weight_map[
                    d : d + actual_size[0], r : r + actual_size[1], c : c + actual_size[2]
                ] += 1
            else:
                reconstructed[r : r + actual_size[1], c : c + actual_size[2]] += tile
                weight_map[r : r + actual_size[1], c : c + actual_size[2]] += 1

    # Handle case where no tiles were found
    if reconstructed is None:
        raise ValueError("No valid tile files found in the output directory")

    # Handle overlapping regions
    if overlap_strategy == "average":
        weight_map[weight_map == 0] = 1
        # In-place divide keeps peak memory at one float32 array + the final uint conversion.
        # (The old `reconstructed / weight_map` created a float32/float64 temporary AND
        # `astype(reconstructed.dtype)` was a no-op copy that together spiked to ~38 GB.)
        np.divide(reconstructed, weight_map, out=reconstructed)
        del weight_map
    elif overlap_strategy == "max":
        # TODO: For max strategy, we need to reconstruct differently
        pass  # Implementation would be more complex

    reconstructed[reconstructed < 0] = 0
    # Convert float32 accumulator back to the original tile dtype (e.g. uint16).
    # This halves the output TIFF size compared to the previous float32 output.
    return reconstructed.astype(output_dtype)


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


#### HEX GRID TILING FUNCTIONS ####
def create_hexagonal_grid(bounds, hex_size, overlap=0.0):
    """
    Create a hexagonal grid covering the given bounds using Red Blob Games standard geometry.
    
    For pointy-top hexagons:
    - Horizontal spacing = sqrt(3) * size
    - Vertical spacing = 3/2 * size
    
    Parameters:
    -----------
    bounds : tuple
        (minx, miny, maxx, maxy) bounding box to cover
    hex_size : float
        Radius of hexagon (distance from center to vertex)
    overlap : float, default 0.0
        Overlap parameter:
        - 0.0: No overlap (disjoint hexagons)
        - 0.5: 50% overlap 
        - 1.0: Complete overlap (same position)
        - Negative values create gaps
    
    Returns:
    --------
    geopandas.GeoDataFrame
        Grid of hexagonal polygons
    """
    from shapely import Polygon, Point, box

    minx, miny, maxx, maxy = bounds
    
    # Red Blob Games standard spacing for pointy-top hexagons
    horizontal_spacing = hex_size * np.sqrt(3)  # sqrt(3) * size
    vertical_spacing = hex_size * 1.5           # 3/2 * size
    
    # Apply overlap: spacing = base_spacing * (1 - overlap)
    # overlap=0.0 -> spacing = base_spacing (no overlap)
    # overlap=0.5 -> spacing = 0.5 * base_spacing (50% overlap)
    # overlap=1.0 -> spacing = 0 (complete overlap)
    actual_horizontal_spacing = horizontal_spacing * (1 - overlap)
    actual_vertical_spacing = vertical_spacing * (1 - overlap)
    
    # Calculate grid dimensions
    width = maxx - minx
    height = maxy - miny
    
    # Number of hexagons needed (with some buffer)
    cols = int(np.ceil(width / actual_horizontal_spacing)) + 2
    rows = int(np.ceil(height / actual_vertical_spacing)) + 2
    
    hexagons = []
    
    for row in range(rows):
        for col in range(cols):
            # Calculate center position
            # Even rows: no horizontal offset
            # Odd rows: offset by half the horizontal spacing
            if row % 2 == 0:
                x = minx + col * actual_horizontal_spacing
            else:
                x = minx + (col + 0.5) * actual_horizontal_spacing
            
            y = miny + row * actual_vertical_spacing
            
            # Create hexagon geometry (pointy-top orientation)
            angles = np.linspace(0, 2 * np.pi, 7)  # 7 points to close the polygon
            # For pointy-top: first vertex at 30° (π/6 radians)
            angles = angles + np.pi/6
            
            hex_x = x + hex_size * np.cos(angles)
            hex_y = y + hex_size * np.sin(angles)
            
            hex_coords = list(zip(hex_x, hex_y))
            hexagon = Polygon(hex_coords)
            
            hexagons.append({
                'geometry': hexagon,
                'row': row,
                'col': col,
                'center_x': x,
                'center_y': y
            })
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(hexagons)
    
    # Filter to only hexagons that intersect with bounds
    bounds_poly = box(minx, miny, maxx, maxy)
    gdf = gdf[gdf.geometry.intersects(bounds_poly)].copy()
    gdf.reset_index(drop=True, inplace=True)
    
    return gdf


def hex_area_from_size(hex_size: float) -> float:
    """
    Calculate area of a regular hexagon given the distance from center to a vertex.

    Parameters:
    - hex_size: float, distance from center to any vertex

    Returns:
    - area: float
    """
    if hex_size <= 0:
        raise ValueError("hex_size must be > 0")
    return (3.0 * math.sqrt(3.0) / 2.0) * (hex_size ** 2)


# ===========================================================================
# OVERLAPPING-TILE GEOMETRY AND TILE-SIZE OPTIMISATION
#
# optimize_grid picks the (tile_y, tile_x) that minimises area for a given tissue mask, to save compute cost.
# ===========================================================================

def bbox_of_mask(mask: np.ndarray) -> tuple[int, int, int, int]:
    """
    (y0, y1, x0, x1) bounding box of any 2D boolean mask, half-open.

    Used to shrink a tiling extent to the region that actually has content.
    Raises if the mask is empty — tiling nothing is never intended.
    """
    rows = np.flatnonzero(mask.any(axis=1))
    cols = np.flatnonzero(mask.any(axis=0))
    if rows.size == 0 or cols.size == 0:
        raise ValueError("mask is empty; nothing to tile")
    return int(rows[0]), int(rows[-1]) + 1, int(cols[0]), int(cols[-1]) + 1


@dataclass(frozen=True)
class Tile:
    """
    One tile of an overlapping grid, in absolute pixel coordinates.

    ``gi``/``gj`` are grid indices, which let a consumer walk neighbours; 
    ``y0``/``x0`` are absolute offsets, they already include the bbox origin.
    """
    gi: int
    gj: int
    y0: int
    x0: int
    height: int
    width: int

    @property
    def y1(self) -> int:
        return self.y0 + self.height

    @property
    def x1(self) -> int:
        return self.x0 + self.width

    @property
    def area(self) -> int:
        return self.height * self.width


def build_grid(bbox: tuple[int, int, int, int], tile, overlap: int) -> list[Tile]:
    """
    Tile a bounding box.

    The layout rules (stride, right-aligned final tile, per-axis solving) live in
    :func:`tile_grid_with_overlap` above.
    """
    y0, y1, x0, x1 = bbox
    height, width = y1 - y0, x1 - x0
    tile_y, tile_x = (tile, tile) if isinstance(tile, (int, np.integer)) else tile
    tile_y, tile_x = min(int(tile_y), height), min(int(tile_x), width)
    if tile_y <= overlap or tile_x <= overlap:
        raise ValueError(
            f"tile ({tile_y}x{tile_x}, clamped to the {height}x{width} extent) "
            f"must exceed overlap ({overlap})"
        )

    info = tile_grid_with_overlap((height, width), (tile_y, tile_x), overlap)

    # tile_info has no grid indices: the shared function emits the main grid
    # first and then appends edge rows/columns, so rank the distinct origins.
    ys = sorted({i["position"][0] for i in info})
    xs = sorted({i["position"][1] for i in info})
    yi = {v: k for k, v in enumerate(ys)}
    xi = {v: k for k, v in enumerate(xs)}
    return [
        Tile(
            gi=yi[i["position"][0]], gj=xi[i["position"][1]],
            y0=y0 + i["position"][0], x0=x0 + i["position"][1],
            height=i["actual_size"][0], width=i["actual_size"][1],
        )
        for i in info
    ]

def redundancy(tiles: list[Tile], bbox: tuple[int, int, int, int]) -> float:
    """Total tile area over bbox area — the fraction of pixels processed twice."""
    y0, y1, x0, x1 = bbox
    return sum(t.area for t in tiles) / ((y1 - y0) * (x1 - x0))


def tile_coverage(mask: np.ndarray, tile: Tile) -> float:
    """Fraction of a tile's pixels that are tissue."""
    window = mask[tile.y0:tile.y1, tile.x0:tile.x1]
    if window.size == 0:
        return 0.0
    return float(window.mean())


# ---------------------------------------------------------------------------
# tile-size optimisation
#
# Tiles need not be square, and the cost surface is not smooth: an arbitrary tile
# does not pack evenly into the extent, so the adaptive grid appends a
# right-aligned final tile whose overlap exceeds `o`. That makes redundancy
# non-monotonic in tile size. There is therefore no analytic optimum to solve for.
# The surface is cheap to evaluate exhaustively.
# ---------------------------------------------------------------------------
# The objective is simply TOTAL AREA: sum of the areas of tiles that contain any tissue.

# Default search bounds per axis, in full-resolution pixels.
# tile_min should be at least poly_tile_size * scale (see `CellposeConfig.min_gpu_tile`)
DEFAULT_TILE_MIN = 1500
DEFAULT_TILE_MAX = 12000

def mask_integral_image(
    mask: np.ndarray,
    ds: int = 8
) -> tuple[np.ndarray, int]:
    """
    Summed-area table of tissue pixels over a ``ds``-downsampled mask.

    Downsampling keeps the table small and costs no meaningful accuracy.

    Returned table is padded by one zero row/column so rectangle sums need no
    boundary special-casing.
    """
    if ds < 1:
        raise ValueError(f"ds must be >= 1, got {ds}")
    h, w = mask.shape
    nh, nw = h // ds, w // ds
    if nh == 0 or nw == 0:
        raise ValueError(f"ds={ds} is too coarse for a {h}x{w} mask")
    reduced = (mask[: nh * ds, : nw * ds]
               .reshape(nh, ds, nw, ds)
               .sum(axis=(1, 3), dtype=np.int64))
    integral = np.zeros((nh + 1, nw + 1), dtype=np.int64)
    np.cumsum(np.cumsum(reduced, axis=0), axis=1, out=integral[1:, 1:])
    return integral, ds


def rect_sums(
    integral: np.ndarray,
    ds: int,
    y0, y1, x0, x1
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorised tissue-pixel sums for rectangles given in full-res coordinates.

    Returns ``(tissue_pixels, sampled_pixels)``. Bounds are snapped outward to
    block edges (floor the lower, ceil the upper), so the sampled rectangle is a
    superset of the requested one. ``sampled_pixels`` is returned so callers can judge.
    """
    limit_y, limit_x = integral.shape[0] - 1, integral.shape[1] - 1
    a = np.clip(np.asarray(y0) // ds, 0, limit_y)
    b = np.clip(-(-np.asarray(y1) // ds), 0, limit_y)   # ceil
    c = np.clip(np.asarray(x0) // ds, 0, limit_x)
    d = np.clip(-(-np.asarray(x1) // ds), 0, limit_x)
    sums = integral[b, d] - integral[a, d] - integral[b, c] + integral[a, c]
    sampled = (b - a) * ds * (d - c) * ds
    return sums, sampled


def score_grid(
    integral: np.ndarray,
    ds: int,
    bbox: tuple[int, int, int, int],
    tile_y: int,
    tile_x: int,
    overlap: int,
) -> dict:
    """
    Cost of one candidate ``(tile_y, tile_x, overlap)``, via the integral image.

    ``cost_area`` is the total pixel area of tiles that contain tissue — the
    quantity to minimise.
    """
    y0, y1, x0, x1 = bbox
    tiles = build_grid(bbox, (tile_y, tile_x), overlap)
    gy = np.array([t.y0 for t in tiles])
    gx = np.array([t.x0 for t in tiles])
    height = tiles[0].height
    width = tiles[0].width

    sums, sampled = rect_sums(integral, ds, gy, gy + height, gx, gx + width)

    n_tiles = gy.size
    area = height * width
    kept_mask = sums > 0
    n_kept = int(kept_mask.sum())

    # `sums` counts tissue pixels (each ds x ds block contributes up to ds^2), so
    # "fully tissue" is `sums == sampled`, both in pixels.
    n_partial = int((kept_mask & (sums < sampled)).sum())
    return {
        "tile_y": tile_y,
        "tile_x": tile_x,
        "overlap": overlap,
        "grid": (len({t.gi for t in tiles}), len({t.gj for t in tiles})),
        "n_tiles": n_tiles,
        "n_kept": n_kept,
        "n_partial": n_partial,
        "tile_area": area,
        "redundancy": (n_tiles * area) / ((y1 - y0) * (x1 - x0)),
        # the objective: total area actually pushed through the GPU
        "cost_area": n_kept * area,
    }


def optimize_grid(
    mask: np.ndarray,
    *,
    overlap: int | list[int] = 500,
    tile_min: int = DEFAULT_TILE_MIN,
    tile_max: int = DEFAULT_TILE_MAX,
    step: int = 250,
    max_tile_area: int | None = None,
    ds: int = 8,
    allow_non_square: bool = True,
) -> tuple[dict, list[dict]]:
    """Search the 2D ``(tile_y, tile_x)`` surface and return the cheapest grid.

    Exhaustive rather than analytic because the surface is discontinuous (see the
    module comment above). Costs one integral-image pass plus one vectorised
    lookup per candidate.

    Parameters:
    mask (np.ndarray): boolean tissue mask at full resolution.
    overlap (int | list[int]): overlap(s) to search. Must stay above roughly one
        cell diameter (~95 px at 0.106 um/px) so a boundary cell is fully
        contained in at least one tile — this is a correctness floor, not a
        tunable, so it is not searched downward by default.
    tile_min / tile_max / step (int): candidate range per axis, in full-res px.
    max_tile_area (int | None): reject candidates whose tile area exceeds this.
        Host RAM is the real tile-size limit and scales with area.
    allow_non_square (bool): if False, restrict to ``tile_y == tile_x``.

    Returns ``(best, all_scored)``; ``all_scored`` is sorted by cost so callers
    can log the runner-up surface or record it for provenance.
    """
    integral, ds = mask_integral_image(mask, ds)
    bbox = bbox_of_mask(mask)
    y_extent, x_extent = bbox[1] - bbox[0], bbox[3] - bbox[2]
    overlaps = [overlap] if isinstance(overlap, int) else list(overlap)

    def candidates(extent: int) -> list[int]:
        top = min(tile_max, extent)
        vals = list(range(tile_min, top + 1, step))
        if not vals or vals[-1] != top:
            vals.append(top)                 # always consider "one tile per axis"
        return vals

    ys_c, xs_c = candidates(y_extent), candidates(x_extent)
    scored: list[dict] = []
    for o in overlaps:
        for ty in ys_c:
            if ty <= o:
                continue
            for tx in (xs_c if allow_non_square else [ty]):
                if tx <= o or tx > x_extent and tx != min(tile_max, x_extent):
                    continue
                if max_tile_area is not None and ty * tx > max_tile_area:
                    continue
                scored.append(score_grid(integral, ds, bbox, ty, tx, o))

    if not scored:
        raise ValueError(
            "no viable tile candidates; check tile_min/tile_max/overlap and "
            "max_tile_area against the tissue bbox "
            f"({y_extent}x{x_extent})"
        )
    scored.sort(key=lambda s: s["cost_area"])
    best = scored[0]
    logger.info(
        "optimize_grid: searched %d candidates -> best tile %dx%d overlap=%d, "
        "grid %dx%d, %d/%d tiles kept (%d partial), redundancy=%.2f, cost=%.0f GPU-s",
        len(scored), best["tile_y"], best["tile_x"], best["overlap"],
        *best["grid"], best["n_kept"], best["n_tiles"], best["n_partial"],
        best["redundancy"], best["cost_area"],
    )
    square_best = next((s for s in scored
                        if s["tile_y"] == s["tile_x"]), None)
    if square_best is not None and square_best is not best:
        logger.info(
            "optimize_grid: best SQUARE tile would be %d (area %.3g px^2, "
            "%.1f%% worse than the non-square optimum)",
            square_best["tile_y"], square_best["cost_area"],
            (square_best["cost_area"] / best["cost_area"] - 1) * 100,
        )
    return best, scored
