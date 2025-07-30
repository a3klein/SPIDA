import os
from dotenv import load_dotenv  # type: ignore

import numpy as np
from pathlib import Path
import warnings
import logging

import multiprocessing as mp
from functools import partial
from tqdm import tqdm
import itertools

from natsort import natsorted
import skimage as ski
from shapely.geometry import Polygon
import geopandas as gpd
from spida.utilities.tiling import (
    tile_image_with_overlap,
    merge_tile_polygons,
)  # , merge_overlapping_polygons

from cellpose import models, core, io  # type: ignore

io.logger_setup()
load_dotenv()

if not core.use_gpu():
    raise ImportError("No GPU access")

logger = logging.getLogger(__package__)


def _load_image(
    image_path: Path,
    image_ext: str = ".tif",
    nuc_stain_name: str = "DAPI",
    cyto_stain_name: str = None,
):
    """
    Load an image from the specified path.

    Parameters:
    image_path (str): Path to the image file.
    image_ext (str): Extension of the image file (default is '.tif').
    nuc_stain_name (str): Name of the nuclear stain (default is 'DAPI').
    cyto_stain_name (str): Name of the cytoplasmic stain (default is 'PolyT').

    Returns:
    np.ndarray: Loaded image as a numpy array.
    """
    files = natsorted(
        [
            f
            for f in image_path.glob("*" + image_ext)
            if "_masks" not in f.name and "_flows" not in f.name
        ]
    )
    logger.info(f"Found {len(files)} images in {image_path}")

    nuc_file = None
    for f in files:
        if nuc_stain_name in f.name:
            nuc_file = f

    if nuc_file is None:
        raise ValueError(
            f"No nuclear stain file found with name '{nuc_stain_name}' in {image_path}"
        )

    nuc_img = io.imread(nuc_file)
    logger.info(f"Loaded images: Nuclear - {nuc_file}, Image Shape: {nuc_img.shape}")

    if cyto_stain_name is not None:
        for f in files:
            if cyto_stain_name in f.name:
                cyto_file = f

        cyto_img = io.imread(cyto_file)
        logger.info(
            f"Loaded images: Cytoplasmic - {cyto_file}, Image Shape: {cyto_img.shape}"
        )
        img = np.stack([nuc_img, cyto_img], axis=0)
        return img

    else:
        logger.info(
            f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}, returning nuclear image only."
        )
        return np.expand_dims(nuc_img, axis=0)


def _load_model(model_name: str = "cpsam", **kwargs):
    """
    Load The Cellpose Model
    """
    model_path = os.getenv(
        "CELLPOSE_MODEL_PATH",
        "/anvil/projects/x-mcb130189/aklein/programs/.cellpose/models",
    )
    model = models.CellposeModel(
        pretrained_model=f"{model_path}/{model_name}", gpu=True
    )
    return model


def _cellpose_wrapper(
    img,
    batch_size=8,
    flow_threshold=0.0,
    cellprob_threshold=-4,
    tile_norm_blocksize=0,
    diameter=None,
    min_size=100,
    **kwargs,
):
    """
    Wrapper for Cellpose model evaluation.

    Parameters:
    img (np.ndarray): Input image.
    flow_threshold (float): Flow threshold for Cellpose.
    cellprob_threshold (float): Cell probability threshold for Cellpose.
    tile_norm_blocksize (int): Normalization block size for Cellpose.
    batch_size (int): Batch size for Cellpose evaluation.

    Returns:
    tuple: Masks, flows, and styles from the Cellpose model.
    """

    model = _load_model(**kwargs)

    return model.eval(
        img,
        batch_size=batch_size,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        channel_axis=0,
        min_size=min_size,
        normalize={"tile_norm_blocksize": tile_norm_blocksize},
    )


def convert_mask(mask_id, masks, tolerance=0.5):
    """Convert a single mask to a polygon.
    Parameters
    ----------
    mask_id : int
        A single mask to convert.
    masks : np.ndarray
        A (Y,X) or (Z,Y,X) numpy array of segmentation masks.
    tolerance : float, optional
        The maximum allowed geometry displacement. The higher this value, the smaller the number of vertices in the resulting geometry.
    Returns
    -------
    shapely.geometry.Polygon
        The converted polygon.
    """

    polygons = []
    cells = []
    layers = []
    mask = masks == mask_id
    # Do this across all z layers
    for m in np.where(mask.any(axis=(1, 2)))[0]:
        # Pad mask to ensure contours are closed
        padded_mask = np.pad(mask[m], 1, mode="constant")
        contours = ski.measure.find_contours(padded_mask, 0.5)
        contours = [contour - 1 for contour in contours]
        if len(contours[0]) > 3:
            polygons.append(
                Polygon(contours[0]).simplify(
                    tolerance=tolerance, preserve_topology=True
                )
            )
            cells.append(mask_id)
            layers.append(m)
    return polygons, cells, layers


def _masks_to_polygons(masks: np.ndarray, tolerance: float = 0.5):
    """Convert Cellpose masks to polygons.

    Parameters:
    masks (np.ndarray): Cellpose segmentation masks.
    tolerance (float): Tolerance for polygon simplification.

    Returns:
    tuple: List of polygons, list of cell IDs, list of z layers.
    """
    mask_ids = np.unique(masks)[1:]

    masks = np.expand_dims(masks, axis=0)
    masks = masks.swapaxes(-1, -2)
    mask_ids = np.unique(masks)[1:]  # Exclude background (0)

    polygons = []
    for mask_id in mask_ids:
        polygons_batch, cells_batch, layers_batch = convert_mask(
            mask_id, masks, tolerance=tolerance
        )
        polygons.append((polygons_batch, cells_batch, layers_batch))

    # Or using itertools.chain for efficiency with large datasets:
    all_polygons = list(
        itertools.chain.from_iterable(
            polygons_batch for polygons_batch, _, _ in polygons
        )
    )
    all_cells = list(
        itertools.chain.from_iterable(cells_batch for _, cells_batch, _ in polygons)
    )
    all_layers = list(
        itertools.chain.from_iterable(layers_batch for _, _, layers_batch in polygons)
    )

    gdf = gpd.GeoDataFrame(
        {"ID": all_cells, "Geometry": all_polygons, "z": all_layers},
        geometry="Geometry",
        crs=None,
    )
    # logging.info(f"Converting {len(mask_ids)} masks to polygons; gdf shape: {gdf.shape}.")
    return gdf


# NEW USING SKIMAGE FIND CONTOURS
def masks_to_geodataframe(
    masks: np.ndarray,
    tolerance: float = 0.5,
):
    """Convert Cellpose masks to a GeoDataFrame with polygons."""
    tiles, tile_info = tile_image_with_overlap(masks, tile_size=1000, overlap=200)
    logger.info(f"Generated {len(tiles)} tiles.")

    parallel_func = partial(
        _masks_to_polygons,
        tolerance=tolerance,  # Adjust tolerance as needed
    )

    with mp.Pool(max_cpu) as pool:
        geo_list = list(tqdm(pool.map(parallel_func, tiles), total=len(tiles)))

    logger.info("Finished Conversion")

    # Merge all tiles into global coordinates
    merged_polygons = merge_tile_polygons(tile_info, geo_list, geom_col="Geometry")

    final_polygons = merged_polygons.dissolve(by="ID").reset_index()
    # # Optional: merge overlapping polygons across tile boundaries
    # final_polygons = merge_overlapping_polygons(merged_polygons)

    return final_polygons


def run_cellposeSAM(
    root_dir: str,
    output_dir: str,
    region: str,
    scale: int = 4,  # value to downsample image by
    min_size: int = 100,  # minimum size of cell to keep
    image_ext: str = ".tif",
    nuc_stain_name: str = "DAPI",
    cyto_stain_name: str = None,
    **kwargs,
):
    """
    Run Cellpose SAM segmentation on a specified region.
    Parameters:
    root_dir (str): The root directory containing the images.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    scale (int): The factor by which to downscale the image (default is 4).
    min_size (int): Minimum size of cells to keep in the segmentation (default is 100).
    image_ext (str): Extension of the image files (default is '.tif').
    nuc_stain_name (str): Name of the nuclear stain (default is 'DAPI').
    cyto_stain_name (str): Name of the cytoplasmic stain (default is 'PolyT').
    """

    logger.info(f"root_dir = {root_dir}")
    logger.info(f"output_dir = {output_dir}")
    logger.info(f"region = {region}")
    logger.info(f"scale = {scale}")
    logger.info(f"min_size = {min_size}")
    logger.info(f"image_ext = {image_ext}")
    logger.info(f"nuc_stain_name = {nuc_stain_name}")
    logger.info(f"cyto_stain_name = {cyto_stain_name}")
    for key, value in kwargs.items():
        logger.info(f"Cellpose SAM parameter: {key} = {value}")

    global max_cpu
    max_cpu = mp.cpu_count()
    logger.info(f"Max CPU count: {max_cpu}")

    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    # data_path =f"{root_dir}/{region}/images/"
    data_path = root_dir

    img = _load_image(
        Path(data_path),
        image_ext=image_ext,
        nuc_stain_name=nuc_stain_name,
        cyto_stain_name=cyto_stain_name,
    )

    og_shape = img.shape
    logger.info(f"OG image shape: {og_shape}")

    downscale = (1,) * (img.ndim - 2) + (scale, scale)
    img = ski.transform.downscale_local_mean(img, downscale)
    logger.info(f"Image shape after scaling: {img.shape}")

    logger.info("STARTING SEGMENTATION")
    masks, flows, styles = _cellpose_wrapper(img, batch_size=16, **kwargs)
    logger.info("SEGMENTATION COMPLETED")
    # ski.io.imsave(f"{output_dir}/{region}/masks.tif", masks.astype(np.uint16))

    logger.info(f"Cellpose segmentation completed for region {region}.")
    logger.info(
        f"Number of masks detected: {len(np.unique(masks)) - 1}; mask shape: {masks.shape}"
    )

    # rescaled_masks = cv2.resize(masks, (og_shape[1], og_shape[0]), interpolation=cv2.INTER_NEAREST)
    # print("Rescaled masks shape:", rescaled_masks.shape)

    # masks to geodataframe
    gdf = masks_to_geodataframe(masks)
    gdf.geometry = gdf.geometry.affine_transform(
        [scale, 0, 0, scale, 0, 0]
    )  # upscale polygons to original image size

    # Filter out small cells by size < min_size:
    logger.info(f"Filtering out cells with area < {min_size}")
    n_pre = len(gdf)
    gdf = gdf[gdf.geometry.area > min_size]
    n_post = len(gdf)
    logger.info(f"Filtered out {n_pre - n_post} cells; remaining cells: {n_post}")

    # Save masks and flows
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf.to_parquet(f"{output_dir}/{region}/polygons.parquet", index=False)
