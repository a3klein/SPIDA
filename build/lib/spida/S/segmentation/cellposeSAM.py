import os
from dotenv import load_dotenv  # type: ignore

import numpy as np
import pandas as pd
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
import cv2
from spida.utilities.tiling import (
    tile_image_with_overlap,
    merge_tile_polygons,
)  # , merge_overlapping_polygons

from cellpose import models, core, io  # type: ignore

# io.logger_setup() # TODO: Maybe put this in a temporary file since this creates race conditions between multiple runs. 
load_dotenv()

logger = logging.getLogger(__package__)


def _require_gpu():
    if not core.use_gpu():
        raise ImportError("No GPU access")


def _load_image(
    image_path: Path,
    image_ext: str = ".tif",
    nuc_stain_name: str = "DAPI",
    cyto_stain_name: str = None,
    downscale: int | None = None,
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

    nuc_files = []
    for f in files:
        if nuc_stain_name in f.name:
            nuc_files.append(f)

    if len(nuc_files) == 0:
        raise ValueError(
            f"No nuclear stain file found with name '{nuc_stain_name}' in {image_path}"
        )
    elif len(nuc_files) == 1: 
        nuc_file = nuc_files[0]
        nuc_img = io.imread(nuc_file)
        if downscale and downscale > 1:
            factors = (1,) * (nuc_img.ndim - 2) + (downscale, downscale)
            nuc_img = ski.transform.downscale_local_mean(nuc_img, factors)
    else: 
        # Stack 3D image
        nuc_files = natsorted(nuc_files)
        stack_3d = []
        for _fn in nuc_files:
            logger.info(f"  {_fn}")
            _img = io.imread(_fn)
            if downscale and downscale > 1:
                _img = ski.transform.downscale_local_mean(
                    _img, (downscale, downscale)
                )
            stack_3d.append(_img)
        nuc_img = np.stack(stack_3d, axis=-1)

    logger.info(f"Loaded images: Nuclear - {nuc_files}, Image Shape: {nuc_img.shape}")

    if cyto_stain_name is not None:
        cyto_files = []
        for f in files:
            if cyto_stain_name in f.name:
                cyto_files.append(f)

        if len(cyto_files) == 0:
            raise ValueError(
                f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}"
            )
        elif len(cyto_files) == 1: 
            cyto_file = cyto_files[0]
            cyto_img = io.imread(cyto_file)
            if downscale and downscale > 1:
                factors = (1,) * (cyto_img.ndim - 2) + (downscale, downscale)
                cyto_img = ski.transform.downscale_local_mean(cyto_img, factors)
        else: 
            # Stack 3D image
            cyto_files = natsorted(cyto_files)
            stack_3d = []
            for _fn in cyto_files:
                logger.info(f"  {_fn}")
                _img = io.imread(_fn)
                if downscale and downscale > 1:
                    _img = ski.transform.downscale_local_mean(
                        _img, (downscale, downscale)
                    )
                stack_3d.append(_img)
            cyto_img = np.stack(stack_3d, axis=-1)

        logger.info(f"Loaded images: Cytoplasmic - {cyto_files}, Image Shape: {cyto_img.shape}")
        return np.stack([nuc_img, cyto_img], axis=0)
    else:
        logger.info(
            f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}, returning nuclear image only."
        )
        return np.expand_dims(nuc_img, axis=0)


def _load_model(model_name: str = "cpsam", **kwargs):
    """
    Load The Cellpose Model
    """
    _require_gpu()
    model_path = os.getenv(
        "CELLPOSE_MODEL_PATH",
        "/anvil/projects/x-mcb130189/aklein/programs/.cellpose/models",
    )
    model = models.CellposeModel(
        pretrained_model=f"{model_path}/{model_name}",
        gpu=True,
        use_bfloat16=False
    )
    return model


def _cellpose_wrapper(
    img,
    batch_size=8,
    flow_threshold=0.0,
    cellprob_threshold=-2,
    tile_norm_blocksize=0,
    diameter=None,
    min_size=100,
    do_3D=False,
    stitch_threshold: float | None = None,
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

    if do_3D:
        if img.ndim != 4:
            raise ValueError(
                f"Unsupported image shape for 3D: {img.shape}. Expected (C, Y, X, Z)."
            )
        # expected input for 3D with channels last: (Z, Y, X, C)
        img = np.transpose(img, (3, 1, 2, 0))
        channel_axis = -1
        z_axis = 0
        logger.info(f"do_3D=True; image dimensions = {img.shape}, channel_axis={channel_axis}, z_axis={z_axis}")

        return model.eval(
            img,
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            diameter=diameter,
            channel_axis=channel_axis,
            min_size=min_size,
            do_3D=True,
            z_axis=z_axis,
            normalize={"tile_norm_blocksize": tile_norm_blocksize},
        )
    else: 
        z_axis = None
        if img.ndim == 4:
            # Treat as 2D stack with channels: (C, Y, X, Z) -> (Z, Y, X, C)
            img = np.transpose(img, (3, 1, 2, 0))
            channel_axis = -1
            z_axis = 0
        elif img.ndim == 3:
            # 2D with channels-first: (C, Y, X) -> (Y, X, C)
            img = np.moveaxis(img, 0, -1)
            channel_axis = -1
            z_axis = None
            stitch_threshold = 0
        elif img.ndim == 2:
            channel_axis = None
            z_axis = None
            stitch_threshold = 0
        else:
            raise ValueError("Unsupported image shape for 2D segmentation.")
        
        logger.info(f"do_3D=False; image dimensions = {img.shape}, channel_axis={channel_axis}, z_axis={z_axis}")

        eval_kwargs = dict(
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            diameter=diameter,
            channel_axis=channel_axis,
            min_size=min_size,
            normalize={"tile_norm_blocksize": tile_norm_blocksize},
            stitch_threshold=stitch_threshold,
            do_3D=False,
        )
        if z_axis is not None:
            eval_kwargs["z_axis"] = z_axis

        logger.info(f"Running Cellpose Model Evaluation with parameters:")
        for key, value in eval_kwargs.items():
            logger.info(f"  {key}: {value}")
        return model.eval(img, **eval_kwargs)


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

    if masks.ndim == 2: 
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
    tile_size: int | tuple = 1000,
    overlap: int | tuple = 200,
):
    """Convert Cellpose masks to a GeoDataFrame with polygons."""
    tiles, tile_info = tile_image_with_overlap(
        masks, tile_size=tile_size, overlap=overlap
    )
    logger.info(f"Generated {len(tiles)} tiles.")

    parallel_func = partial(
        _masks_to_polygons,
        tolerance=tolerance,  # Adjust tolerance as needed
    )

    with mp.Pool(max_cpu) as pool:
        geo_list = list(tqdm(pool.map(parallel_func, tiles), total=len(tiles)))

    # Merge all tiles into global coordinates
    merged_polygons = merge_tile_polygons(tile_info, geo_list, geom_col="Geometry")

    logger.info("Finished Conversion")

    return merged_polygons

    # Took the final dissolve step outside of the merge tiles because it is N-D dependent
    # final_polygons = merged_polygons.dissolve(by="ID").reset_index()
    
    # Optional: merge overlapping polygons across tile boundaries
    # --> Does not work with 3D stacks, however dissolve is already doing a good job of merging across tile boundaries, so may not be necessary.
    # final_polygons = merge_overlapping_polygons(merged_polygons)

    # return final_polygons

def run_cellposeSAM(
    root_dir: str,
    output_dir: str,
    region: str,
    scale: int = 4,  # value to downsample image by
    min_size: int = 100,  # minimum size of cell to keep
    image_ext: str = ".tif",
    nuc_stain_name: str = "DAPI",
    cyto_stain_name: str = None,
    project_3d_to_2d: bool = False,
    stitch_threshold: float | None = 0.2,
    apply_clahe: bool = False,
    do_3D: bool | None = None,
    micron_per_z : float = 1.5,
    min_z : int = 3,  # minimum number of z layers a cell must appear in to be kept (only for 3D stacks)
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
        stitch_threshold (float | None): Cellpose stitch threshold for 2D stacks (default 0.4).
        apply_clahe (bool): Whether to apply CLAHE contrast enhancement (default False).
        do_3D (bool | None): Whether to perform 3D segmentation in cellpose (default None, will use 2D segmentation + stitching).
        project_3d_to_2d (bool): Whether to project 3D image to 2D using max intensity projection before segmentation (default False).
        micron_per_z (float): Microns per z layer for 3D data (default 1.5).
        min_z (int): Minimum number of z layers a cell must appear in to be kept (only for 3D stacks; default 3).
    **kwargs: Additional keyword arguments to pass to the Cellpose model.
    Returns:
        bool: True if 3D segmentation was performed, False if 2D segmentation was performed.
    """

    logger.info(f"root_dir = {root_dir}")
    logger.info(f"output_dir = {output_dir}")
    logger.info(f"region = {region}")
    logger.info(f"scale = {scale}")
    logger.info(f"min_size = {min_size}")
    logger.info(f"image_ext = {image_ext}")
    logger.info(f"nuc_stain_name = {nuc_stain_name}")
    logger.info(f"cyto_stain_name = {cyto_stain_name}")
    logger.info(f"project_3d_to_2d = {project_3d_to_2d}")
    logger.info(f"stitch_threshold = {stitch_threshold}")
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
        downscale=scale,
    )

    og_shape = img.shape
    logger.info(f"OG image shape: {og_shape}")

    if project_3d_to_2d and img.ndim == 4:
        logger.info("Converting 3D image to 2D by maximum intensity projection.")
        img = np.max(img, axis=1)
        logger.info(f"Image shape after 3D to 2D conversion: {img.shape}")

    if apply_clahe and img.ndim == 3:
        clahe_clipLimit = kwargs.pop("clahe_clipLimit", 20)
        clahe_tileGridSize = kwargs.pop("clahe_tileGridSize", 8)
        logger.info(
            f"Applying CLAHE with clipLimit={clahe_clipLimit} and tileGridSize={clahe_tileGridSize}."
        )
        clahe = cv2.createCLAHE(
            clipLimit=clahe_clipLimit,
            tileGridSize=(clahe_tileGridSize, clahe_tileGridSize)
        )
        img = img.astype(np.uint16)
        for _channels in range(img.shape[0]):
            img[_channels] = clahe.apply((img[_channels]))
        img = (img/256).astype('uint8')

    logger.info(f"Image shape after scaling: {img.shape}")

    do_3d = do_3D
    if do_3d is None:
        do_3d = kwargs.pop("do_3D", False)
    if project_3d_to_2d:
        do_3d = False
    logger.info(f"STARTING SEGMENTATION (do_3D={do_3d})")
    if do_3d:
        masks, flows, styles = _cellpose_wrapper(
            img, batch_size=16, do_3D=True, **kwargs
        )
    else:
        masks, flows, styles = _cellpose_wrapper(
            img,
            batch_size=16,
            do_3D=False,
            stitch_threshold=stitch_threshold,
            **kwargs,
        )
    logger.info("SEGMENTATION COMPLETED")
    # ski.io.imsave(f"{output_dir}/{region}/masks.tif", masks.astype(np.uint16))

    logger.info(f"Cellpose segmentation completed for region {region}.")
    logger.info(
        f"Number of masks detected: {len(np.unique(masks)) - 1}; mask shape: {masks.shape}"
    )

    # rescaled_masks = cv2.resize(masks, (og_shape[1], og_shape[0]), interpolation=cv2.INTER_NEAREST)
    # print("Rescaled masks shape:", rescaled_masks.shape)

    is_3d = img.ndim == 4
    logger.info(f"Image dimensionality: {img.ndim}D; treating as {'3D' if is_3d else '2D'} for polygon conversion.")

    tile_size = kwargs.pop("tile_size", 1000)
    overlap = kwargs.pop("overlap", 200)

    # masks to geodataframe
    gdf = masks_to_geodataframe(masks, tile_size=tile_size, overlap=overlap)
    if is_3d: 
        logger.info("Adding ZIndex, ZLevel, and EntityID columns to GeoDataFrame for 3D data.")
        gdf['ZIndex'] = gdf['z']
        gdf['ZLevel'] = gdf['z'] * micron_per_z  # Assuming 1.5 microns per z layer
        gdf['UID'] = gdf['tile_id'].astype(str) + '-' + gdf['ID'].astype(str)
        gdf['EntityID'] = pd.factorize(gdf['UID'])[0] + 1  # Unique integer ID for each cell, starting from 1
        # dissolving to remove overlap artifacts
        gdf['gz_id'] = gdf['ID'].astype(str) + "_" + gdf['z'].astype(str)
        gdf = gdf.dissolve(by="gz_id").reset_index().drop(columns=['gz_id'])
    else: 
        gdf = gdf.dissolve(by="ID").reset_index()

    gdf.geometry = gdf.geometry.affine_transform(
        [scale, 0, 0, scale, 0, 0]
    )  # upscale polygons to original image size

    # Filter out small cells by size < min_size:
    if is_3d: # TODO: 3D case use volume?  
        logger.info(f"Filtering out cell in stack with area < {min_size}")
        n_pre = len(gdf)
        gdf = gdf[gdf.geometry.area > min_size]
        n_post = len(gdf)
        logger.info(f"Filtered out {n_pre - n_post} cells; remaining cells: {n_post}")

        logger.info(f"Filtering out cells that appear in fewer than {min_z} z layers")
        eid_vc = gdf.groupby("EntityID")['z'].count()
        keep_eids = eid_vc[eid_vc >= min_z].index
        n_pre = len(gdf)
        gdf = gdf[gdf["EntityID"].isin(keep_eids)].copy()
        n_post = len(gdf)
        logger.info(f"Filtered out {n_pre - n_post} cells; remaining cells: {n_post}")
    else: 
        logger.info(f"Filtering out cells with area < {min_size}")
        n_pre = len(gdf)
        gdf = gdf[gdf.geometry.area > min_size]
        n_post = len(gdf)
        logger.info(f"Filtered out {n_pre - n_post} cells; remaining cells: {n_post}")

    # Save masks and flows
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf.to_parquet(f"{output_dir}/{region}/polygons.parquet", index=False)

    return is_3d
