import os
from dotenv import load_dotenv # type: ignore
load_dotenv()

import numpy as np
from pathlib import Path
import warnings

from natsort import natsorted
import cv2

from cellpose import models, core, io # type: ignore

io.logger_setup

if core.use_gpu()==False: 
    raise ImportError("No GPU access")


import skimage as ski
import supervision as spv
from shapely.geometry import Polygon
import geopandas as gpd


def _load_image(image_path:Path, image_ext:str='.tif', nuc_stain_name:str='DAPI', cyto_stain_name:str ='PolyT'): 
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
    files = natsorted([f for f in image_path.glob("*"+image_ext) if "_masks" not in f.name and "_flows" not in f.name])
    print(f"Found {len(files)} images in {image_path}")
    
    nuc_file = None
    for f in files: 
        if nuc_stain_name in f.name: 
            nuc_file = f
    
    if nuc_file is None: 
        raise ValueError(f"No nuclear stain file found with name '{nuc_stain_name}' in {image_path}")

    nuc_img = io.imread(nuc_file)
    print(f"Loaded images: Nuclear - {nuc_file}, Image Shape: {nuc_img.shape}")

    if cyto_stain_name is not None: 
        for f in files: 
            if cyto_stain_name in f.name: 
                cyto_file = f
        
        cyto_img = io.imread(cyto_file)
        print(f"Loaded images: Cytoplasmic - {cyto_file}, Image Shape: {cyto_img.shape}")
        img = np.stack([nuc_img, cyto_img], axis=-1)
        return img

    else: 
        print(f"No cytoplasmic stain file found with name '{cyto_stain_name}' in {image_path}, returning nuclear image only.")
        return np.expand_dims(nuc_img, axis=-1)

def _load_model(model_name: str = "cpsam", **kwargs): 
    """
    Load The Cellpose Model 
    """
    model_path = os.getenv("CELLPOSE_MODEL_PATH", "/anvil/projects/x-mcb130189/aklein/programs/.cellpose/models/") # Set the GPU to use
    model=models.CellposeModel(pretrained_model = f"{model_path}/{model_name}", 
                               gpu=True)
    return model

def _cellpose_wrapper(img, 
                      batch_size=8, 
                      flow_threshold=0.0, 
                      cellprob_threshold=-4, 
                      tile_norm_blocksize=0,
                      **kwargs):
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

    return model.eval(img,
                      batch_size=batch_size,
                      flow_threshold=flow_threshold, 
                      cellprob_threshold=cellprob_threshold, 
                      normalize={"tile_norm_blocksize": tile_norm_blocksize},)


def _masks_to_geodataframe(masks):
    """
    Convert Cellpose masks to a GeoDataFrame with polygons.
    Parameters:
    masks (np.ndarray): Cellpose segmentation masks.
    Returns:
    gpd.GeoDataFrame: GeoDataFrame containing polygons of the masks.
    """

    polygons = spv.mask_to_polygons(masks)
    sh_polygons = [Polygon(p) for p in polygons]

    gdf = gpd.GeoDataFrame(
        {
            "Geometry": sh_polygons,
            "ID": [i for i in range(len(sh_polygons))],
        },
        geometry="Geometry",
    )
    return gdf

def run_cellposeSAM(
        root_dir:str,
        output_dir:str,
        region:str,
        scale:int = 4, # value to downsample image by
        min_size:int = 100, # minimum size of cell to keep
        image_ext:str = ".tif",
        nuc_stain_name:str = "DAPI",
        cyto_stain_name:str = "PolyT",
        **kwargs):
    
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

    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    data_path =f"{root_dir}/{region}/images/"
    
    img = _load_image(Path(data_path),
                      image_ext=image_ext,
                      nuc_stain_name=nuc_stain_name,
                      cyto_stain_name=cyto_stain_name,
                      )
    
    og_shape = img.shape
    print("OG image shape:", og_shape)

    downscale = (scale, scale) + (1,) * (img.ndim - 2)
    img = ski.transform.downscale_local_mean(img, downscale)
    
    # img = img[::scale, ::scale, :]
    print("Image shape after scaling:", img.shape)
    
    print("STARTING SEGMENTATION")
    masks, flows, styles = _cellpose_wrapper(img, batch_size=16, **kwargs)
    print("SEGMENTATION COMPLETED")

    print(f"Cellpose segmentation completed for region {region}.")
    print(f"Number of masks detected: {len(np.unique(masks)) - 1}; mask shape: {masks.shape}")

    # rescaled_masks = cv2.resize(masks, (og_shape[1], og_shape[0]), interpolation=cv2.INTER_NEAREST)
    # print("Rescaled masks shape:", rescaled_masks.shape)

    # masks to geodataframe
    gdf = _masks_to_geodataframe(masks)
    gdf.geometry = gdf.geometry.affine_transform([scale, 0, 0, scale, 0, 0]) # upscale polygons to original image size

    # Filter out small cells by size < min_size: 
    print(f"Filtering out cells with area < {min_size}")
    n_pre = len(gdf)
    gdf = gdf[gdf.geometry.area > min_size]
    n_post = len(gdf)
    print(f"Filtered out {n_pre - n_post} cells; remaining cells: {n_post}")
    
    # Save masks and flows
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf.to_parquet(f"{output_dir}/{region}/polygons.parquet", index=False)
