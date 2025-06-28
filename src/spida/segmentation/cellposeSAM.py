
import numpy as np
from pathlib import Path

from natsort import natsorted
import cv2

from cellpose import models, core, io # type: ignore

io.logger_setup

if core.use_gpu()==False: 
    raise ImportError("No GPU access")

model=models.CellposeModel(gpu=True)

import supervision as spv
from shapely.geometry import Polygon
import geopandas as gpd


def _load_image(image_path:Path, image_ext:str='.tif', nuc_stain_name:str='DAPI', cyto_stain_name:str='PolyT'): 
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
    
    for f in files: 
        if nuc_stain_name in f.name: 
            nuc_file = f
        if cyto_stain_name in f.name: 
            cyto_file = f

    nuc_img = io.imread(nuc_file)
    cyto_img = io.imread(cyto_file)
    print(f"Loaded images: Nuclear - {nuc_file}, Cytoplasmic - {cyto_file}")
    print(f"Image shape: Nuclear - {nuc_img.shape}, Cyto - {cyto_img.shape}")

    img = np.stack([nuc_img, cyto_img], axis=-1)
    return img


def _cellpose_wrapper(img, 
                      batch_size=8, 
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
    return model.eval(img,
                      batch_size=batch_size,
                      flow_threshold=kwargs.get("flow_threshold", 0.6), 
                      cellprob_threshold=kwargs.get("cellprob_threshold", 0.0), 
                      normalize={"tile_norm_blocksize":kwargs.get("tile_norm_blocksize", 0)})


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

def run_cellposeSAM(root_dir:str, output_dir:str, region:str, **kwargs):
    """
    Run Cellpose SAM segmentation on a specified region.
    Parameters:
    root_dir (str): The root directory containing the images.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    """
    # scale for scaling image to allow for larger images in memory
    scale = kwargs.get("scale", 2)
    
    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    data_path =f"{root_dir}/{region}/images/"
    img = _load_image(Path(data_path),
                      image_ext=kwargs.get("image_ext", ".tif"),
                      nuc_stain_name=kwargs.get("nuc_stain_name", "DAPI"),
                      cyto_stain_name=kwargs.get("cyto_stain_name", "PolyT"),
                      )
    
    og_shape = img.shape
    img = img[::scale, ::scale, :]
    print("Image shape after scaling:", img.shape)
    
    print("STARTING SEGMENTATION")
    masks, flows, styles = _cellpose_wrapper(img, batch_size=16, kwargs=kwargs)
    print("SEGMENTATION COMPLETED")

    print(f"Cellpose segmentation completed for region {region}.")
    print(f"Number of masks detected: {len(np.unique(masks)) - 1}; mask shape: {masks.shape}")

    rescaled_masks = cv2.resize(masks, (og_shape[1], og_shape[0]), interpolation=cv2.INTER_NEAREST)
    print("Rescaled masks shape:", rescaled_masks.shape)
    
    # Save masks and flows
    gdf = _masks_to_geodataframe(rescaled_masks)
    gdf.to_parquet(f"{output_dir}/{region}/polygons.parquet", index=False)
