from pathlib import Path
from natsort import natsorted
from dotenv import load_dotenv # type: ignore
import logging

from PIL import Image
import numpy as np

import supervision as spv
from shapely.geometry import Polygon
import geopandas as gpd

from deepcell.applications import Mesmer # type: ignore

Image.MAX_IMAGE_PIXELS = None  # Disable the limit on image size
load_dotenv()
logger = logging.getLogger(__package__)


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
    logger.info(f"Found {len(files)} images in {image_path}")
    
    for f in files: 
        if nuc_stain_name in f.name: 
            nuc_file = f
        if cyto_stain_name in f.name: 
            cyto_file = f

    nuc_img = Image.open(nuc_file)
    nuc_img = np.array(nuc_img) 

    cyto_img = Image.open(cyto_file)
    cyto_img = np.array(cyto_img) 

    logger.info(f"Loaded images: Nuclear - {nuc_file}, Cytoplasmic - {cyto_file}")
    logger.info(f"Image shape: Nuclear - {nuc_img.shape}, Cyto - {cyto_img.shape}")

    img = np.stack([nuc_img, cyto_img], axis=-1)
    return np.expand_dims(img, axis=0)

def _mesmer_wrapper(img:np.array, model, batch_size=16, compartment='nuclear', image_mpp=0.2): 
    """
    Wrapper for Mesmer model evaluation.
    Parameters:
    img (np.ndarray): Input image.
    model (Mesmer): Pre-trained Mesmer model.
    batch_size (int): Batch size for Mesmer evaluation.
    compartment (str): Compartment to segment ('nuclear' or 'cytoplasmic').
    image_mpp (float): Microns per pixel for the image.
    Returns:
    np.ndarray: Predicted segmentation masks.
    """

    return model.predict(img, batch_size=batch_size, compartment=compartment, image_mpp=image_mpp)

# TODO: fix this like I fixed the Cellpose SAM one
def _masks_to_geodataframe(masks):
    """
    Convert segmentation masks to a GeoDataFrame with polygons.
    Parameters:
    masks (np.ndarray): segmentation masks.
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


def run_mesmer(root_dir:Path, output_dir:Path, region:str, **kwargs): 
    """
    Run Mesmer segmentation on a specified region.
    Parameters:
    root_dir (str): The root directory containing the images.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    kwargs (dict): Additional keyword arguments for Mesmer configuration.
    """

    # Ensure the output directory exists
    Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    # Load the image
    data_path =f"{root_dir}/{region}/images/"
    img = _load_image(Path(data_path),
                      image_ext=kwargs.get("image_ext", ".tif"),
                      nuc_stain_name=kwargs.get("nuc_stain_name", "DAPI"),
                      cyto_stain_name=kwargs.get("cyto_stain_name", "PolyT"))
    
    app = Mesmer()
    logger.info("STARTING SEGMENTATION")
    masks = _mesmer_wrapper(img, app)
    logger.info("SEGMENTATION COMPLETED")

    logger.info(f"Mesmer segmentation completed for region {region}.")
    logger.info(f"Number of masks detected: {len(np.unique(masks)) - 1}")

    # Save masks and flows
    gdf = _masks_to_geodataframe(masks[0, ..., 0])
    gdf.to_parquet(f"{output_dir}/{region}/polygons.parquet", index=False)