import os
import sys
import glob
import pathlib
import warnings
import logging

from dotenv import load_dotenv # type: ignore
load_dotenv()

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

from spida._utilities import _gen_keys
from spida._constants import *
from .ingest_exp import read_merscope, load_vpt_segmentation, load_proseg_segmentation
from spida.pl import plot_images, plot_shapes, plot_points, plot_overlap

from matplotlib.backends.backend_pdf import PdfPages


logger = logging.getLogger(__package__)


def ingest_region(exp_name:str, 
                  reg_name:str,
                  type:str="merscope",  
                  prefix_name:str="default", 
                  source:str="machine", 
                  plot:bool=False, 
                  **kwargs):
    """
    Entry point for ingesting spatial data / segmentation outputs into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    type (str): Type of the data to ingest (default is "merscope").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").

    """

    logger.info("INGESTING REGION; EXPERIMENT %s, REGION %s " %(exp_name, reg_name) )
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    if type == "merscope" and source == "machine":
        root_path = os.getenv("PROCESSED_ROOT_PATH", "/ceph/cephatlas/merscope_data/processed")
        processed_path = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
        input_path = f"{root_path}/{exp_name}/out/{reg_name}"
        zarr_path = f"{processed_path}/{exp_name}/{reg_name}"

        sdata = read_merscope(input_path, zarr_path, exp_name=exp_name, reg_name=reg_name, prefix_name=prefix_name)

        image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

        logging.info(sdata.tables.keys())
        if plot: 
            image_path = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
            image_path = f"{image_path}/{exp_name}/default/{reg_name}/pixi-ing.pdf"
            pathlib.Path(image_path).parent.mkdir(parents=True, exist_ok=True)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                with PdfPages(image_path) as pdf:
                    plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                    plot_shapes(sdata, KEYS[SHAPES_KEY], table_name=KEYS[TABLE_KEY], cs="pixel", pdf_file=pdf)
                    plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                    plot_overlap(sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[POINTS_KEY], KEYS[TABLE_KEY], image_scale_keys, cs="pixel", pdf_file=pdf)


def ingest_all(
        exp_name:str, 
        type:str="merscope", 
        prefix_name:str="default",
        source:str="machine",
        plot:bool=False, 
        **kwargs):
    """
    Ingest all regions of an experiment into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    type (str): Type of the data to ingest (default is "merscope").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    source (str): Source of the data (default is "machine").
    """
    
    logger.info("INGESTING ALL REGIONS; EXPERIMENT %s" %(exp_name) )

    if type == "merscope" and source == "machine":
        root_path = os.getenv("PROCESSED_ROOT_PATH", "/ceph/cephatlas/merscope_data/processed")
        input_path = f"{root_path}/{exp_name}/out"
        regions = glob.glob(f"{input_path}/region_*")
        for reg in regions: 
            ingest_region(exp_name, reg.split("/")[-1], type=type, prefix_name=prefix_name, source=source, plot=plot)



def load_segmentation_region(
        exp_name:str, 
        reg_name:str, 
        seg_dir:str,
        type:str="vpt", 
        prefix_name:str="default", 
        plot:bool=False,
        **load_kwargs
        ):
    """
    Load segmentation data into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    type (str): Type of the segmentation data to load (default is "vpt").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    """

    logger.info("LOADING SEGMENTATION; EXPERIMENT %s, REGION %s, SEGMENTATION %s" %(exp_name, reg_name, type) )


    zarr_root = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_root}/{exp_name}/{reg_name}"
    logger.info(zarr_path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata = sd.read_zarr(zarr_path)

    if type == "vpt":
        sdata = load_vpt_segmentation(sdata, exp_name, reg_name, vpt_path=seg_dir, prefix_name=prefix_name, **load_kwargs)
    elif type == "proseg":
        sdata = load_proseg_segmentation(sdata, exp_name, reg_name, proseg_path=seg_dir, prefix_name=prefix_name, **load_kwargs)

    if plot: 
        # plot params
        KEYS = _gen_keys(prefix_name, exp_name, reg_name)
        image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())
        
        # define plot pdf
        image_root = os.getenv("IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images")
        image_path = f"{image_root}/{exp_name}/{prefix_name}/{reg_name}/pixi-load.pdf"
        pathlib.Path(image_path).parent.mkdir(parents=True, exist_ok=True)
        
        # plot
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with PdfPages(image_path) as pdf:
                plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                plot_shapes(sdata, KEYS[SHAPES_KEY], table_name=KEYS[TABLE_KEY], cs="pixel", pdf_file=pdf)
                plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                plot_overlap(sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[POINTS_KEY], KEYS[TABLE_KEY], image_scale_keys, cs="pixel", pdf_file=pdf)


def load_segmentation_all(
        exp_name:str,
        seg_dir:str,
        type:str="vpt", 
        prefix_name:str="default",
        plot:bool=False,
        **load_kwargs,
        ):
    """
    Load segmentation data for all regions of an experiment into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    seg_dir (str): Directory containing the segmentation data.
    type (str): Type of the segmentation data to load (default is "vpt").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    """

    regions = glob.glob(f"{seg_dir}/region_*")
    for reg_dir in regions: 
        reg_name = reg_dir.split("/")[-1]
        load_segmentation_region(exp_name, reg_name, reg_dir, type=type, prefix_name=prefix_name, plot=plot, **load_kwargs)