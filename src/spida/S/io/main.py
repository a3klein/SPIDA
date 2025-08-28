from spida._utilities import _gen_keys
from spida._constants import *
from .ingest_exp import (
    read_merscope,
    load_vpt_segmentation,
    load_proseg_segmentation_v2,
    load_proseg_segmentation_v3,
    load_decon_images,
)
from spida.pl import plot_images, plot_shapes, plot_points, plot_overlap, plot_seg_load

import os
import glob
from pathlib import Path
import warnings
import logging

from dotenv import load_dotenv  # type: ignore

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

from matplotlib.backends.backend_pdf import PdfPages

load_dotenv()
logger = logging.getLogger(__package__)


def ingest_region(
    exp_name: str,
    reg_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    root_path : str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Entry point for ingesting spatial data / segmentation outputs into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    type (str): Type of the data to ingest (default is "merscope").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    source (str): Source of the data (default is "machine").
    plot (bool): Whether to generate plots (default is False).
    root_path (str|Path): Root path for the data (default is None, uses environment variable).
    zarr_store (str|Path): Path to store the zarr files (default is None, uses environment variable).
    image_store (str|Path): Path to store the images (default is None, uses environment variable).
    kwargs: Additional keyword arguments for the ingestion functions.
    """

    logger.info("INGESTING REGION; EXPERIMENT %s, REGION %s " % (exp_name, reg_name))
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    if type == "merscope" and source == "machine":
        if root_path is None:
            root_path = os.getenv("PROCESSED_ROOT_PATH")
        if zarr_store is None:
            zarr_store = os.getenv("ZARR_STORAGE_PATH")
        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
        input_path = f"{root_path}/{exp_name}/out/{reg_name}"
        zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"

        sdata = read_merscope(
            input_path,
            zarr_path,
            exp_name=exp_name,
            reg_name=reg_name,
            prefix_name=prefix_name,
        )

        image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

        logger.info(sdata.tables.keys())
        if plot:
            image_path = f"{image_store}/{exp_name}/default/{reg_name}/pixi-ing.pdf"
            Path(image_path).parent.mkdir(parents=True, exist_ok=True)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                with PdfPages(image_path) as pdf:
                    plot_images(
                        sdata,
                        KEYS[IMAGE_KEY],
                        image_scale_keys,
                        image_channels,
                        cs="global",
                        pdf_file=pdf,
                    )
                    plot_shapes(
                        sdata,
                        KEYS[SHAPES_KEY],
                        table_name=KEYS[TABLE_KEY],
                        cs="pixel",
                        pdf_file=pdf,
                    )
                    plot_points(
                        sdata,
                        KEYS[POINTS_KEY],
                        KEYS[TABLE_KEY],
                        cs="pixel",
                        cmap="tab10",
                        pdf_file=pdf,
                    )
                    plot_overlap(
                        sdata,
                        KEYS[IMAGE_KEY],
                        KEYS[SHAPES_KEY],
                        KEYS[POINTS_KEY],
                        KEYS[TABLE_KEY],
                        image_scale_keys,
                        cs="pixel",
                        pdf_file=pdf,
                    )


def ingest_all(
    exp_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    root_path: str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Ingest all regions of an experiment into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    type (str): Type of the data to ingest (default is "merscope").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    source (str): Source of the data (default is "machine").
    plot (bool): Whether to generate plots (default is False).
    root_path (str|Path): Root path for the data (default is None, uses environment variable).
    zarr_store (str|Path): Path to store the zarr files (default is None, uses environment variable).
    image_store (str|Path): Path to store the images (default is None, uses environment variable).
    kwargs: Additional keyword arguments for the ingestion functions
    """

    logger.info("INGESTING ALL REGIONS; EXPERIMENT %s" % (exp_name))

    if type == "merscope" and source == "machine":
        if root_path is None:
            root_path = os.getenv("PROCESSED_ROOT_PATH")
        input_path = f"{root_path}/{exp_name}/out"
        regions = glob.glob(f"{input_path}/region_*")
        for reg in regions:
            ingest_region(
                exp_name,
                reg.split("/")[-1],
                type=type,
                prefix_name=prefix_name,
                source=source,
                plot=plot,
                root_path=root_path,
                zarr_store=zarr_store,
                image_store=image_store,
            )


def load_segmentation_region(
    exp_name: str,
    reg_name: str,
    seg_dir: str,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    zarr_store: str | Path | None = None,
    image_store: str | Path | None = None,
    **load_kwargs,
):
    """
    Load segmentation data into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    type (str): Type of the segmentation data to load (default is "vpt").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    """

    load_kwargs = load_kwargs.get("load_kwargs", load_kwargs)

    logger.info(
        "LOADING SEGMENTATION; EXPERIMENT %s, REGION %s, SEGMENTATION %s"
        % (exp_name, reg_name, type)
    )

    if zarr_store is None: 
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    logger.info(zarr_path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata = sd.read_zarr(zarr_path)

    if type == "vpt":
        sdata = load_vpt_segmentation(
            sdata,
            exp_name,
            reg_name,
            vpt_path=seg_dir,
            prefix_name=prefix_name,
            **load_kwargs,
        )
    elif type == "proseg":
        from spida.S.segmentation.proseg import _determine_proseg_version
        proseg_version = _determine_proseg_version()
        if proseg_version[0] == '2':
            sdata = load_proseg_segmentation_v2(
                sdata,
                exp_name,
                reg_name,
                proseg_path=seg_dir,
                prefix_name=prefix_name,
                **load_kwargs,
            )
        elif proseg_version[0] == '3':
            sdata = load_proseg_segmentation_v3(
                sdata,
                exp_name,
                reg_name,
                proseg_path=seg_dir,
                prefix_name=prefix_name,
                **load_kwargs,
            )

    # TODO: add a function to filter the spatialdata object to something very permissive (like at least 5 transcripts per cell?)

    if plot:
        # plot params
        KEYS = _gen_keys(prefix_name, exp_name, reg_name)
        # image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

        # define plot pdf
        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
        image_path = f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi-load.pdf"
        Path(image_path).parent.mkdir(parents=True, exist_ok=True)

        # plot
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with PdfPages(image_path) as pdf:
                # plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                # plot_shapes(
                #     sdata,
                #     KEYS[SHAPES_KEY],
                #     table_name=KEYS[TABLE_KEY],
                #     cs="pixel",
                #     pdf_file=pdf,
                # )
                plot_seg_load(
                    sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], cs="global", pdf_file=pdf
                )
                # plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                # plot_overlap(
                #     sdata,
                #     KEYS[IMAGE_KEY],
                #     KEYS[SHAPES_KEY],
                #     KEYS[POINTS_KEY],
                #     KEYS[TABLE_KEY],
                #     image_scale_keys,
                #     cs="pixel",
                #     pdf_file=pdf,
                # )


def load_segmentation_all(
    exp_name: str,
    seg_dir: str,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
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
        load_segmentation_region(
            exp_name,
            reg_name,
            reg_dir,
            type=type,
            prefix_name=prefix_name,
            plot=plot,
            **load_kwargs,
        )


def load_deconvolution_region(
    exp_name: str,
    reg_name: str,
    image_dir: str,
    image_name: str = "decon_image",
    suffix: str = ".decon.tif",
    plot: bool = False,
    root_path: str | Path | None = None,
    zarr_store: str | Path | None = None,
    image_store: str | Path | None = None,
    **load_kwargs,
):
    """
    Load deconvolution images into spatialdata objects.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    image_dir (str|Path): Directory containing the deconvolution images.
    image_name (str): Name of the image key in the spatialdata object (default is "decon_image").
    suffix (str): Suffix of the deconvolution image files (default is ".decon.tif").
    plot (bool): Whether to generate plots (default is False).
    root_path (str|Path): Root path for the data (default is None, uses environment variable).
    zarr_store (str|Path): Path to store the zarr files (default is
    image_store (str|Path): Path to store the images (default is None, uses environment variable).
    load_kwargs: Additional keyword arguments for the load_decon_images function.
    """

    logger.info(
        "LOADING DECONVOLUTION IMAGES; EXPERIMENT %s, REGION %s" % (exp_name, reg_name)
    )

    if image_dir is None:
        if root_path is None:
            image_dir = os.getenv("PROCESSED_ROOT_PATH")
        else: 
            image_dir = root_path
    image_dir = f"{image_dir}/{exp_name}/out/{reg_name}/images"
    logger.info(f"Loading deconvolution images from {image_dir}")

    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    logger.info(zarr_path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata = sd.read_zarr(zarr_path)

    sdata = load_decon_images(sdata, image_dir, image_name, suffix, **load_kwargs)

    if plot:
        image_channels = sd.models.get_channel_names(sdata[image_name])
        image_scale_keys = list(sdata[image_name].keys())

        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")

        image_path = f"{image_store}/{exp_name}/default/{reg_name}/pixi-decon.pdf"
        Path(image_path).parent.mkdir(parents=True, exist_ok=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with PdfPages(image_path) as pdf:
                plot_images(
                    sdata,
                    image_name,
                    image_scale_keys,
                    image_channels,
                    cs="global",
                    pdf_file=pdf,
                )
