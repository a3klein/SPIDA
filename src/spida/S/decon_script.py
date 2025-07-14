import os
import sys 

import argparse
import logging
from pathlib import Path
from tqdm import tqdm

import multiprocessing as mp
import warnings
from functools import partial

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import skimage as ski
import imageio.v3 as iio

from spida.utilities.tiling import tile_image_with_overlap, save_tiles, reconstruct_image_from_tiles
from spida.utilities.script_utils import parse_path, parse_list, parse_dict
from spida.utilities.read_raw import read_info
from spida.S.filters import deconwolf

logger = logging.getLogger("decon_large_image")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

def project_down_2D(input_file, output_file: str | Path = None): 
    """
    Project a 3D image to 2D by taking the maximum projection along the z-axis.
    """
    if output_file is None:
        output_file = input_file.with_suffix(".2d.tif")
    decon_img = iio.imread(input_file)

    if decon_img.ndim == 3 and decon_img.shape[0] > 1:
        # Remove first and last z slices if there are multiple z-slices
        if decon_img.shape[0] > 2:
            decon_img = decon_img[1:-1]
        # Z-project to 2D using max projection
        decon_img_2d = decon_img.max(axis=0)
    else:
        decon_img_2d = decon_img.squeeze()

    # Save the 2D projected tile
    iio.imwrite(output_file, decon_img_2d)

def continue_stalled_run(_file, **filter_args): 
    """
    Continue processing stalled runs by checking for existing files.
    """
    decon_file = _file.with_suffix(".decon.tif")
    projected_file = _file.with_suffix(".decon.2d.tif")
    
    if projected_file.exists(): 
        logger.info(f"{projected_file} previously processed, continuing")
        return projected_file
    
    if decon_file.exists(): 
        logger.info(f"Loading existing deconvolved file: {decon_file}")
        pass
    
    else:
        logger.info(f"Deconvolving {_file}")
        decon_file = deconwolf(_file, 
                                decon_path=decon_file,
                                **filter_args
                                )
        
    # Project to 2D and save immediately to avoid memory overhead
    project_down_2D(decon_file, projected_file)
    logger.info(f"Saved 2D projected tile to {projected_file}")        
    return projected_file

def fresh_run(_file, **filter_args):
    """ 
    Perform a fresh run of deconvolution and projection for a given file.
    """
    decon_file = _file.with_suffix(".decon.tif")
    projected_file = _file.with_suffix(".decon.2d.tif")

    decon_file = deconwolf(_file, 
                           decon_path=decon_file,
                           **filter_args
                           )
    # Project to 2D and save immediately to avoid memory overhead
    project_down_2D(decon_file, projected_file)
    logger.info(f"Saved 2D projected tile to {projected_file}")        
    return projected_file


def get_parser():
    """Get parser for decon_large_image"""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-i", "--image_path", type=parse_path, required=True, help="Path to the image file or directory")
    parser.add_argument("--data_org_path", type=str, required=True, help="Path to data organization file")
    parser.add_argument(
        "-o", "--output_dir", type=parse_path, default="tiles_output", help="Output directory for tiles"
    )
    parser.add_argument(
        "--channels", type=parse_list, required=True, help="Channel for segmentation (e.g., DAPI or PolyT,DAPI)"
    )
    parser.add_argument("-ts", "--tile_size", type=int, default=2960, help="Tile size in pixels (default: 2960)")
    parser.add_argument("--overlap", type=int, default=100, help="Overlap between tiles in pixels (default: 100)")
    parser.add_argument("--visualize_grid", action="store_true", help="Visualize the tiling grid")
    
    
    parser.add_argument("--z_step", type=float, default=1.5, help="axial(z) step size in micrometers")
    parser.add_argument("--filter", type=str, default=None, help="Filter to apply to the image before segmentation")
    parser.add_argument(
        "--filter_args", type=parse_dict, default={}, help="Additional filter arguments (e.g., key1=val1,key2=val2)"
    )
    parser.add_argument("--gpu", type=bool, default=False, help="Use GPU")
    parser.add_argument("--continue_stalled", type=bool, default=False, help="Continue processing if some tiles already processed")
    parser.add_argument("--plot_thr", type=bool, default=False, help="Plot thresholding histogram")
    parser.set_defaults(func=decon_large_image)
    return parser


def decon_large_image(
    image_path: str | Path,
    data_org_path: str | Path = "{input}/dataorganization.csv",
    channels: str | list[str] = "DAPI",
    tile_size: int | tuple = 2960,
    overlap: int | tuple = 100,
    output_dir: str | Path = "tiles_output",
    filter_args: dict = None,
    filter: str = "deconwolf",
    gpu: bool = True,
    z_step: float = 1.5,
    continue_stalled : bool = False,
    plot_thr : bool = False, 
    **kwargs
):
    """
    Deconvolve a large image by tiling it into smaller sections.
    """

    logger.info(f"Starting deconvolution of image: {image_path}")

    if isinstance(image_path, str):
        image_path = Path(image_path)

    if "{input}" in data_org_path:
        data_org_path = data_org_path.format(input=input)

    df_info = read_info(data_org_path)
    use_channels = channels if isinstance(channels, list) else [channels]
    channels = df_info.set_index("channel").loc[use_channels, :].reset_index()

    for _channel in channels['channel']: 
        print(_channel)
        # Load image
        channel_image_path = [fn for fn in image_path.glob("*.tif") if _channel in fn.name][0]
        large_img = iio.imread(channel_image_path)
        logger.info(f"Loaded Image with shape: {large_img.shape} from {channel_image_path}")

        # Tile Images
        tiles, tile_info = tile_image_with_overlap(
            large_img,
            tile_size=tile_size,
            overlap=overlap,
            color=_channel,
        )
        logger.info(f"Tiled image into {len(tiles)} tiles with tile_size: {tile_size}, overlap: {overlap}")

        ### Filtering the Tiles to remove empty tiles 
        tile_maxes = [tile['max_intensity'] for tile in tile_info]
        thr = ski.filters.threshold_minimum(np.asarray(tile_maxes))
        for _ti in tile_info: 
            _ti['thresholded'] = _ti['max_intensity'] > thr
        logging.info(f"Applied Minimum threshold: {thr} to tile max intensities")
        if plot_thr: 
            fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
            sns.histplot(tile_maxes, bins=100, kde=True, ax=ax)
            ax.axvline(thr, c="red", linestyle="--", label="Threshold")
            fig.savefig(output_dir / f"tile_maxes_{_channel}.png")
        ret = [(t, ti) for t, ti in zip(tiles, tile_info) if ti['thresholded']]
        tiles = [t for t, _ in ret]
        tile_info = [ti for _, ti in ret]
        logger.info(f"Filtered tiles to {len(tiles)} with max intensity above threshold: {thr}")

        # Save Tiles (applying 7 slice expansion): 
        to_3D = lambda x: np.asarray([x]*7)
        saved_files = save_tiles(tiles, tile_info, output_dir, func=to_3D)
        logger.info(f"Saved {len(saved_files)} tiles to {output_dir}")

        del tiles # Clear memory

        logger.info(f"Applying {filter} filter from fishtank")
        # Setting the deconwolf args
        filter_args = filter_args if filter_args else {}
        filter_args['colors'] = channels['color'].values
        filter_args['z_step'] = z_step
        filter_args['iter'] = 100
        filter_args['gpu'] = gpu

        l = mp.Lock() # for locking the deconwolf subprocesses
        
        if continue_stalled:
            logger.info("Continuing stalled run...")
            
            with mp.Manager() as manager: 
                lock = manager.Lock()
                filter_args['lock'] = lock
                # Deconvolve each tile
                parallel_func = partial(
                    continue_stalled_run,
                    **filter_args
                )
                with mp.Pool(mp.cpu_count()) as pool:
                    deconed_paths = list(tqdm(pool.imap_unordered(parallel_func, saved_files), total=len(saved_files)))

            # continue_stalled_run(saved_files, **filter_args)
        else:
            logger.info("Starting new deconvolution run...")

            with mp.Manager() as manager: 
                lock = manager.Lock()
                filter_args['lock'] = lock
            
                # Deconvolve each tile
                parallel_func = partial(
                    fresh_run,
                    **filter_args
                )
                # with mp.Pool(mp.cpu_count(), initializer=init, initargs=(l, )) as pool:
                with mp.Pool(mp.cpu_count()) as pool:
                    deconed_paths = list(tqdm(pool.imap_unordered(parallel_func, saved_files), total=len(saved_files)))
    
        logger.info(f"Completed deconvolution and 2D projection for all tiles.")

        # Reconstruct deconvolved image from saved 2D tiles
        deconed_image = reconstruct_image_from_tiles(
            output_dir=output_dir, 
            tile_info=tile_info, 
            original_shape=large_img.shape,
            suffix=".decon.2d"
        )
        logger.info(f"Reconstructed deconvolved image with shape: {deconed_image.shape}")

        iio.imsave(channel_image_path.with_suffix(".decon.tif"), deconed_image)
        logger.info(f"Saved deconvolved image to {channel_image_path.with_suffix('.decon.tif')}")

        return deconed_image
