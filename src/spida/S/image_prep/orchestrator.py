
import logging
from pathlib import Path
from tqdm import tqdm

import multiprocessing as mp
from dataclasses import replace
from functools import partial

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import skimage as ski
import imageio.v3 as iio
import tifffile

from spida.io.schema import tissue_mask

# TODO: Implement some image normalization pre and post tiling to normalize the images (right now specific
# patches of the reuslting images are very dark, and others very bright, which doesn't affect cellpose becuase it looks
# at local contrast, but it does affect visualization and it might affect other segmentation methods.)

from spida.utilities.tiling import (
    tile_image_with_overlap,
    save_tiles,
    reconstruct_image_from_tiles,
    reconstruct_image_from_tile_files,
    visualize_tiling_grid,
)

# from spida.utilities.script_utils import parse_path, parse_list, parse_dict
from spida.utilities.read_raw import read_info
from spida.S.image_prep.tile_io import project_down_2D, project_to_3D, write_tiff
from spida.S.image_prep.tile_selection import select_tiles
from spida.S.image_prep.transforms import TileContext, resolve_transform

logger = logging.getLogger(__package__)

from spida.pl.decon_plots import save_threshold_plots_pdf


# unknown output location fails at startup.
def _stitch_dir(channel_image_path: Path, stitched_dir: str | Path | None) -> Path:
    """Where stitched output goes: alongside the input by default, or an explicit directory."""
    if stitched_dir is None:
        return channel_image_path.parent
    d = Path(stitched_dir)
    d.mkdir(parents=True, exist_ok=True)
    return d

def continue_stalled_run(_file, _transform_name=None, _ctx=None, **transform_args):
    """
    Continue processing stalled runs by checking for existing files.
    """
    _transform = resolve_transform(_transform_name)
    decon_file = _file.with_suffix(".decon.tif")
    if is_3d and _zstackexpand:
        projected_file = _file.with_suffix(".decon.3d.tif")
    else:
        projected_file = _file.with_suffix(".decon.2d.tif")

    if projected_file.exists():
        logger.info(f"{projected_file} previously processed, continuing")
        return projected_file

    if decon_file.exists():
        logger.info(f"Loading existing deconvolved file: {decon_file}")
        pass

    else:
        logger.info(f"Transforming {_file}")
        decon_file = _transform.run(_file, decon_file, _ctx, **transform_args)

    # Project to 2D and save immediately to avoid memory overhead
    if not is_3d: 
        project_down_2D(decon_file, projected_file)
        logger.info(f"Saved 2D projected tile to {projected_file}")

    if is_3d and _zstackexpand:
        project_to_3D(decon_file, projected_file, z_stack_size=z_stack_size)
        decon_file.unlink()  # Padded decon tile no longer needed — free disk space immediately.
        logger.info(f"Saved 3D projected tile to {projected_file}")

    return projected_file


def fresh_run(_file, _transform_name=None, _ctx=None, **transform_args):
    """
    Perform a fresh run of the transform and projection for a given file.
    """
    _transform = resolve_transform(_transform_name)
    decon_file = _file.with_suffix(".decon.tif")

    decon_file = _transform.run(_file, decon_file, _ctx, **transform_args)

    if not is_3d:
        # Project to 2D and save immediately to avoid memory overhead
        projected_file = _file.with_suffix(".decon.2d.tif")
        project_down_2D(decon_file, projected_file)
        logger.info(f"Saved 2D projected tile to {projected_file}")

    if is_3d and _zstackexpand:
        projected_file = _file.with_suffix(".decon.3d.tif")
        project_to_3D(decon_file, projected_file, z_stack_size=z_stack_size)
        decon_file.unlink()  # Padded decon tile no longer needed — free disk space immediately.
        logger.info(f"Saved 3D projected tile to {projected_file}")

    return projected_file


def transform_image(
    image_path: str | Path,
    data_org_path: str | Path = "{input}/dataorganization.csv",
    channels: str | list[str] = "DAPI",
    tile_size: int | tuple = 2960,
    overlap: int | tuple = 100,
    output_dir: str | Path = "tiles_output",
    transform_args: dict = None,
    transform: str = "deconwolf",
    gpu: bool = True,
    z_step: float = 1.5,
    z_stack_expand: bool = True,
    z_stack_expand_slices: int = 3,
    continue_stalled: bool = False,
    thr_tiles: bool = True,
    plot_thr: bool = False,
    match_pre: bool = False,
    image_store: str | Path | None = None,
    clear_tiles : bool = True,
    select_tiles_method: str | None = None,
    mask_path: str | Path | None = None,
    stitched_dir: str | Path | None = None,
    deconwolf_config: str | Path | None = None,
    **kwargs,
):
    """
    Apply a per-tile transform to a large image by tiling it into smaller sections.

    stitched_dir
        Where the stitched `*_z{n}.decon.tif` outputs go. Default None next to image_path
    deconwolf_config
        The deconwolf .config.ini (binary, PSF dir, [Optics]); required by deconwolf.
    """

    logger.info(f"Starting {transform} on image: {image_path}")

    _transform_fn = resolve_transform(transform)

    if image_store is None:
        image_store = Path(output_dir)
    exp_name = Path(data_org_path).parents[1].name
    reg_name = Path(image_path).parents[0].name
    image_save_dir = Path(image_store) / exp_name / "default" / reg_name
    image_save_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(image_path, str):
        image_path = Path(image_path)

    if "{input}" in data_org_path:
        inp = image_path.parents[2] / 'raw'
        data_org_path = data_org_path.format(input=inp)

    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df_info = read_info(data_org_path)
    use_channels = channels if isinstance(channels, list) else [channels]
    channels = df_info.set_index("channel").loc[use_channels, :].reset_index()

    for _channel in channels["channel"]:
        logger.info(f"Loading channel: {_channel}")
        channel_image_path = [
            fn for fn in image_path.glob("*.tif") if (_channel in fn.name) & ("decon" not in fn.name)
        ]
        global is_3d
        is_3d = len(channel_image_path) > 1

        if is_3d:
            global z_stack_size
            global _zstackexpand
            z_stack_size = len(channel_image_path)
            _zstackexpand = z_stack_expand
            channel_image_path = sorted(channel_image_path)
            # Pre-allocate the full 3D volume and read each z-slice directly into it.
            # This avoids the 2× peak that np.stack creates (list-of-slices + new array simultaneously).
            logger.info(f"..{channel_image_path[0]}")
            first_slice = iio.imread(channel_image_path[0])
            large_img = np.empty((z_stack_size, *first_slice.shape), dtype=first_slice.dtype)
            large_img[0] = first_slice
            del first_slice
            for i, _fn in enumerate(channel_image_path[1:], 1):
                logger.info(f"..{_fn}")
                large_img[i] = iio.imread(_fn)
            channel_image_path = channel_image_path[0].parent / f"{channel_image_path[0].stem[:-1]}_stack.tif"
        else:
            channel_image_path = channel_image_path[0]
            large_img = iio.imread(channel_image_path)
        logger.info(
            f"Loaded Image with shape: {large_img.shape} from {channel_image_path}"
        )

        # Tile Images
        tiles, tile_info = tile_image_with_overlap(
            large_img,
            tile_size=tile_size,
            overlap=overlap,
            color=_channel,
        )
        logger.info(
            f"Tiled image into {len(tiles)} tiles with tile_size: {tile_size}, overlap: {overlap}"
        )

        ### Drop empty tiles (strategy dispatch; see image_prep.tile_selection)
        # select_tiles_method overrides thr_tiles (backwards compatible: True->otsu, False->none)
        method = select_tiles_method or ("otsu" if thr_tiles else "none")
        select_kw = {}
        if method == "tissue_mask":
            mask_tif = Path(mask_path) if mask_path else Path(image_path) / tissue_mask.FILENAME
            select_kw["mask_tif"] = mask_tif
            logger.info(f"tissue_mask tile filter using mask: {mask_tif}")
        tiles, sub_tile_info, tile_maxes, thr = select_tiles(
            tiles, tile_info, method=method, **select_kw,
        )

        # plotting threshold decision
        if plot_thr:
            logger.info("Plotting thresholding histogram and tile grid visualizations")
            save_threshold_plots_pdf(
                image_save_dir=image_save_dir,
                channel=_channel,
                grid_tile_info=tile_info,
                recon_tile_info=sub_tile_info,
                tiles=tiles,
                large_img=large_img,
                tile_maxes=tile_maxes,
                thr=thr,
                is_3d=is_3d,
            )
        
        # save tiles 
        if is_3d:
            logger.info("Saving 3D tiles")
            # Adding padding on the edges of the z-stack for deconwolf
            def expand_z_stack(x):
                exp_x = []
                for z in range(len(x)): 
                    if (z == 0) or (z == len(x) - 1): 
                        for _ in range(z_stack_expand_slices):
                            exp_x.append(x[0])
                    else: 
                        exp_x.append(x[z])
                return np.asarray(exp_x)
            # expand_z_stack: prepend z_stack_expand_slices copies of the first slice and
            # append z_stack_expand_slices copies of the last slice, then include ALL original
            # slices in between.  For 7 slices + expand=3 the output is 13 slices:
            #   [z0,z0,z0, z0,z1,z2,z3,z4,z5,z6, z6,z6,z6]
            # project_to_3D then extracts the middle 7: start=(13-7)//2=3 → [z0..z6].
            # Previous code was buggy: it produced 11 slices, dropped z0 and z6 as data,
            # and used x[0] (not x[-1]) for end padding, corrupting the last two output slices.
            def expand_z_stack(x):
                front = [x[0]] * z_stack_expand_slices
                data  = [x[z] for z in range(len(x))]
                back  = [x[-1]] * z_stack_expand_slices
                return np.stack(front + data + back)
            if z_stack_expand:
                saved_files = save_tiles(tiles, sub_tile_info, output_dir, func=expand_z_stack)
            else:
                saved_files = save_tiles(tiles, sub_tile_info, output_dir)
        else:
            logger.info("Saving 2D tiles")
            # applying 7 slice expansion
            def to_3D(x):
                return np.asarray([x] * 7)
            saved_files = save_tiles(tiles, sub_tile_info, output_dir, func=to_3D)
        logger.info(f"Saved {len(saved_files)} tiles to {output_dir}")

        del tiles  # Tiles are NumPy views into large_img, not copies; this only removes the view objects.
        # Capture shape metadata then release the 3D volume from RAM — all tile data is on disk
        # and large_img is not needed again until reconstruction (where only the shape is required).
        _img_shape = large_img.shape
        if is_3d:
            del large_img  # Release the full 3D volume (z_stack_size × FOV) during deconvolution.

        logger.info(f"Applying {transform} transform")
        # Defaults and validation belong to the transform, not here: it is the only
        # thing that knows which parameters it accepts and how they interact.
        _params = {**_transform_fn.DEFAULTS, **(transform_args or {})}
        _params = _transform_fn.validate(_params)

        _ctx = TileContext(
            channel=_channel,
            colors=channels["color"].values,
            z_step_um=z_step,
            gpu=gpu,
            optics_config=deconwolf_config,
        )

        if continue_stalled:
            logger.info("Continuing stalled run...")

            with mp.Manager() as manager:
                # Only GPU transforms take the lock; CPU ones run the pool in parallel.
                ctx = (replace(_ctx, gpu_lock=manager.Lock())
                       if _transform_fn.USES_GPU else _ctx)
                # Transform each tile
                parallel_func = partial(continue_stalled_run,
                                        _transform_name=transform, _ctx=ctx, **_params)
                with mp.Pool(mp.cpu_count()) as pool:
                    deconed_paths = list(  # noqa: F841
                        tqdm(
                            pool.imap_unordered(parallel_func, saved_files),
                            total=len(saved_files),
                        )
                    )

        else:
            logger.info("Starting new run...")

            with mp.Manager() as manager:
                ctx = (replace(_ctx, gpu_lock=manager.Lock())
                       if _transform_fn.USES_GPU else _ctx)

                # Transform each tile
                parallel_func = partial(fresh_run,
                                        _transform_name=transform, _ctx=ctx, **_params)
                # with mp.Pool(mp.cpu_count(), initializer=init, initargs=(l, )) as pool:
                with mp.Pool(mp.cpu_count()) as pool:
                    deconed_paths = list( # noqa: F841
                        tqdm(
                            pool.imap_unordered(parallel_func, saved_files),
                            total=len(saved_files),
                        )
                    ) 

        logger.info("Completed deconvolution and 2D projection for all tiles.")

        if is_3d:
            depth = _img_shape[0]
            suffix = ".decon.3d" if _zstackexpand else ".decon"
            logger.info(
                f"Reconstructing deconvolved 3D image one z-slice at a time (depth={depth})"
            )
            written = []
            for z in tqdm(range(depth), desc=f"Stitching z-slices ({_channel})"):
                z_slice = reconstruct_image_from_tile_files(
                    output_dir=output_dir,
                    tile_info=sub_tile_info,
                    original_shape=_img_shape,
                    suffix=suffix,
                    match_pre=match_pre,
                    z_index=z,
                )
                out_path = (
                    _stitch_dir(channel_image_path, stitched_dir) / f"{channel_image_path.stem}_z{z:01d}.decon.tif"
                )
                write_tiff(out_path, z_slice)
                written.append(out_path)
                del z_slice

            logger.info(
                f"Saved {len(written)} stitched z-slices under {_stitch_dir(channel_image_path, stitched_dir)}"
            )
            
            # Optionally clear tile files after stitching to free up space
            if clear_tiles: 
                logger.info("Clearing tile files to free up space...")
                for f in output_dir.glob(f"*{_channel}*"):
                    f.unlink()
                logger.info("Tile files cleared.")
            
            return written

        # 2D path (single reconstructed image)
        deconed_image = reconstruct_image_from_tile_files(
            output_dir=output_dir,
            tile_info=sub_tile_info,
            original_shape=_img_shape,
            suffix=".decon.2d",
            match_pre=match_pre,
        )
        logger.info(f"Reconstructed deconvolved image with shape: {deconed_image.shape}")

        _out_2d = (
            _stitch_dir(channel_image_path, stitched_dir) / f"{channel_image_path.stem}.decon.tif"
        )
        with tifffile.TiffWriter(_out_2d, bigtiff=True) as tiff_writer:
            tiff_writer.write(deconed_image)
        logger.info(f"Saved deconvolved image to {_out_2d}")

        # Optionally clear tile files after stitching to free up space
        if clear_tiles: 
            logger.info("Clearing tile files to free up space...")
            for f in output_dir.glob(f"*{_channel}*"):
                f.unlink()
            logger.info("Tile files cleared.")

        return deconed_image
