from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from spida.utilities.tiling import visualize_tiling_grid


def choose_downscale(image_shape: tuple[int, ...], target_pixels: int = 2000) -> int:
    max_dim = max(image_shape[-2:])
    return max(1, int(math.ceil(max_dim / target_pixels)))


def downsample_original_slice(
    image: np.ndarray,
    downscale: int,
    z_index: int | None = None,
) -> np.ndarray:
    if image.ndim == 3:
        if z_index is None:
            z_index = image.shape[0] // 2
        image = image[z_index]
    return image[::downscale, ::downscale].astype(np.float32, copy=False)


def reconstruct_downsampled_from_tiles(
    tiles: list[np.ndarray],
    tile_info: list[dict],
    original_shape: tuple[int, ...],
    downscale: int,
    z_index: int | None = None,
) -> np.ndarray:
    height, width = original_shape[-2:]
    height_ds = (height + downscale - 1) // downscale
    width_ds = (width + downscale - 1) // downscale

    reconstructed = np.zeros((height_ds, width_ds), dtype=np.float32)
    weight_map = np.zeros((height_ds, width_ds), dtype=np.float32)

    for tile, info in zip(tiles, tile_info):
        pos = info["position"]
        row, col = pos[-2], pos[-1]
        tile_slice = tile
        if z_index is not None and tile.ndim == 3:
            tile_slice = tile[z_index]

        row_offset = (-row) % downscale
        col_offset = (-col) % downscale
        tile_ds = tile_slice[row_offset::downscale, col_offset::downscale].astype(
            np.float32, copy=False
        )

        row_ds = (row + row_offset) // downscale
        col_ds = (col + col_offset) // downscale

        row_end = min(row_ds + tile_ds.shape[0], height_ds)
        col_end = min(col_ds + tile_ds.shape[1], width_ds)
        tile_ds = tile_ds[: row_end - row_ds, : col_end - col_ds]

        reconstructed[row_ds:row_end, col_ds:col_end] += tile_ds
        weight_map[row_ds:row_end, col_ds:col_end] += 1

    weight_map[weight_map == 0] = 1
    reconstructed = reconstructed / weight_map
    return reconstructed


def save_threshold_plots_pdf(
    *,
    image_save_dir: Path,
    channel: str,
    grid_tile_info: list[dict],
    recon_tile_info: list[dict],
    tiles: list[np.ndarray],
    large_img: np.ndarray,
    tile_maxes: list[float],
    thr: float,
    is_3d: bool,
) -> Path:
    image_save_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = image_save_dir / f"plot_tile_thresholds_{channel}.pdf"

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
        sns.histplot(tile_maxes, bins=100, kde=True, ax=ax)
        ax.axvline(thr, c="red", linestyle="--", label="Threshold")
        ax.legend()
        ax.set_title(f"Tile max intensities ({channel})")
        pdf.savefig(fig)
        plt.close(fig)

        fig = visualize_tiling_grid(
            grid_tile_info, large_img.shape, color_prop="max_intensity"
        )
        pdf.savefig(fig)
        plt.close(fig)

        fig = visualize_tiling_grid(
            grid_tile_info, large_img.shape, color_prop="thresholded"
        )
        pdf.savefig(fig)
        plt.close(fig)

        downscale = choose_downscale(large_img.shape)
        z_index = large_img.shape[0] // 2 if is_3d else None

        original_ds = downsample_original_slice(large_img, downscale, z_index)
        reconstructed_ds = reconstruct_downsampled_from_tiles(
            tiles,
            recon_tile_info,
            large_img.shape,
            downscale,
            z_index=z_index,
        )
        diff = np.abs(original_ds - reconstructed_ds)

        fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
        ax.imshow(diff, cmap="Grays_r")
        title = f"Downsampled diff (1/{downscale} scale)"
        if is_3d:
            title += f", z={z_index}"
        ax.set_title(title)
        ax.axis("off")
        pdf.savefig(fig)
        plt.close(fig)

    return pdf_path
