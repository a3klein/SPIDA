from __future__ import annotations

from pathlib import Path

import numpy as np
import imageio.v3 as iio
import tifffile

from spida.utilities.tiling import (
    tile_image_with_overlap,
    save_tiles,
    reconstruct_image_from_tile_files,
)
from spida.S.decon_script import write_tiff
from spida.pl.decon_plots import save_threshold_plots_pdf


def _make_2d_image(h: int = 64, w: int = 72) -> np.ndarray:
    y = np.linspace(0, 1, h, dtype=np.float32)[:, None]
    x = np.linspace(0, 1, w, dtype=np.float32)[None, :]
    img = (1000 * (0.3 * y + 0.7 * x)).astype(np.uint16)
    return img


def _make_3d_image(d: int = 5, h: int = 48, w: int = 56) -> np.ndarray:
    base = _make_2d_image(h, w).astype(np.uint16)
    stack = np.stack([base + (i * 10) for i in range(d)], axis=0)
    return stack


def test_reconstruct_image_from_tile_files_2d_roundtrip(tmp_path: Path):
    img = _make_2d_image(64, 72)

    tiles, tile_info = tile_image_with_overlap(
        img, tile_size=32, overlap=8, color="DAPI"
    )
    saved = save_tiles(tiles, tile_info, tmp_path, prefix="tile", format="tif")

    # Create "deconvolved" 2D tile outputs (identity transform)
    for src in saved:
        dst = src.with_suffix(".decon.2d.tif")
        iio.imwrite(dst, iio.imread(src))

    recon = reconstruct_image_from_tile_files(
        output_dir=tmp_path,
        tile_info=tile_info,
        original_shape=img.shape,
        suffix=".decon.2d",
        match_pre=False,
    )

    assert recon.shape == img.shape
    assert np.array_equal(recon, img)


def test_reconstruct_image_from_tile_files_3d_by_z_slice(tmp_path: Path):
    img = _make_3d_image(d=5, h=48, w=56)

    tiles, tile_info = tile_image_with_overlap(
        img, tile_size=(2, 32, 32), overlap=(0, 8, 8), color="DAPI"
    )
    saved = save_tiles(tiles, tile_info, tmp_path, prefix="tile", format="tif")

    # Create "deconvolved" 3D tile outputs (identity transform)
    for src in saved:
        dst = src.with_suffix(".decon.tif")
        iio.imwrite(dst, iio.imread(src))

    # Reconstruct one z-slice at a time (streaming) and compare
    depth = img.shape[0]
    for z in range(depth):
        z_slice = reconstruct_image_from_tile_files(
            output_dir=tmp_path,
            tile_info=tile_info,
            original_shape=img.shape,
            suffix=".decon",
            match_pre=False,
            z_index=z,
        )
        assert z_slice.shape == img.shape[1:]
        assert np.array_equal(z_slice, img[z])


def test_write_tiff_uses_bigtiff_when_threshold_exceeded(tmp_path: Path):
    img = np.zeros((8, 8), dtype=np.uint16)
    out_path = tmp_path / "bigtiff.tif"

    write_tiff(out_path, img, bigtiff_threshold=1)

    with tifffile.TiffFile(out_path) as tif:
        assert tif.is_bigtiff is True


def test_write_tiff_uses_classic_tiff_below_threshold(tmp_path: Path):
    img = np.zeros((8, 8), dtype=np.uint16)
    out_path = tmp_path / "classic.tif"

    write_tiff(out_path, img, bigtiff_threshold=10**12)

    with tifffile.TiffFile(out_path) as tif:
        assert tif.is_bigtiff is False


def test_threshold_plots_pdf_created_for_3d(tmp_path: Path):
    img = _make_3d_image(d=3, h=64, w=72)

    tiles, tile_info = tile_image_with_overlap(
        img, tile_size=(3, 32, 32), overlap=(0, 8, 8), color="DAPI"
    )

    for tile, info in zip(tiles, tile_info):
        info["max_intensity"] = float(np.max(tile))
    tile_maxes = [info["max_intensity"] for info in tile_info]
    thr = 0.0
    for info in tile_info:
        info["thresholded"] = True

    pdf_path = save_threshold_plots_pdf(
        image_save_dir=tmp_path,
        channel="DAPI",
        grid_tile_info=tile_info,
        recon_tile_info=tile_info,
        tiles=tiles,
        large_img=img,
        tile_maxes=tile_maxes,
        thr=thr,
        is_3d=True,
    )

    assert pdf_path.exists()
    assert pdf_path.stat().st_size > 0
