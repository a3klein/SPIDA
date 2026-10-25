"""
Image reading. windowed / downscaled / masked reads of 2D TIFF planes, stacked into (C, Y, X, Z).
read_stains always returns 4D; Z = 1 for a single plane.
"""
from __future__ import annotations
import logging
from pathlib import Path
import numpy as np

logger = logging.getLogger(__package__)

# Rows of a plane read at full resolution at once when streaming a downscale read. bounds peak memory.
DEFAULT_STRIP_ROWS = 1024


def probe_mosaic_shape(path: str | Path) -> tuple[int, int]:
    """Returns (rows, cols) of a mosaic TIFF, read from the header only."""
    import tifffile

    with tifffile.TiffFile(str(path)) as tf:
        shape = tuple(int(n) for n in tf.series[0].shape)
    if len(shape) != 2:
        raise ValueError(f"{path}: expected a 2D mosaic plane, got shape {shape}")
    return shape  # type: ignore[return-value]


def read_plane_window(
    path: Path,
    y0: int,
    y1: int,
    x0: int,
    x1: int,
    downscale: int | None = None,
    strip_rows: int = DEFAULT_STRIP_ROWS,
    mask: np.ndarray | None = None
) -> np.ndarray:
    """
    A plane's [y0:y1, x0:x1] window as a 2D array.
    Downscaling happens during the read.
    Mask is applied before downscale to each strip.
    Falls back to a full read if the plane is not memmap-able, with warning.

    Parameters: 
        path - mosaic TIFF path
        y0, y1, x0, x1 - window coordinates
        downscale - xy reduction factor, applied during the read
        strip_rows - Number of rows held at once while streaming a downscale
        mask - full-resolution, mosaic-aligned boolean tissue mask; out-of-tissue pixels.
    
    Returns: 
        2D array of the window, downscaled if requested.
    """
    import tifffile
    import skimage as ski

    try:
        src = tifffile.memmap(str(path), mode="r")
    except (ValueError, MemoryError):
        logger.warning(
            "%s is not memmap-able (compressed?); falling back to a full read, "
            "which needs the whole plane in RAM.", Path(path).name,
        )
        src = tifffile.imread(str(path))

    try:
        if not downscale or downscale <= 1:
            # copy when masking: memmap slices are read-only
            win = np.array(src[y0:y1, x0:x1]) if mask is not None else np.asarray(src[y0:y1, x0:x1])
            if mask is not None:
                win[~mask[y0:y0 + win.shape[0], x0:x0 + win.shape[1]]] = 0
            return win

        # Strip height must be a multiple of downscale.
        # only the final strip may carry a partial block.
        step = max(downscale, (strip_rows // downscale) * downscale)
        chunks = []
        row = y0
        while row < y1:
            end = min(row + step, y1)
            if y1 - end < downscale: # absorb a short tail into this strip
                end = y1
            if mask is None:
                chunk = np.asarray(src[row:end, x0:x1])
            else:
                chunk = np.array(src[row:end, x0:x1])
                sub = mask[row:row + chunk.shape[0], x0:x0 + chunk.shape[1]]
                chunk[~sub] = 0
            chunks.append(ski.transform.downscale_local_mean(
                chunk, (downscale, downscale)))
            row = end
        return np.concatenate(chunks, axis=0) if len(chunks) > 1 else chunks[0]
    finally:
        del src


def read_stains(
    stain_planes: dict[str, list[Path]],
    window: tuple[int, int, int, int] | None = None,
    downscale: int | None = None,
    mask: np.ndarray | None = None,
    strip_rows: int = DEFAULT_STRIP_ROWS
) -> np.ndarray:
    """Read every stain into one (C, Y, X, Z) array.

    Parameters
    ----------
    stain_planes
        stain name -> z-ordered plane paths, as returned by :func:`find_stain_planes`.
    window
        (y0, y1, x0, x1) for a windowed read; None reads whole planes.
    downscale
        xy reduction factor, applied during the read.
    mask
        full-resolution, mosaic-aligned boolean tissue mask.
    strip_rows
        full-resolution rows held at once while streaming a downscale.

    Returns 4D (C, Y, X, Z), with Z = 1 for a single-plane stain.
    """
    if not stain_planes:
        raise ValueError("stain_planes is empty; nothing to read")

    if window is None:
        first = next(iter(stain_planes.values()))[0]
        rows, cols = probe_mosaic_shape(first)
        y0, y1, x0, x1 = 0, rows, 0, cols
    else:
        y0, y1, x0, x1 = window

    channels = []
    for stain, planes in stain_planes.items():
        planes = list(planes)
        if not planes:
            raise ValueError(f"stain {stain!r} has no planes")
        stack = [
            read_plane_window(p, y0, y1, x0, x1, downscale=downscale,
                              strip_rows=strip_rows, mask=mask)
            for p in planes
        ]
        channels.append(np.stack(stack, axis=-1))          # (Y, X, Z)

    shapes = {c.shape for c in channels}
    if len(shapes) != 1:
        raise ValueError(
            "stains disagree on shape: "
            + ", ".join(f"{s}={c.shape}" for s, c in zip(stain_planes, channels))
            + ". Every stain must have the same plane dimensions and z-count."
        )
    return np.stack(channels, axis=0)                      # (C, Y, X, Z)
