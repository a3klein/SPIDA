"""Tiled image preparation: tile a mosaic, select tiles, transform each, stitch back.

`tile_selection` decides which tiles are processed; `transforms` decides what happens
to each one.
"""
from spida.S.image_prep.config import get_optics
from spida.S.image_prep.orchestrator import transform_image
from spida.S.image_prep.tile_selection import select_tiles
from spida.S.image_prep.transforms import TRANSFORMS, resolve_transform

__all__ = [
    "transform_image",
    "get_optics",
    "select_tiles",
    "TRANSFORMS",
    "resolve_transform",
]
