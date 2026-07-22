"""Normalize a segmenter's raw boundary output into the segmentation schema:
``boundaries_micron.parquet``.

Single spec-dispatched entry point, ``ingest_polygons``. Lives at the segmentation
top level (not under ``backends/``) because it is the bridge between a backend's
raw output and the method-agnostic segmentation schema — it relies on the backend
registry/specs but is not itself a segmenter.

Segmentation schema boundary layout (VPT names): ``ID, Geometry, EntityID, ParentID,
ParentType, Type, ZIndex, ZLevel`` in micron space, one polygon per cell per z-plane.
``ZLevel`` is 1-based: ``(z + 1) * micron_per_z`` (every plane owns a full slab).
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.affinity import affine_transform

from .backends.base import SegmentationClass
from .segmentation_utils import (
    convert_geometry,
    _read_micron_to_mosaic_transform,
    _to_multipoly,
    _get_valid,
    _set_process_id,
    _get_id,
    _GEOM,
    _ENTITY_COL,
    _ZINDEX_COL,
)

logger = logging.getLogger(__package__)

_SCHEMA_COLS = ["ID", _GEOM, _ENTITY_COL, "ParentID", "ParentType", "Type",
             _ZINDEX_COL, "ZLevel"]


def _finalize(g: gpd.GeoDataFrame, entity_ids, entity_type: str) -> gpd.GeoDataFrame:
    """Attach the segmentation-schema bookkeeping columns and order them."""
    g = g.copy()
    g[_ENTITY_COL] = np.asarray(entity_ids, dtype=np.int64)
    g["ID"] = np.arange(len(g), dtype=np.int64)
    g["Type"] = entity_type
    for c in ("ParentID", "ParentType", "Name"):
        g[c] = None
    return g[[*_SCHEMA_COLS]]


def _ingest_3d_pixel(pq, micron_to_mosaic, micron_per_z, entity_type, zcol):
    """cellpose 3D: pixel-space, already multi-plane. Micron-transform + assign
    ZIndex/ZLevel/EntityID (no harmonize, no replicate)."""
    geoms = gpd.GeoSeries([_to_multipoly(_get_valid(s)) for s in pq.geometry.values])
    g = gpd.GeoDataFrame(
        {"_cell": pq["ID"].to_numpy(), _ZINDEX_COL: pq[zcol].to_numpy(), _GEOM: geoms},
        geometry=_GEOM,
    )
    if micron_to_mosaic is not None:
        mat = np.linalg.inv(_read_micron_to_mosaic_transform(micron_to_mosaic))
        flat = [*mat[:2, :2].flatten(), *mat[:2, 2].flatten()]
        g[_GEOM] = g[_GEOM].apply(lambda x: _to_multipoly(affine_transform(x, flat)))
        g = g.set_geometry(_GEOM)
    g["ZLevel"] = (g[_ZINDEX_COL] + 1) * micron_per_z
    # one EntityID per cell (consistent across its z-planes)
    _set_process_id()
    id_map = {c: _get_id() for c in np.sort(g["_cell"].unique())}
    entity_ids = g["_cell"].map(id_map).to_numpy()
    return _finalize(g.drop(columns=["_cell"]), entity_ids, entity_type)


def _ingest_proseg(boundaries_path, micron_per_z, entity_type, cell_col, z_col):
    """proseg: micron-space per-layer geojson. Assign ZIndex/ZLevel/EntityID; no
    transform. EntityID stays tied to proseg's cell id (counts/tx are keyed by it)."""
    path = str(boundaries_path)
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            shapes = gpd.read_file(f)
    else:
        shapes = gpd.read_file(path)
    geoms = gpd.GeoSeries([_to_multipoly(_get_valid(s)) for s in shapes.geometry.values])
    g = gpd.GeoDataFrame(
        {"_cell": shapes[cell_col].to_numpy(), _ZINDEX_COL: shapes[z_col].to_numpy(),
         _GEOM: geoms},
        geometry=_GEOM,
    )
    g["ZLevel"] = (g[_ZINDEX_COL] + 1) * micron_per_z
    # EntityID IS proseg's own cell id, so the segmentation schema boundaries align with proseg's
    # native counts/metadata/transcripts (all keyed by `cell`) without a remap.
    entity_ids = g["_cell"].astype(np.int64)
    return _finalize(g.drop(columns=["_cell"]), entity_ids, entity_type)


def ingest_polygons(
    spec: SegmentationClass,
    boundaries_path: str | Path,
    output_path: str | Path,
    *,
    micron_to_mosaic: str | Path | None = None,
    micron_per_z: float = 1.5,
    n_z_planes: int = 7,
    entity_type: str = "cell",
) -> gpd.GeoDataFrame:
    """Normalize ``boundaries_path`` (raw backend output for ``spec``) to the segmentation schema
    and write ``boundaries_micron.parquet`` at ``output_path``.

    Dispatch:
      * cellpose/mesmer 2D  -> convert_geometry (harmonize + replicate + transform)
      * cellpose 3D         -> micron transform + ZIndex/ZLevel/EntityID
      * proseg (v2/v3)      -> geojson layers -> parquet (already micron)
    """
    logger.info("ingest_polygons: %s (v%s)", spec.name, spec.version)

    if spec.name == "proseg":
        out = _ingest_proseg(
            boundaries_path, micron_per_z, entity_type,
            cell_col=spec.columns.cell_id, z_col=spec.columns.z_index,
        )
    elif spec.name in ("cellpose", "mesmer"):
        pq = gpd.read_parquet(boundaries_path)
        zcol = "z" if "z" in pq.columns else _ZINDEX_COL
        is_3d = zcol in pq.columns and pq[zcol].nunique() > 1
        if is_3d:
            out = _ingest_3d_pixel(pq, micron_to_mosaic, micron_per_z, entity_type, zcol)
        else:
            # 2D: the validated convert_geometry does harmonize + replicate + transform
            return convert_geometry(
                boundaries_path, output_path, micron_to_mosaic,
                n_z_planes=n_z_planes, spacing_z=micron_per_z, entity_type=entity_type,
            )
    else:
        raise ValueError(f"ingest_polygons: unsupported method {spec.name!r}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output_path)
    logger.info("ingest_polygons: wrote %d rows -> %s", len(out), output_path)
    return out
