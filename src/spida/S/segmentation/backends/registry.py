"""Registry of segmentation-method specifications.

One ``SegmentationClass`` instance per method; ``get_spec(name, version)`` resolves
a method (+ optional version, for proseg) to its spec.
"""

from __future__ import annotations

from spida.io.formats.segmenters import (
    PROSEG_CELL_METADATA, PROSEG_POLYGON_LAYERS, PROSEG_SPATIALDATA, PROSEG_V2_COUNTS, PROSEG_V2_TRANSCRIPTS,
)
from .base import ColumnMap, SegmentationClass

# --- boundary-only segmenters (pixel-space masks -> need partition/derive/sum) ---
CELLPOSE = SegmentationClass(
    name="cellpose",
    requires_gpu=True,
    env="cellpose",
    dim="auto",                                   # 2D or 3D detected from input
    boundaries_file="polygons.parquet",
    native_format="cellpose",
    columns=ColumnMap(cell_id="ID"),
    legacy_boundaries_file="cellpose_micron_space.parquet",
)

MESMER = SegmentationClass(
    name="mesmer",
    requires_gpu=True,
    env="cellpose",                               # Change to deepcell env when that one becomes active
    dim="2d",
    boundaries_file="polygons.parquet",
    native_format="cellpose",
    columns=ColumnMap(cell_id="ID"),
    legacy_boundaries_file="cellpose_micron_space.parquet",
)

# --- assignment-native segmenter (micron-space; provides counts/metadata/tx) ---
PROSEG_V2 = SegmentationClass(
    name="proseg",
    version="2",
    requires_gpu=False,
    env="preprocessing",                          # just shells the rust binary
    dim="3d",
    boundaries_file=PROSEG_POLYGON_LAYERS,
    boundaries_space="micron",
    native_format="proseg",
    counts_file=PROSEG_V2_COUNTS,
    metadata_file=PROSEG_CELL_METADATA,
    transcripts_file=PROSEG_V2_TRANSCRIPTS,
    columns=ColumnMap(
        cell_id="cell", z_index="layer",
        tx_x="x", tx_y="y", tx_z="z", tx_assignment="assignment",
    ),
    legacy_boundaries_file="proseg_polygons.parquet",
)

PROSEG_V3 = SegmentationClass(
    name="proseg",
    version="3",
    requires_gpu=False,
    env="preprocessing",
    dim="3d",
    boundaries_file=PROSEG_POLYGON_LAYERS,
    boundaries_space="micron",
    native_format="proseg",
    metadata_file=PROSEG_CELL_METADATA,
    spatialdata_file=PROSEG_SPATIALDATA,       # bundles counts + assigned transcripts
    columns=ColumnMap(
        cell_id="cell", z_index="layer",
        tx_x="x", tx_y="y", tx_z="z", tx_assignment="assignment",
    ),
    legacy_boundaries_file="proseg_polygons.parquet",
)

REGISTRY: dict[tuple[str, str | None], SegmentationClass] = {
    ("cellpose", None): CELLPOSE,
    ("mesmer", None): MESMER,
    ("proseg", "2"): PROSEG_V2,
    ("proseg", "3"): PROSEG_V3,
}


def get_spec(name: str, version: str | None = None) -> SegmentationClass:
    """Resolve a method (+ optional version) to its spec.

    For proseg, ``version`` selects v2/v3 (defaults to v3, the current release).
    """
    if name == "proseg" and version is None:
        version = "3"
    try:
        return REGISTRY[(name, version)]
    except KeyError:
        raise KeyError(
            f"no segmentation spec for name={name!r} version={version!r}; "
            f"known: {sorted(REGISTRY)}"
        )
