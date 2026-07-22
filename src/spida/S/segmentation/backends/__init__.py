"""Segmentation backends: per-method specs + (later) backend implementations."""

from .base import (
    ColumnMap,
    SegmentationClass,
    SCHEMA_BOUNDARIES,
    SCHEMA_CELL_BY_GENE,
    SCHEMA_CELL_METADATA,
    SCHEMA_TRANSCRIPTS,
    SCHEMA_SUM_SIGNALS,
)
from .registry import (
    REGISTRY,
    get_spec,
    CELLPOSE,
    MESMER,
    PROSEG_V2,
    PROSEG_V3,
)

__all__ = [
    "ColumnMap",
    "SegmentationClass",
    "SCHEMA_BOUNDARIES",
    "SCHEMA_CELL_BY_GENE",
    "SCHEMA_CELL_METADATA",
    "SCHEMA_TRANSCRIPTS",
    "SCHEMA_SUM_SIGNALS",
    "REGISTRY",
    "get_spec",
    "CELLPOSE",
    "MESMER",
    "PROSEG_V2",
    "PROSEG_V3",
]
