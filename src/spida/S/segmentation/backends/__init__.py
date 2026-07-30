"""Segmentation backends: per-method specs + (later) backend implementations."""

from .base import (
    ColumnMap,
    SegmentationClass,
    resolve_region_dir,
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

_CONFIG_CLASSES = {
    "cellpose": (".cellpose_config", "CellposeConfig"),
}


def get_config_class(method: str):
    """Return the ``{Method}Config`` class for ``method``, or None if it has none yet."""
    entry = _CONFIG_CLASSES.get(method)
    if entry is None:
        return None
    from importlib import import_module

    module_name, class_name = entry
    return getattr(import_module(module_name, __package__), class_name)


def describe_method_params(method: str) -> str:
    """Render the run-parameter reference for ``method`` (``--list-params``)."""
    try:
        get_spec(method)
    except KeyError:
        known = sorted({name for name, _ in REGISTRY})
        return (f"Unknown segmentation method {method!r}. Known methods: "
                f"{', '.join(known)}.")

    config_class = get_config_class(method)
    if config_class is None:
        have = ", ".join(sorted(_CONFIG_CLASSES)) or "none"
        return (
            f"'{method}' does not have a run-configuration class yet, so its "
            f"parameters are not self-describing. Methods with one: {have}.\n"
            f"See the {method} backend module for the parameters it accepts."
        )
    return config_class.describe()


__all__ = [
    "ColumnMap",
    "SegmentationClass",
    "get_config_class",
    "describe_method_params",
    "resolve_region_dir",
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
