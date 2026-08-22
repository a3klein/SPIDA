# Main pulls in anndata / spatialdata, which are absent from the cellpose pixi environment,

_LAZY = {
    "ingest_region": ".main",
    "load_segmentation_region": ".main",
}
__all__ = list(_LAZY)

def __getattr__(name: str):
    if name in _LAZY:
        from importlib import import_module
        module = import_module(_LAZY[name], __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(set(globals()) | set(_LAZY))
