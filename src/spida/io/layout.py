"""Where a dataset lives, and what the elements are called."""

from __future__ import annotations
import re
import json
import logging
from pathlib import Path
from dataclasses import dataclass, field

# SPIDA imports
from spida._constants import IMAGE_KEY, POINTS_KEY, SHAPES_KEY, TABLE_KEY

logger = logging.getLogger(__name__)

# TODO: cleaner naming for segmentation dirs only? 
def resolve_region_dir(
    store,
    exp_name,
    reg_name,
    label,
    *,
    must_exist: bool = False
) -> Path:
    """
    One region's segmentation directory: ``{store}/{exp}/{reg}/{label}``.
    Legacy ok: ``{store}/{exp}/{label}/{reg}``

    If ``must_exist``, falls back to the legacy when only it exists.
    """
    current = Path(store) / exp_name / reg_name / label
    if must_exist and not current.exists():
        legacy = Path(store) / exp_name / label / reg_name
        if legacy.exists():
            logger.info("resolve_region_dir: using legacy layout %s", legacy)
            return legacy
    return current


def element_keys(
    prefix: str, # segmentatio prefix, e.g. "default" or "cellpose"
    namespace: str | None # the <exp>_<reg> part of element names, e.g. "exp1_regionA"
) -> dict:
    """Element names for one prefix. Empty/None ``namespace`` leaves out the ``<exp>_<reg>`` part."""
    ns = f"_{namespace}" if namespace else ""
    return {
        SHAPES_KEY: f"{prefix}{ns}_polygons",
        POINTS_KEY: f"{prefix}{ns}_transcripts",
        TABLE_KEY: f"{prefix}_table",
        IMAGE_KEY: f"default{ns}_z3",
    }

def table_key(prefix: str, suffix: str = "") -> str:
    """A table's name; tables carry no ``<exp>_<reg>`` part, so no store is needed."""
    return f"{element_keys(prefix, None)[TABLE_KEY]}{suffix}"

# the store records its <exp>_<reg> in attrs; fall back to the element names.
def read_namespace(
    zarr_path: Path
) -> str:
    """The ``<exp>_<reg>`` a store was ingested with: ``attrs["spida"]["namespace"]``, else from its element names."""
    zarr_path = Path(zarr_path)
    for meta in (zarr_path / "zarr.json", zarr_path / ".zattrs"):
        if meta.exists():
            attrs = json.loads(meta.read_text())
            attrs = attrs.get("attributes", attrs)
            ns = attrs.get("spida", {}).get("namespace")
            if ns is not None:
                return ns
    shapes = zarr_path / "shapes"
    names = sorted(p.name for p in shapes.iterdir()) if shapes.is_dir() else []
    found = [m.group(1) or "" for n in names if (m := re.fullmatch(r"default(?:_(.+))?_polygons", n))]
    if len(found) == 1:
        return found[0]
    raise ValueError(
        f"cannot tell which experiment/region {zarr_path} was ingested as: no attrs['spida'] and "
        f"default_*_polygons candidates {found or 'none'}; pass --exp/--reg"
    )

def _recorded_namespace(zarr_path: Path | None) -> str | None:
    """The ``<exp>_<reg>`` an existing store was ingested with, or None if there is no store yet."""
    if zarr_path is None or not Path(zarr_path).exists():
        return None
    try:
        return read_namespace(zarr_path)
    except ValueError:
        return None


# Class to store filepaths for a dataset
def region_label(
    region: str,
    labels: dict | None = None,
) -> str:
    """
    A region's label (its donor): ``labels[region]`` if given, else the name without ``region_``.
    Parameters:
    region (str): Region name, e.g. region_UCI5224.
    labels (dict, optional): {region name: label}.
    """
    if labels and region in labels:
        return labels[region]
    return region[len("region_"):] if region.startswith("region_") else region


def load_region_labels(
    labels,
) -> dict | None:
    """
    A region-label mapping from a dict, a JSON file path, or None.
    Parameters:
    labels (dict | str | Path | None): The mapping or its JSON file.
    """
    if labels is None or isinstance(labels, dict):
        return labels
    import json

    return json.loads(Path(labels).read_text())


@dataclass(frozen=True)
class RegionMeta:
    """Names a compute step writes into its outputs (obs columns, uns, figure titles)."""

    experiment: str = ""
    region: str = ""
    segmentation: str = ""
    donor: str | None = None


@dataclass(frozen=True)
class Dataset:
    exp: str | None = None
    reg: str | None = None
    zarr: Path | None = None # zarr store location
    raw: Path | None = None # where raw data lives
    seg: Path | None = None # this dataset's segmentations, one dir per label
    figures: Path | None = None # this dataset's figures, one dir per label
    anndata: Path | None = None # h5ad copies of this dataset's tables, one file per table key
    annotations: Path | None = None # annotation intermediates, {table_key}/{method}/{reference}
    seg_legacy: Path | None = None # read-only fallback {seg_store}/{exp}/{label}/{reg}
    region_labels: dict | None = field(default=None, compare=False, hash=False)  # {region name: label}

    @property
    def named(self) -> bool: # if both exp and reg are set.
        return bool(self.exp and self.reg)

    @property
    def namespace(self) -> str:
        if not self.named:
            return read_namespace(self.require("zarr"))
        ns = f"{self.exp}_{self.reg}"
        recorded = _recorded_namespace(self.zarr)
        if recorded is not None and recorded != ns:
            raise ValueError(f"{self.zarr} was ingested as {recorded!r}, but --exp/--reg give {ns!r}")
        return ns

    def keys(self, prefix: str) -> dict:
        return element_keys(prefix, self.namespace)

    def require(self, name: str) -> Path:
        value = getattr(self, name)
        if value is None:
            raise ValueError(f"dataset has no {name} path. Pass --{name}_path or configure its store root")
        return value

    def meta(self, segmentation: str, labels: dict | None = None) -> RegionMeta:
        """Names for a compute step on segmentation ``segmentation``; ``labels`` overrides the dataset's mapping."""
        reg = self.reg or ""
        return RegionMeta(self.exp or "", reg, segmentation, region_label(reg, labels or self.region_labels))

    def raw_images(self) -> Path: # input images for dataset
        return self.require("raw") / "images"

    def seg_dir(self, label: str, must_exist: bool = False) -> Path:
        path = self.require("seg") / label
        if must_exist and not path.exists() and self.seg_legacy is not None:
            legacy = self.seg_legacy / label / self.reg
            if legacy.exists():
                logger.info("seg_dir: using legacy layout %s", legacy)
                return legacy
        return path

    def figures_dir(self, label: str) -> Path:
        return self.require("figures") / label

    def table_h5ad(self, table_key: str) -> Path:
        return self.require("anndata") / f"{table_key}.h5ad"

    def annotation_dir(self, table_key: str, method: str, reference: str) -> Path:
        return self.require("annotations") / table_key / method / reference


def dataset_from(
    config,
    exp=None,
    reg=None,
    zarr=None,
    raw=None,
    seg=None,
    figures=None,
    anndata=None,
    annotations=None,
    region_labels=None,
) -> Dataset:
    """
    Entry point for building a Dataset.

    Names resolve paths from the config store roots.
    Explicit paths override location.
    Names set the <exp>_<reg> part of element names in the zarr store.
    """
    if bool(exp) != bool(reg):
        raise ValueError("pass both --exp and --reg, or neither")
    named = bool(exp)
    if not named and all(p is None for p in (zarr, raw, seg, figures, anndata, annotations)):
        raise ValueError("pass --exp and --reg, or the locations the step needs (--zarr_path / --raw_path / "
                         "--seg_path / --figures_path / --anndata_path / --annotations_path)")

    def _tree(root, *parts): # chain together pathlib Paths 
        return Path(root).joinpath(*parts) if named and root is not None else None

    return Dataset(
        exp=exp or None,
        reg=reg or None,
        zarr=Path(zarr) if zarr is not None else _tree(config.zarr_store, exp, reg),
        raw=Path(raw) if raw is not None else _tree(config.root_path, exp, "out", reg), # The "out" here is quirk for my specific filepaths
        seg=Path(seg) if seg is not None else _tree(config.segmentation_store, exp, reg),
        figures=Path(figures) if figures is not None else _tree(config.image_store, exp, reg),
        anndata=Path(anndata) if anndata is not None else _tree(config.anndata_store, exp, reg),
        annotations=Path(annotations) if annotations is not None else _tree(config.annotation_store, exp, reg),
        seg_legacy=None if seg is not None else _tree(config.segmentation_store, exp),
        region_labels=load_region_labels(region_labels if region_labels is not None
                                         else getattr(config, "region_labels", None)),
    )


class TreeLayout:
    """The ``{root}/{exp}/{reg}`` tree."""

    def __init__(self, config):
        self.config = config

    def dataset(self, exp, reg, **paths) -> Dataset:
        return dataset_from(self.config, exp, reg, **paths)

    def units(self, exp: str, under: str = "zarr", label: str | None = None) -> list[str]:
        """Region names of an experiment, found under the zarr store, the raw root, or a segmentation label."""
        if under == "zarr":
            base = Path(self.config.zarr_store) / exp
        elif under == "raw":
            base = Path(self.config.root_path) / exp / "out" # again, "out" here is a personal file quirk that may need resolution. 
        elif under == "seg":
            base = Path(self.config.segmentation_store) / exp
            regs = sorted(p.name for p in base.glob("region_*") if (p / label).is_dir())
            return regs or sorted(p.name for p in (base / label).glob("region_*"))
        else:
            raise ValueError(f"under must be zarr, raw or seg, got {under!r}")
        return sorted(p.name for p in base.glob("region_*") if p.is_dir())
