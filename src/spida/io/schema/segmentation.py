"""
The segmentation schema: five files per segmentation (VPT names, keyed on ``EntityID``), the
standard every backend's output is normalized to and the loader reads.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

logger = logging.getLogger(__package__)

# The loader resolves these first, then falls back to a method's legacy name for older SPIDA data.
SCHEMA_BOUNDARIES = "boundaries_micron.parquet"
SCHEMA_CELL_BY_GENE = "cell_by_gene.csv"
SCHEMA_CELL_METADATA = "cell_metadata.csv"
SCHEMA_TRANSCRIPTS = "detected_transcripts.csv"
SCHEMA_SUM_SIGNALS = "sum_signals.csv"
SCHEMA_FILES = (SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_TRANSCRIPTS,
                SCHEMA_CELL_METADATA, SCHEMA_SUM_SIGNALS)
LEGACY_VPT_BOUNDARIES = "cellpose_micron_space.parquet"   # what the VPT backend writes



def files_written(
    spec=None,
    backend: str = "native",
    vpt_boundaries: str = LEGACY_VPT_BOUNDARIES,
) -> list[str]:
    """
    Schema files a processing step may write into a segmentation's folder (its overwrite guard).
    Parameters:
    spec (SegmentationClass, optional): The method; None for a user-provided segmentation (all five).
    backend (str): "native" or "vpt".
    vpt_boundaries (str): Boundaries name the VPT backend writes.
    """
    if backend == "vpt":
        return [vpt_boundaries, SCHEMA_CELL_BY_GENE, SCHEMA_TRANSCRIPTS, SCHEMA_CELL_METADATA, SCHEMA_SUM_SIGNALS]
    if spec is None:
        return list(SCHEMA_FILES)
    out = [SCHEMA_BOUNDARIES, SCHEMA_SUM_SIGNALS, SCHEMA_CELL_METADATA]
    if "partition_transcripts" in spec.needs:
        out += [SCHEMA_CELL_BY_GENE, SCHEMA_TRANSCRIPTS]
    return out

def resolve(
    seg_dir: str | Path,
    name: str,
    legacy: str | None = None,
) -> Path:
    """
    ``seg_dir / name`` if it exists, else ``seg_dir / legacy`` if that exists, else ``seg_dir / name``.
    Parameters:
    seg_dir (str | Path): One segmentation's folder.
    name (str): The schema file name.
    legacy (str, optional): The older name to fall back to.
    """
    seg_dir = Path(seg_dir)
    path = seg_dir / name
    if not path.exists() and legacy is not None and (seg_dir / legacy).exists():
        logger.info("using legacy file %r (%r absent)", legacy, name)
        return seg_dir / legacy
    return path


def read_boundaries(
    path: str | Path,
):
    """
    A boundaries file (``.parquet``, ``.geojson`` or ``.geojson.gz``) as a GeoDataFrame.
    Parameters:
    path (str | Path): The boundaries file.
    """
    import geopandas as gpd

    path = str(path)
    if path.endswith(".parquet"):
        return gpd.read_parquet(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return gpd.read_file(f)
    return gpd.read_file(path)


def write_boundaries(
    gdf,
    path: str | Path,
) -> Path:
    """
    Write schema boundaries (``boundaries_micron.parquet``).
    Parameters:
    gdf (GeoDataFrame): The boundaries, in schema columns.
    path (str | Path): The boundaries file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(path)
    return path


def write_cell_metadata(
    df,
    seg_dir: str | Path,
) -> Path:
    """
    Write ``cell_metadata.csv`` (indexed by EntityID) into a segmentation's folder.
    Parameters:
    df (pd.DataFrame): Per-cell metadata.
    seg_dir (str | Path): The segmentation's folder.
    """
    path = Path(seg_dir) / SCHEMA_CELL_METADATA
    df.to_csv(path)
    return path
