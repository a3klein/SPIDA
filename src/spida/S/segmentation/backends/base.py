"""Segmentation backend specification.

A ``SegmentationClass`` describes one segmentation method: which pixi env its
backend runs in, which raw artifacts it produces, how those map to the
(VPT-named) segmentation schema, and — derived from what it produces —
which post-processing steps it still needs.

The segmentation schema uses VPT's names on disk
(``EntityID``, ``ZIndex``, ``ZLevel``, ``Geometry``; transcript ``global_x/y/z``,
``gene``, ``barcode_id``) so the VPT binary remains usable as a secondary backend.
"""

from __future__ import annotations

import os
import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__package__)

# segmentation schema output filenames (identical across methods; written by the
# process-segmentation steps). The loader resolves these first, then falls back
# to each method's legacy name for read-compatibility with older SPIDA data.
SCHEMA_BOUNDARIES = "boundaries_micron.parquet"
SCHEMA_CELL_BY_GENE = "cell_by_gene.csv"
SCHEMA_CELL_METADATA = "cell_metadata.csv"
SCHEMA_TRANSCRIPTS = "detected_transcripts.csv"
SCHEMA_SUM_SIGNALS = "sum_signals.csv"


@dataclass
class ColumnMap:
    """Names a method uses in its RAW output for the segmentation schema fields.

    The ingest/loader steps rename these to the segmentation schema standard. Only the fields a
    method actually emits are meaningful (e.g. transcript fields are irrelevant
    for a boundary-only segmenter).
    """

    cell_id: str = "ID"                # -> EntityID  (cellpose/mesmer "ID", proseg "cell")
    geometry: str = "Geometry"         # -> Geometry
    z_index: str | None = None         # -> ZIndex    (proseg "layer"; None => 2D / assigned)
    # transcript-table columns (only if the method emits assigned transcripts)
    tx_gene: str = "gene"
    tx_x: str = "global_x"
    tx_y: str = "global_y"
    tx_z: str = "global_z"
    tx_assignment: str | None = None   # -> cell_id   (proseg "assignment")


@dataclass
class SegmentationClass:
    """Specification for one segmentation method (see module docstring)."""

    name: str
    version: str | None = None

    # ---- backend behavior ----
    requires_gpu: bool = False
    env: str = "cellpose"              # pixi env the *backend* (segment step) runs in
    dim: str = "auto"                  # "2d" | "3d" | "auto" (detect from input)

    # ---- artifacts the method PRODUCES (region-relative names; None == absent) ----
    boundaries_file: str | None = None    # raw polygons (all methods)
    counts_file: str | None = None        # cell-by-gene, if the method computes it
    metadata_file: str | None = None      # per-cell metadata (volume…), if native
    transcripts_file: str | None = None   # assigned transcripts, if native
    spatialdata_file: str | None = None   # e.g. proseg v3 zarr (bundles counts+tx)

    # ---- raw->segmentation schema name mapping ----
    columns: ColumnMap = field(default_factory=ColumnMap)

    # ---- transcript source when the method does NOT assign transcripts ----
    raw_transcripts_file: str = "detected_transcripts.parquet"

    # ---- legacy segmentation schema boundary name (for loader read back-compat) ----
    legacy_boundaries_file: str | None = None

    # ---- capabilities, derived from what's populated ----
    @property
    def provides_assignment(self) -> bool:
        return self.transcripts_file is not None or self.spatialdata_file is not None

    @property
    def provides_counts(self) -> bool:
        return self.counts_file is not None or self.spatialdata_file is not None

    @property
    def provides_metadata(self) -> bool:
        return self.metadata_file is not None or self.spatialdata_file is not None

    @property
    def needs(self) -> list[str]:
        """Post-processing steps this method still requires, inferred from what
        it does NOT provide. ``sum_signals`` always runs (no segmenter computes
        image intensity)."""
        steps = ["ingest_polygons"]                  # raw -> segmentation schema boundaries (per-method)
        if not self.provides_assignment:
            steps.append("partition_transcripts")    # build cell-by-gene
        if not self.provides_metadata:
            steps.append("derive_entity_metadata")   # volume/centroid/bounds/…
        steps.append("sum_signals")                  # image intensity: always
        return steps

    def require_env(self, step: str = "segment") -> None:
        """Assert the current pixi env matches the one this step needs.

        The ``segment`` step runs in ``self.env``; all post-processing runs in
        ``preprocessing``. No-op when ``PIXI_ENVIRONMENT_NAME`` is unset (e.g.
        outside pixi).
        """
        required = self.env if step == "segment" else "preprocessing"
        current = os.environ.get("PIXI_ENVIRONMENT_NAME")
        if current is not None and current != required:
            raise RuntimeError(
                f"{self.name} '{step}' step must run in the '{required}' pixi env "
                f"(got '{current}'). Try: pixi run -e {required} ..."
            )
