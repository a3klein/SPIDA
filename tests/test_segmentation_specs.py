"""Tests for the segmentation backend specs (SegmentationClass registry).

The `needs`/capabilities contract is what drives the whole pipeline: it must
report the right post-processing steps per method based on what each provides.
"""

import pytest

from spida.S.segmentation.backends import (
    ColumnMap,
    SegmentationClass,
    get_spec,
    REGISTRY,
)


def test_boundary_only_needs_full_postprocess():
    """cellpose/mesmer provide nothing -> need partition + derive + sum."""
    for name in ("cellpose", "mesmer"):
        spec = get_spec(name)
        assert spec.needs == [
            "ingest_polygons", "partition_transcripts",
            "derive_entity_metadata", "sum_signals",
        ]
        assert not spec.provides_counts
        assert not spec.provides_metadata
        assert not spec.provides_assignment
        assert spec.env == "cellpose"


def test_proseg_v2_needs_only_ingest_and_signals():
    spec = get_spec("proseg", "2")
    assert spec.provides_counts
    assert spec.provides_metadata
    assert spec.provides_assignment
    assert spec.needs == ["ingest_polygons", "sum_signals"]
    assert spec.env == "preprocessing"


def test_proseg_v3_bundles_via_zarr():
    spec = get_spec("proseg", "3")
    # counts + assignment come from the zarr, metadata from csv (+zarr)
    assert spec.spatialdata_file == "proseg_outputs.zarr"
    assert spec.provides_counts and spec.provides_assignment and spec.provides_metadata
    assert spec.needs == ["ingest_polygons", "sum_signals"]


def test_proseg_default_version_is_v3():
    assert get_spec("proseg") is REGISTRY[("proseg", "3")]


def test_sum_signals_always_needed():
    for spec in REGISTRY.values():
        assert spec.needs[-1] == "sum_signals"
        assert spec.needs[0] == "ingest_polygons"


def test_unknown_spec_raises():
    with pytest.raises(KeyError):
        get_spec("does-not-exist")


def test_capabilities_derive_from_populated_fields():
    """A hand-built spec: populating an artifact flips the derived capability."""
    bare = SegmentationClass(name="x", columns=ColumnMap())
    assert not bare.provides_counts and "partition_transcripts" in bare.needs

    with_counts = SegmentationClass(name="x", counts_file="c.csv",
                                    transcripts_file="t.csv", metadata_file="m.csv")
    assert with_counts.provides_counts
    assert with_counts.needs == ["ingest_polygons", "sum_signals"]


def test_require_env(monkeypatch):
    spec = get_spec("cellpose")
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "preprocessing")
    with pytest.raises(RuntimeError):
        spec.require_env("segment")          # cellpose backend needs cellpose env
    spec.require_env("process")              # post-processing wants preprocessing -> ok
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "cellpose")
    spec.require_env("segment")              # ok
