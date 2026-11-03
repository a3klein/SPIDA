"""Phase-2 guards on the segmentation files: each step guards only the files it writes.

segment-region and process-segmentation-region share {seg}/{exp}/{reg}/{label}/, so a guard on the
folder would block process after segment. Backends and heavy post-processing are stubbed (they write
the same file names); process_custom_segmentation runs for real.
"""
import os
import json
import sys
import types
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
import shapely.geometry as geom
from click.testing import CliRunner

from spida.io.layout import Dataset
from spida.S.segmentation import main as seg_main
from spida.io.schema.segmentation import (
    SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_CELL_METADATA, SCHEMA_SUM_SIGNALS, SCHEMA_TRANSCRIPTS,
)


@pytest.fixture(autouse=True)
def _any_pixi_env(monkeypatch):
    monkeypatch.delenv("PIXI_ENVIRONMENT_NAME", raising=False)  # spec.require_env no-ops


@pytest.fixture
def ds(tmp_path):
    raw = tmp_path / "raw" / "E" / "out" / "region_1"
    (raw / "images").mkdir(parents=True)
    for f in ("images/mosaic_DAPI_z0.tif", "images/micron_to_mosaic_pixel_transform.csv",
              "detected_transcripts.parquet", "detected_transcripts.csv"):      # inputs the steps require
        (raw / f).write_text("")
    for prefix in ("cellpose_nuc", "cellpose_cell"):                            # align reads these shapes
        (tmp_path / "z.zarr" / "shapes" / f"{prefix}_E_region_1_polygons").mkdir(parents=True)
    return Dataset(exp="E", reg="region_1", raw=raw, seg=tmp_path / "seg" / "E" / "region_1",
                   zarr=tmp_path / "z.zarr")


@pytest.fixture
def backends(monkeypatch):
    """Stub cellpose / proseg backends: record calls, write their raw output names."""
    calls = []

    def cellpose(images_dir, output_dir, reg_name, **k):
        calls.append(("cellpose", k))
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        (Path(output_dir) / "polygons.parquet").write_text("x")

    def proseg(transcripts_path, output_dir, **k):
        calls.append(("proseg", k))
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for f in ("cell-polygons-layers.geojson.gz", "cell-metadata.csv.gz"):
            (Path(output_dir) / f).write_text("x")

    for name, fn in (("cellpose", cellpose), ("proseg", proseg)):
        mod = types.ModuleType(f"spida.S.segmentation.backends.{name}")
        setattr(mod, f"run_{name}", fn)
        monkeypatch.setitem(sys.modules, mod.__name__, mod)
    return calls


@pytest.fixture
def post(monkeypatch):
    """Stub the native post-processing steps: each writes the schema file it owns."""
    import spida.S.segmentation.ingest as ing
    import spida.S.segmentation.postprocess as pp
    calls = []

    def ingest(spec, src, dst, **k):
        calls.append("ingest")
        Path(dst).write_text("x")

    def partition(b, t, output_entity_by_gene, output_transcripts, **k):
        Path(output_entity_by_gene).write_text("x")
        Path(output_transcripts).write_text("x")

    def sum_signals(b, images, m2m, output_csv, **k):
        Path(output_csv).write_text("x")
        return pd.DataFrame(index=pd.Index([1], name="EntityID"))

    monkeypatch.setattr(ing, "ingest_polygons", ingest)
    monkeypatch.setattr(pp, "partition_transcripts", partition)
    monkeypatch.setattr(pp, "derive_entity_metadata", lambda *a, **k: pd.DataFrame(index=pd.Index([1], name="EntityID")))
    monkeypatch.setattr(pp, "sum_signals", sum_signals)
    return calls


# ---- segment-region ----

def test_segment_rerun_refused_before_backend(ds, backends):
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    assert len(backends) == 1
    with pytest.raises(FileExistsError, match="polygons.parquet"):
        seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    assert len(backends) == 1  # backend never started
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell", overwrite=True)
    assert len(backends) == 2


def _cellpose_prior(ds):
    (ds.seg / "cellpose_cell").mkdir(parents=True, exist_ok=True)
    (ds.seg / "cellpose_cell" / "detected_transcripts.csv").write_text("")  # what segment-region proseg reads


def test_segment_guard_uses_the_method_spec(ds, backends):
    _cellpose_prior(ds)
    (ds.seg / "proseg_cell").mkdir(parents=True)
    (ds.seg / "proseg_cell" / "polygons.parquet").write_text("x")  # cellpose's name, not proseg's
    seg_main.segment_region("proseg", ds, prefix_name="proseg_cell")
    with pytest.raises(FileExistsError, match="cell-polygons-layers.geojson.gz"):
        seg_main.segment_region("proseg", ds, prefix_name="proseg_cell")


def test_other_labels_do_not_block(ds, backends):
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_nuc")
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    assert len(backends) == 2


# ---- segment -> process in the same folder ----

def test_process_after_segment_not_blocked_then_both_refused(ds, backends, post):
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell")
    out = ds.seg / "cellpose_cell"
    for f in (SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_TRANSCRIPTS, SCHEMA_CELL_METADATA, SCHEMA_SUM_SIGNALS):
        assert (out / f).is_file(), f
    with pytest.raises(FileExistsError) as e:
        seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell")
    assert SCHEMA_BOUNDARIES in str(e.value) and "polygons.parquet" not in str(e.value)
    assert post == ["ingest"]  # re-run stopped before ingest
    with pytest.raises(FileExistsError, match="polygons.parquet"):
        seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell", overwrite=True)
    assert post == ["ingest", "ingest"]


def test_proseg_process_guards_only_what_it_writes(ds, backends, post, monkeypatch):
    # proseg brings its own counts / transcripts: those schema names are not its outputs
    _cellpose_prior(ds)
    seg_main.segment_region("proseg", ds, prefix_name="proseg_cell")
    out = ds.seg / "proseg_cell"
    (out / SCHEMA_CELL_BY_GENE).write_text("x")
    monkeypatch.setattr(pd, "read_csv", lambda *a, **k: pd.DataFrame({"cell": [1]}))
    seg_main.process_segmentation_region("proseg", ds, prefix_name="proseg_cell")
    with pytest.raises(FileExistsError) as e:
        seg_main.process_segmentation_region("proseg", ds, prefix_name="proseg_cell")
    assert SCHEMA_CELL_BY_GENE not in str(e.value)


def test_vpt_backend_guards_its_own_names(ds, backends, monkeypatch):
    import spida.S.segmentation.vpt as vpt
    ran = []
    monkeypatch.setattr(vpt, "seg_to_vpt", lambda *a, **k: ran.append(1))
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    out = ds.seg / "cellpose_cell"
    seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell", backend="vpt")
    (out / "cellpose_micron_space.parquet").write_text("x")
    with pytest.raises(FileExistsError, match="cellpose_micron_space.parquet"):
        seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell", backend="vpt")
    assert ran == [1]


# ---- process-custom-segmentation (real) ----

def _boxes(tmp_path):
    g = gpd.GeoDataFrame([{"my_cell": f"c{i}", "Geometry": geom.box(i * 20, 0, i * 20 + 15, 15)} for i in range(3)],
                         geometry="Geometry")
    p = tmp_path / "bounds.parquet"
    g.to_parquet(p)
    return p


def test_custom_real_rerun(tmp_path):
    out = tmp_path / "out"
    seg_main.process_custom_segmentation(_boxes(tmp_path), out, cell_id_col="my_cell")
    first = gpd.read_parquet(out / SCHEMA_BOUNDARIES)
    t0 = (out / SCHEMA_BOUNDARIES).stat().st_mtime_ns
    with pytest.raises(FileExistsError, match=SCHEMA_BOUNDARIES):
        seg_main.process_custom_segmentation(_boxes(tmp_path), out, cell_id_col="my_cell")
    assert (out / SCHEMA_BOUNDARIES).stat().st_mtime_ns == t0
    os.utime(out / SCHEMA_BOUNDARIES, ns=(t0 - 10**10, t0 - 10**10))  # coarse clock: make "newer" unambiguous
    t0 = (out / SCHEMA_BOUNDARIES).stat().st_mtime_ns
    seg_main.process_custom_segmentation(_boxes(tmp_path), out, cell_id_col="my_cell", overwrite=True)
    again = gpd.read_parquet(out / SCHEMA_BOUNDARIES)
    assert (out / SCHEMA_BOUNDARIES).stat().st_mtime_ns > t0
    pd.testing.assert_frame_equal(pd.DataFrame(first.drop(columns="Geometry")),
                                  pd.DataFrame(again.drop(columns="Geometry")))


# ---- align-segmentation ----

def test_align_refused_before_aligning(ds, monkeypatch):
    import spida.S.segmentation.align as align
    calls = []

    def fake_align(**k):
        calls.append("align")
        p = Path(k["output_dir"]) / f"aligned_{k['prefix1']}_{k['prefix2']}_polygons.parquet"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")
        return {"output_path": str(p), "n_aligned": 1}

    monkeypatch.setattr(align, "align_segmentations", fake_align)
    monkeypatch.setattr(seg_main, "process_custom_segmentation", lambda *a, **k: calls.append(("custom", k["overwrite"])))
    seg_main.align_geometries(ds)
    with pytest.raises(FileExistsError, match="aligned_cellpose_nuc_cellpose_cell_polygons.parquet"):
        seg_main.align_geometries(ds)
    seg_main.align_geometries(ds, overwrite=True)
    assert calls == ["align", ("custom", False), "align", ("custom", True)]


# ---- CLI: --overwrite is not forwarded as a backend parameter ----

def test_segment_cli_overwrite_not_a_backend_param(monkeypatch, tmp_path):
    from spida.S import cli as s_cli
    seen = {}
    monkeypatch.setattr(seg_main, "segment_region", lambda method, d, **k: seen.update(k))
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"PROCESSED_ROOT_PATH": str(tmp_path / "raw"), "SEGMENTATION_OUT_PATH": str(tmp_path / "seg")}))
    res = CliRunner().invoke(s_cli.cli, ["--config", str(cfg), "segment-region", "proseg", "--exp", "E",
                                         "--reg", "region_1", "--overwrite", "--voxel-size", "2"])
    assert res.exit_code == 0, res.output
    assert seen["overwrite"] is True and "voxel-size" in str(seen) and list(seen).count("overwrite") == 1


# ---- drift: what a step writes is what its guard lists (phase 3d) ----

def test_files_written_per_method():
    from spida.io.schema.segmentation import SCHEMA_FILES, files_written
    from spida.S.segmentation.backends import get_spec
    assert set(files_written(get_spec("cellpose"))) == set(SCHEMA_FILES)       # partitions transcripts
    assert set(files_written(get_spec("proseg"))) == {SCHEMA_BOUNDARIES, SCHEMA_SUM_SIGNALS, SCHEMA_CELL_METADATA}
    assert files_written(get_spec("cellpose"), "vpt")[0] == "cellpose_micron_space.parquet"
    assert set(files_written()) == set(SCHEMA_FILES)


def test_custom_writes_nothing_outside_its_guard(tmp_path):
    from spida.io.schema.segmentation import files_written
    tx = tmp_path / "tx.csv"
    pd.DataFrame([{"global_x": i * 20 + 7, "global_y": 7, "global_z": 0, "gene": f"g{i % 2}", "barcode_id": i % 2}
                  for i in range(3)]).to_csv(tx, index=False)
    out = tmp_path / "out"
    seg_main.process_custom_segmentation(_boxes(tmp_path), out, cell_id_col="my_cell", transcripts_path=tx)
    written = {p.name for p in out.iterdir() if p.is_file()}
    assert written and written <= set(files_written())  # a re-run would be refused for every file it wrote


# ---- inputs checked before any work (phase 3e) ----

def test_process_missing_images_stops_before_ingest(ds, backends, post):
    seg_main.segment_region("cellpose", ds, prefix_name="cellpose_cell")
    (ds.raw / "images" / "mosaic_DAPI_z0.tif").unlink()
    (ds.raw / "images" / "micron_to_mosaic_pixel_transform.csv").unlink()     # images dir now empty
    with pytest.raises(FileNotFoundError, match="images"):
        seg_main.process_segmentation_region("cellpose", ds, prefix_name="cellpose_cell")
    assert post == []  # nothing ran


def test_load_missing_schema_file_stops_before_loading(ds, monkeypatch):
    from spida.S.load import main as io_main
    called = []
    monkeypatch.setattr(io_main, "read_segmentation", lambda *a, **k: called.append(1))
    (ds.seg / "cellpose_cell").mkdir(parents=True)
    (ds.zarr / "shapes" / "default_E_region_1_polygons").mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="cell_by_gene.csv"):
        io_main.load_segmentation_region(ds, type="cellpose", prefix_name="cellpose_cell")
    assert not called
