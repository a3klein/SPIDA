import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from spida.S import cli as cli_module


@pytest.fixture
def cfg(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({k: str(tmp_path / v) for k, v in {
        "ZARR_STORAGE_PATH": "z", "PROCESSED_ROOT_PATH": "raw",
        "SEGMENTATION_OUT_PATH": "seg", "IMAGE_STORE_PATH": "img"}.items()}))
    return p


def _run(monkeypatch, target, cfg, *args):
    seen = {}
    monkeypatch.setattr(target, lambda *a, **k: seen.update(args=a, kwargs=k))
    res = CliRunner().invoke(cli_module.cli, ["--config", str(cfg), *args])
    return res, seen


def test_ingest_names_resolve_the_tree(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.io.main.ingest_region", cfg, "ingest-region", "--exp", "E", "--reg", "region_1")
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert d.zarr == tmp_path / "z" / "E" / "region_1"
    assert d.raw == tmp_path / "raw" / "E" / "out" / "region_1"
    assert d.figures == tmp_path / "img" / "E" / "region_1"


def test_ingest_paths_only_is_nameless(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.io.main.ingest_region", cfg, "ingest-region",
                     "--zarr_path", str(tmp_path / "x.zarr"), "--raw_path", str(tmp_path / "r"))
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert not d.named and d.zarr == tmp_path / "x.zarr" and d.raw == tmp_path / "r"


def test_names_with_path_override(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.io.main.load_segmentation_region", cfg, "load-segmentation-region",
                     "--exp", "E", "--reg", "region_1", "--zarr_path", str(tmp_path / "x.zarr"), "--prefix_name", "cellpose_cell")
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert d.zarr == tmp_path / "x.zarr" and d.seg == tmp_path / "seg" / "E" / "region_1"


def test_segment_region_path_mode(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.segmentation.main.segment_region", cfg, "segment-region", "proseg",
                     "--raw_path", str(tmp_path / "r"), "--seg_path", str(tmp_path / "s"))
    assert res.exit_code == 0, res.output
    method, d = seen["args"]
    assert method == "proseg" and d.raw == tmp_path / "r" and d.seg == tmp_path / "s"


@pytest.mark.parametrize("args", [["--exp", "E"], []])
def test_bad_dataset_arguments_are_usage_errors(monkeypatch, cfg, args):
    res, _ = _run(monkeypatch, "spida.S.io.main.ingest_region", cfg, "ingest-region", *args)
    assert res.exit_code == 2


def test_drop_namespace_renames_and_repoints_table():
    import warnings
    warnings.filterwarnings("ignore")
    import spatialdata as sd  # before dask.dataframe
    import anndata as ad
    import dask.dataframe as dd
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    from shapely.geometry import Point
    from spatialdata.models import PointsModel, ShapesModel, TableModel

    from spida.S.io.ingest_exp import _drop_namespace

    ds = "default_region"
    shapes = ShapesModel.parse(gpd.GeoDataFrame({"EntityID": [0, 1]}, geometry=[Point(0, 0).buffer(1), Point(3, 3).buffer(1)], index=["0", "1"]))
    pts = PointsModel.parse(dd.from_pandas(pd.DataFrame({"x": [0.0], "y": [0.0], "gene": ["A"]}), npartitions=1))
    obs = pd.DataFrame({"cells_region": pd.Categorical([f"{ds}_polygons"] * 2), "EntityID": ["0", "1"]}, index=["0", "1"])
    table = TableModel.parse(ad.AnnData(np.ones((2, 1), dtype=np.float32), obs=obs),
                             region=f"{ds}_polygons", region_key="cells_region", instance_key="EntityID")
    s = sd.SpatialData(shapes={f"{ds}_polygons": shapes}, points={f"{ds}_transcripts": pts}, tables={"table": table})

    _drop_namespace(s, dataset_id=ds, prefix_name="default")
    assert set(s.shapes) == {"default_polygons"} and set(s.points) == {"default_transcripts"}
    assert s["table"].uns["spatialdata_attrs"]["region"] == "default_polygons"
    assert set(s["table"].obs["cells_region"]) == {"default_polygons"}


def test_align_segmentation_uses_dataset(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.segmentation.main.align_geometries", cfg, "align-segmentation",
                     "--exp", "E", "--reg", "region_1", "--prefix1", "nuc", "--prefix2", "cell")
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert d.seg_dir("align") == tmp_path / "seg" / "E" / "region_1" / "align"
    assert seen["kwargs"]["prefix1"] == "nuc" and seen["kwargs"]["out_dir_name"] == "align"


def test_align_geometries_goes_native(monkeypatch, tmp_path):
    """align -> process_custom_segmentation (no VPT), pixel-space boundaries, outputs in seg/{label}."""
    from spida.io.layout import Dataset
    from spida.S.segmentation import main as m
    import spida.S.segmentation.align as align

    seen = {}
    monkeypatch.setattr(align, "align_segmentations", lambda **k: seen.update(align=k) or
                        {"output_path": str(k["output_dir"] / "aligned.parquet"), "n_aligned": 1})
    monkeypatch.setattr(m, "process_custom_segmentation", lambda b, o, **k: seen.update(custom=(b, o, k)))
    d = Dataset(exp="E", reg="R", zarr=tmp_path / "z", raw=tmp_path / "raw", seg=tmp_path / "seg")
    m.align_geometries(d, prefix1="nuc", prefix2="cell")
    assert seen["align"]["shape_key1"] == "nuc_E_R_polygons" and seen["align"]["output_dir"] == tmp_path / "seg" / "align"
    b, o, k = seen["custom"]
    assert o == tmp_path / "seg" / "align" and k["boundaries_space"] == "pixel"
    assert k["images_dir"] == tmp_path / "raw" / "images"
