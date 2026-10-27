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
    res, seen = _run(monkeypatch, "spida.S.load.main.ingest_region", cfg, "ingest-region", "--exp", "E", "--reg", "region_1")
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert d.zarr == tmp_path / "z" / "E" / "region_1"
    assert d.raw == tmp_path / "raw" / "E" / "out" / "region_1"
    assert d.figures == tmp_path / "img" / "E" / "region_1"


def test_ingest_paths_only_is_nameless(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.load.main.ingest_region", cfg, "ingest-region",
                     "--zarr_path", str(tmp_path / "x.zarr"), "--raw_path", str(tmp_path / "r"))
    assert res.exit_code == 0, res.output
    d = seen["args"][0]
    assert not d.named and d.zarr == tmp_path / "x.zarr" and d.raw == tmp_path / "r"


def test_names_with_path_override(monkeypatch, cfg, tmp_path):
    res, seen = _run(monkeypatch, "spida.S.load.main.load_segmentation_region", cfg, "load-segmentation-region",
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
    res, _ = _run(monkeypatch, "spida.S.load.main.ingest_region", cfg, "ingest-region", *args)
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

    from spida.io.load import _drop_namespace

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
    (tmp_path / "raw" / "images").mkdir(parents=True)                     # inputs align requires
    (tmp_path / "raw" / "images" / "mosaic_DAPI_z0.tif").write_text("")
    (tmp_path / "raw" / "detected_transcripts.csv").write_text("")
    for p in ("nuc", "cell"):
        (tmp_path / "z" / "shapes" / f"{p}_E_R_polygons").mkdir(parents=True)
    m.align_geometries(d, prefix1="nuc", prefix2="cell")
    assert seen["align"]["shape_key1"] == "nuc_E_R_polygons" and seen["align"]["output_dir"] == tmp_path / "seg" / "align"
    b, o, k = seen["custom"]
    assert o == tmp_path / "seg" / "align" and k["boundaries_space"] == "pixel"
    assert k["images_dir"] == tmp_path / "raw" / "images"


# ---- ingest --overwrite (phase 2) ----

def _store(path):
    path.mkdir(parents=True)
    (path / "zarr.json").write_text("{}")


def test_ingest_region_refuses_existing_store(monkeypatch, tmp_path):
    from spida.S.load import main as m
    from spida.io.layout import Dataset
    called = []
    monkeypatch.setattr(m, "read_merscope", lambda *a, **k: called.append(1))
    d = Dataset(zarr=tmp_path / "x.zarr", raw=tmp_path / "r")
    (d.raw).mkdir()
    (d.raw / "detected_transcripts.csv").write_text("")
    _store(d.zarr)
    with pytest.raises(FileExistsError, match="x.zarr"):
        m.ingest_region(d)
    assert not called  # stops before reading the raw data


def test_ingest_region_overwrite_proceeds(monkeypatch, tmp_path):
    from spida.S.load import main as m
    from spida.io.layout import Dataset
    monkeypatch.setattr(m, "read_merscope", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("reached")))
    d = Dataset(zarr=tmp_path / "x.zarr", raw=tmp_path / "r")
    (d.raw).mkdir()
    (d.raw / "detected_transcripts.csv").write_text("")
    _store(d.zarr)
    with pytest.raises(RuntimeError, match="reached"):
        m.ingest_region(d, overwrite=True)


def test_ingest_all_skips_existing(monkeypatch, cfg, tmp_path):
    from spida.S.load import main as m
    for r in ("region_1", "region_2"):
        (tmp_path / "raw" / "E" / "out" / r).mkdir(parents=True)
    _store(tmp_path / "z" / "E" / "region_1")
    seen = []
    monkeypatch.setattr(m, "ingest_region", lambda d, **k: seen.append((d.reg, k["overwrite"])))
    res = CliRunner().invoke(cli_module.cli, ["--config", str(cfg), "ingest-all", "--exp", "E"])
    assert res.exit_code == 0, res.output
    assert seen == [("region_2", False)]
    res = CliRunner().invoke(cli_module.cli, ["--config", str(cfg), "ingest-all", "--exp", "E", "--overwrite"])
    assert sorted(seen[1:]) == [("region_1", True), ("region_2", True)]


# ---- load-segmentation --overwrite (phase 2) ----

def _seg_store(path):
    import anndata as ad
    import dask.dataframe as dd
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import spatialdata as sd
    from shapely.geometry import Point
    from spatialdata.models import PointsModel, ShapesModel, TableModel
    from spida._constants import POINTS_KEY, SHAPES_KEY, TABLE_KEY
    from spida.io.layout import element_keys

    keys = element_keys("seg", None)
    shapes = gpd.GeoDataFrame({"EntityID": [0, 1]}, geometry=[Point(0, 0).buffer(1), Point(5, 5).buffer(1)], index=["0", "1"])
    pts = dd.from_pandas(pd.DataFrame({"x": [0.0, 5.0], "y": [0.0, 5.0], "gene": ["A", "B"]}), npartitions=1)
    obs = pd.DataFrame({"cells_region": pd.Categorical([keys[SHAPES_KEY]] * 2), "EntityID": [0, 1], "v": [1, 1]}, index=["0", "1"])
    table = TableModel.parse(ad.AnnData(np.ones((2, 2), dtype=np.float32), obs=obs),
                             region=keys[SHAPES_KEY], region_key="cells_region", instance_key="EntityID")
    s = sd.SpatialData(shapes={keys[SHAPES_KEY]: ShapesModel.parse(shapes)},
                       points={keys[POINTS_KEY]: PointsModel.parse(pts)}, tables={keys[TABLE_KEY]: table})
    s.attrs["spida"] = {"namespace": "", "exp": None, "reg": None}  # as ingest-region records it
    s.write(path)
    return keys


def _reload_with_v2(path, keys):
    import numpy as np
    import spatialdata as sd
    from spatialdata.transformations import Affine, set_transformation
    from spida._constants import SHAPES_KEY, TABLE_KEY

    s = sd.read_zarr(path)
    t = s[keys[TABLE_KEY]].copy()
    t.obs["v"] = 2
    s[keys[TABLE_KEY]] = t
    set_transformation(s[keys[SHAPES_KEY]], Affine(np.diag([2.0, 2.0, 1.0]), ("x", "y"), ("x", "y")), "global")
    return s


def test_finalize_rerun_without_overwrite_raises(tmp_path):
    import spatialdata as sd
    from spida.io.load import write_segmentation
    from spida._constants import TABLE_KEY

    path = tmp_path / "s.zarr"
    keys = _seg_store(path)
    with pytest.raises(FileExistsError):
        write_segmentation(_reload_with_v2(path, keys), keys)
    assert (sd.read_zarr(path)[keys[TABLE_KEY]].obs["v"] == 1).all()


def test_finalize_rerun_with_overwrite_replaces_and_keeps_transforms(tmp_path):
    import numpy as np
    import spatialdata as sd
    from spatialdata.transformations import get_transformation
    from spida.io.load import write_segmentation
    from spida._constants import POINTS_KEY, SHAPES_KEY, TABLE_KEY

    path = tmp_path / "s.zarr"
    keys = _seg_store(path)
    out = write_segmentation(_reload_with_v2(path, keys), keys, overwrite=True)
    r = sd.read_zarr(path)
    assert (r[keys[TABLE_KEY]].obs["v"] == 2).all() and (out[keys[TABLE_KEY]].obs["v"] == 2).all()
    for k in (POINTS_KEY, SHAPES_KEY):
        assert "pixel" in get_transformation(r[keys[k]], get_all=True)
    g = get_transformation(r[keys[SHAPES_KEY]], to_coordinate_system="global")
    assert np.allclose(g.to_affine_matrix(("x", "y"), ("x", "y")), np.diag([2.0, 2.0, 1.0]))
    assert not any("__tmp" in p for p in r.elements_paths_on_disk())


def test_load_segmentation_region_refuses_before_loading(monkeypatch, tmp_path):
    from spida.S.load import main as m
    from spida.io.layout import Dataset

    d = Dataset(zarr=tmp_path / "s.zarr", seg=tmp_path / "seg")
    _seg_store(d.zarr)
    (d.seg / "seg").mkdir(parents=True)
    for f in ("boundaries_micron.parquet", "cell_by_gene.csv", "cell_metadata.csv", "detected_transcripts.csv"):
        (d.seg / "seg" / f).write_text("")                     # the schema files load-segmentation requires
    (d.zarr / "shapes" / "default_polygons").mkdir(parents=True)
    called = []
    monkeypatch.setattr(m, "read_segmentation", lambda sdata, *a, **k: called.append("read") or sdata)
    monkeypatch.setattr(m, "write_segmentation", lambda sdata, keys, overwrite=False: called.append(overwrite))
    with pytest.raises(FileExistsError, match="seg_table"):
        m.load_segmentation_region(d, type="cellpose", prefix_name="seg")
    assert not called
    m.load_segmentation_region(d, type="cellpose", prefix_name="seg", overwrite=True)
    assert called == ["read", True]


def _qc_load(monkeypatch, tmp_path, filt, mask=True):
    """load_segmentation_region with stubbed read / filter / write; returns the call log."""
    from spida.S.load import main as m
    from spida.io.layout import Dataset

    d = Dataset(zarr=tmp_path / "s.zarr", seg=tmp_path / "seg")
    _seg_store(d.zarr)
    (d.seg / "seg").mkdir(parents=True)
    for f in ("boundaries_micron.parquet", "cell_by_gene.csv", "cell_metadata.csv", "detected_transcripts.csv"):
        (d.seg / "seg" / f).write_text("")
    (d.zarr / "shapes" / "default_polygons").mkdir(parents=True)
    if mask:
        (d.zarr / "shapes" / "qc_mask").mkdir(parents=True)
    log = []
    monkeypatch.setattr(m, "read_segmentation", lambda sdata, *a, **k: log.append(("read", k)) or sdata)
    monkeypatch.setattr(m, "apply_qc_filter_to_segmentation", lambda sdata, **k: filt(log, k) or sdata)
    monkeypatch.setattr(m, "write_segmentation", lambda sdata, keys, overwrite=False: log.append(("write", overwrite)))
    m.load_segmentation_region(d, type="cellpose", prefix_name="seg", overwrite=True,
                               qc_shapes_key="qc_mask", qc_pass_value=True)
    return log


def test_load_segmentation_qc_filter_between_read_and_write(monkeypatch, tmp_path):
    log = _qc_load(monkeypatch, tmp_path, lambda log, k: log.append(("filter", k)))
    assert [s for s, _ in log] == ["read", "filter", "write"]
    assert "qc_shapes_key" not in log[0][1]                      # the reader never sees the QC options
    assert log[1][1]["qc_shapes_key"] == "qc_mask" and log[1][1]["qc_pass_value"] is True
    assert log[1][1]["shapes_key"] == "seg_polygons"


def test_load_segmentation_qc_filter_failure_writes_nothing(monkeypatch, tmp_path):
    def filt(log, k):
        raise ValueError("QC filter removed all cells")
    with pytest.raises(ValueError, match="removed all"):
        _qc_load(monkeypatch, tmp_path, filt)


def test_load_segmentation_missing_qc_mask_stops_before_reading(monkeypatch, tmp_path):
    with pytest.raises(FileNotFoundError, match="qc_mask"):
        _qc_load(monkeypatch, tmp_path, lambda log, k: pytest.fail("filtered"), mask=False)


def test_load_segmentation_all_skips_loaded(monkeypatch, cfg, tmp_path):
    from spida.S.load import main as m
    import spida.io.layout as layout

    for r in ("region_1", "region_2"):
        (tmp_path / "seg" / "E" / r / "seg").mkdir(parents=True)
    keys = layout.element_keys("seg", "E_region_1")
    from spida._constants import TABLE_KEY
    (tmp_path / "z" / "E" / "region_1" / "tables" / keys[TABLE_KEY]).mkdir(parents=True)
    seen = []
    monkeypatch.setattr(m, "load_segmentation_region", lambda d, **k: seen.append((d.reg, k["overwrite"])))
    res = CliRunner().invoke(cli_module.cli, ["--config", str(cfg), "load-segmentation-all", "--exp", "E",
                                              "--type", "cellpose", "--prefix_name", "seg"])
    assert res.exit_code == 0, res.output
    assert seen == [("region_2", False)]


# ---- load-deconvolution --overwrite (phase 2) ----

def test_load_decon_refuses_then_swaps(monkeypatch, tmp_path):
    import numpy as np
    import spatialdata as sd
    from spatialdata.models import Image2DModel
    import spida.io.load as ie
    from spida.S.load import main as m
    from spida.io.layout import Dataset

    d = Dataset(zarr=tmp_path / "s.zarr", raw=tmp_path / "r")
    _seg_store(d.zarr)
    img = lambda v: Image2DModel.parse(np.full((1, 64, 64), v, np.uint8), dims=("c", "y", "x"), scale_factors=[2], chunks=32)
    monkeypatch.setattr(ie, "_read_decon_images", lambda *a, **k: img(1))
    m.load_deconvolution_region(d, image_dir=tmp_path)
    called = []
    monkeypatch.setattr(ie, "_read_decon_images", lambda *a, **k: called.append(1) or img(2))
    with pytest.raises(FileExistsError, match="images/decon_image"):
        m.load_deconvolution_region(d, image_dir=tmp_path)
    assert not called  # stops before reading the images
    m.load_deconvolution_region(d, image_dir=tmp_path, overwrite=True)
    r = sd.read_zarr(d.zarr)
    assert int(r["decon_image"]["scale0"]["image"].values.min()) == 2 and "decon_image__tmp" not in r.images
