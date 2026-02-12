from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity
from shapely.geometry import Polygon

from spida.S.segmentation.align import align_segmentations
from spida.S.segmentation import vpt as vpt_module


def _square(x0: float, y0: float, size: float = 1.0) -> Polygon:
    return Polygon(
        [(x0, y0), (x0 + size, y0), (x0 + size, y0 + size), (x0, y0 + size)]
    )


@pytest.fixture
def segmentation_fixture(tmp_path: Path):
    exp_name = "test_experiment"
    reg_name = "region_test"
    prefix1 = "cellpose_nuc"
    prefix2 = "cellpose_cell"

    zarr_path = tmp_path / "zarr_store"
    seg_out_dir = tmp_path / "seg_out"
    output_dir = seg_out_dir / reg_name
    root_dir = tmp_path / "root"

    shape_key1 = f"{prefix1}_{exp_name}_{reg_name}_polygons"
    shape_key2 = f"{prefix2}_{exp_name}_{reg_name}_polygons"

    shapes1 = gpd.GeoDataFrame(
        {
            "EntityID": ["c1", "c2", "c3"],
            "geometry": [_square(0, 0), _square(2, 0), _square(4, 0)],
        },
        index=["c1", "c2", "c3"],
        crs="EPSG:4326",
    )
    shapes2 = gpd.GeoDataFrame(
        {
            "EntityID": ["d1", "d2", "d3"],
            "geometry": [_square(0.2, 0.2), _square(2.2, 0.2), _square(4.2, 0.2)],
        },
        index=["d1", "d2", "d3"],
        crs="EPSG:4326",
    )

    shapes1 = ShapesModel.parse(shapes1, transformations={"global": Identity()})
    shapes2 = ShapesModel.parse(shapes2, transformations={"global": Identity()})

    sdata = sd.SpatialData(shapes={shape_key1: shapes1, shape_key2: shapes2})
    sdata.write(zarr_path)

    region_root = root_dir / reg_name
    (region_root / "images").mkdir(parents=True, exist_ok=True)
    (region_root / "images" / "micron_to_mosaic_pixel_transform.csv").write_text(
        "x,y\n0,0\n"
    )
    pd.DataFrame({"gene": ["G1"], "x": [0.0], "y": [0.0]}).to_csv(
        region_root / "detected_transcripts.csv", index=False
    )

    return {
        "exp_name": exp_name,
        "reg_name": reg_name,
        "prefix1": prefix1,
        "prefix2": prefix2,
        "output_dir": output_dir,
        "seg_out_dir": seg_out_dir,
        "zarr_path": zarr_path,
        "root_dir": root_dir,
    }


def test_align_segmentations_geometry_modes(segmentation_fixture):
    exp_name = segmentation_fixture["exp_name"]
    reg_name = segmentation_fixture["reg_name"]
    prefix1 = segmentation_fixture["prefix1"]
    prefix2 = segmentation_fixture["prefix2"]
    output_dir = segmentation_fixture["output_dir"]
    zarr_path = segmentation_fixture["zarr_path"]
    geometry_modes = ['larger', 'prefix1', 'prefix2', 'intersection']

    for mode in geometry_modes:
        result = align_segmentations(
            zarr_path=zarr_path,
            exp_name=exp_name,
            reg_name=reg_name,
            prefix1=prefix1,
            prefix2=prefix2,
            geometry_mode=mode,
            output_dir=output_dir
        )
        assert result['geometry_mode'] == mode


def test_align_segmentations_file_formats(segmentation_fixture):
    exp_name = segmentation_fixture["exp_name"]
    reg_name = segmentation_fixture["reg_name"]
    prefix1 = segmentation_fixture["prefix1"]
    prefix2 = segmentation_fixture["prefix2"]
    output_dir = segmentation_fixture["output_dir"]
    zarr_path = segmentation_fixture["zarr_path"]
    # Test for valid file formats
    result = align_segmentations(
        zarr_path=zarr_path,
        exp_name=exp_name,
        reg_name=reg_name,
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir=output_dir
    )
    # Check if output files are created in the expected format
    output_path = Path(result["output_path"])
    assert output_path.exists()
    assert output_path.suffix == ".parquet"


def test_no_cells_missing_or_introduced(segmentation_fixture):
    exp_name = segmentation_fixture["exp_name"]
    reg_name = segmentation_fixture["reg_name"]
    prefix1 = segmentation_fixture["prefix1"]
    prefix2 = segmentation_fixture["prefix2"]
    output_dir = segmentation_fixture["output_dir"]
    zarr_path = segmentation_fixture["zarr_path"]
    # Run the alignment
    result = align_segmentations(
        zarr_path=zarr_path,
        exp_name=exp_name,
        reg_name=reg_name,
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir=output_dir
    )
    aligned_cells = pd.read_parquet(result["output_path"])
    assert result["n_aligned"] == 3
    assert len(aligned_cells) == 3


def test_seg_to_vpt_on_aligned_output(segmentation_fixture, monkeypatch):
    exp_name = segmentation_fixture["exp_name"]
    reg_name = segmentation_fixture["reg_name"]
    prefix1 = segmentation_fixture["prefix1"]
    prefix2 = segmentation_fixture["prefix2"]
    output_dir = segmentation_fixture["output_dir"]
    seg_out_dir = segmentation_fixture["seg_out_dir"]
    zarr_path = segmentation_fixture["zarr_path"]
    root_dir = segmentation_fixture["root_dir"]

    result = align_segmentations(
        zarr_path=zarr_path,
        exp_name=exp_name,
        reg_name=reg_name,
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir=output_dir,
    )

    aligned_name = Path(result["output_path"]).name

    calls: dict[str, dict[str, object]] = {}

    def _record(name):
        def _inner(*args, **kwargs):
            calls[name] = {"args": args, "kwargs": kwargs}
            return None

        return _inner

    monkeypatch.setattr(vpt_module, "_add_vpt_binary", lambda *args, **kwargs: None)
    monkeypatch.setattr(vpt_module, "_convert_geometry", _record("convert"))
    monkeypatch.setattr(vpt_module, "_cli_partition_transcripts", _record("partition"))
    monkeypatch.setattr(vpt_module, "_cli_get_metadata", _record("metadata"))

    vpt_module.seg_to_vpt(
        root_dir=str(root_dir),
        seg_out_dir=str(seg_out_dir),
        region=reg_name,
        input_boundaries=aligned_name,
    )

    assert "convert" in calls
    assert "partition" in calls
    assert "metadata" in calls
    assert calls["convert"]["kwargs"]["input_boundaries"] == aligned_name
    assert (
        calls["partition"]["kwargs"]["input_boundaries"]
        == calls["convert"]["kwargs"]["output_boundaries"]
    )
    assert (
        calls["metadata"]["kwargs"]["input_boundaries"]
        == calls["convert"]["kwargs"]["output_boundaries"]
    )
