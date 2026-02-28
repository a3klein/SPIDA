import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon

from spida.S.io.filt_to_qc import apply_qc_filter_to_segmentation


class _FakeSData(dict):
    def __init__(self):
        super().__init__()
        self.tables = {}


def _square(x0: float, y0: float, size: float = 1.0) -> Polygon:
    return Polygon(
        [(x0, y0), (x0 + size, y0), (x0 + size, y0 + size), (x0, y0 + size)]
    )


@pytest.fixture
def base_sdata():
    sdata = _FakeSData()

    shapes = gpd.GeoDataFrame(
        {
            "EntityID": ["c1", "c2"],
            "geometry": [_square(0, 0, 1.0), _square(10, 10, 1.0)],
        },
        index=["c1", "c2"],
        geometry="geometry",
    )

    table_obs = pd.DataFrame(
        {
            "EntityID": ["c1", "c2"],
            "other": [1, 2],
        },
        index=["c1", "c2"],
    )
    table = ad.AnnData(
        X=np.array([[1, 0], [0, 1]]),
        obs=table_obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    table.uns["spatialdata_attrs"] = {"instance_key": "EntityID"}

    points = pd.DataFrame(
        {
            "x": [0.1, 10.2, 5.0],
            "y": [0.1, 10.2, 5.0],
            "gene": ["A", "B", "C"],
            "cell_id": ["c1", "c2", "-1"],
        }
    )

    sdata["seg_shapes"] = shapes
    sdata["seg_table"] = table
    sdata["seg_points"] = points
    return sdata


def test_apply_qc_filter_to_segmentation_with_qc_regions(base_sdata):
    qc_regions = gpd.GeoDataFrame(
        {"geometry": [_square(-0.5, -0.5, 3.0)]},
        geometry="geometry",
    )

    out = apply_qc_filter_to_segmentation(
        base_sdata,
        table_key="seg_table",
        shapes_key="seg_shapes",
        points_key="seg_points",
        qc_regions=qc_regions,
    )

    # Keep only c1 shape and c1 table row
    assert list(out["seg_shapes"].index) == ["c1"]
    assert list(out["seg_table"].obs_names) == ["c1"]

    # Points from filtered-out cell (c2) should be unassigned
    pts = out["seg_points"]
    assert pts.loc[0, "cell_id"] == "c1"
    assert pts.loc[1, "cell_id"] == "-1"
    assert pts.loc[2, "cell_id"] == "-1"


def test_apply_qc_filter_uses_qc_table_key_with_geometry_wkt(base_sdata):
    qc_obs = pd.DataFrame(
        {
            "geometry_wkt": [
                _square(-0.5, -0.5, 3.0).wkt,
                _square(20.0, 20.0, 2.0).wkt,
            ],
            # convention in transcript_qc: filtered=False means passed QC
            "filtered": [False, True],
        },
        index=["h1", "h2"],
    )
    qc_adata = ad.AnnData(
        X=np.zeros((2, 1)),
        obs=qc_obs,
        var=pd.DataFrame(index=["dummy"]),
    )
    base_sdata.tables["qc_hex"] = qc_adata

    out = apply_qc_filter_to_segmentation(
        base_sdata,
        table_key="seg_table",
        shapes_key="seg_shapes",
        points_key="seg_points",
        qc_table_key="qc_hex",
        qc_filter_col="filtered",
        qc_pass_value=False,
    )

    assert list(out["seg_shapes"].index) == ["c1"]
    assert list(out["seg_table"].obs_names) == ["c1"]
    assert out["seg_points"].loc[1, "cell_id"] == "-1"


def test_apply_qc_filter_numeric_points_assignment(base_sdata):
    # Switch points to numeric assignment column
    base_sdata["seg_points"] = pd.DataFrame(
        {
            "x": [0.1, 10.2, 5.0],
            "y": [0.1, 10.2, 5.0],
            "gene": ["A", "B", "C"],
            "assignment": [1, 2, -1],
        }
    )
    # Match numeric ids in table/shapes
    base_sdata["seg_shapes"].index = pd.Index(["1", "2"])
    base_sdata["seg_shapes"]["EntityID"] = ["1", "2"]
    base_sdata["seg_table"].obs.index = pd.Index(["1", "2"])
    base_sdata["seg_table"].obs["EntityID"] = ["1", "2"]

    qc_regions = gpd.GeoDataFrame(
        {"geometry": [_square(-0.5, -0.5, 3.0)]},
        geometry="geometry",
    )

    out = apply_qc_filter_to_segmentation(
        base_sdata,
        table_key="seg_table",
        shapes_key="seg_shapes",
        points_key="seg_points",
        qc_regions=qc_regions,
    )

    pts = out["seg_points"]
    assert pts.loc[0, "assignment"] == 1
    assert pts.loc[1, "assignment"] == -1
    assert pts.loc[2, "assignment"] == -1


def test_apply_qc_filter_raises_for_missing_points_cell_column(base_sdata):
    base_sdata["seg_points"] = pd.DataFrame(
        {
            "x": [0.1],
            "y": [0.1],
            "gene": ["A"],
            "not_a_cell_col": ["c1"],
        }
    )

    qc_regions = gpd.GeoDataFrame(
        {"geometry": [_square(-0.5, -0.5, 3.0)]},
        geometry="geometry",
    )

    with pytest.raises(KeyError, match="Could not infer points cell-id column"):
        apply_qc_filter_to_segmentation(
            base_sdata,
            table_key="seg_table",
            shapes_key="seg_shapes",
            points_key="seg_points",
            qc_regions=qc_regions,
        )
