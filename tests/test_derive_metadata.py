"""Tests for the pure-Python VPT `derive-entity-metadata` reimplementation.

Ground truth (expected_cell_metadata.parquet) was produced by the real
`vpt derive-entity-metadata` binary on the committed fixture. All columns are
asserted equal EXCEPT ``anisotropy``, which is shapely-version dependent (see
segmentation_utils.derive_entity_metadata docstring).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.segmentation_utils import derive_entity_metadata

DATA = Path(__file__).parent / "data" / "segmentation_utils"
BOUNDARIES = DATA / "cellpose_micron_space.parquet"
CELL_BY_GENE = DATA / "expected_cell_by_gene.parquet"  # read + written as csv below
EXPECTED_META = DATA / "expected_cell_metadata.parquet"

# columns compared for equality (anisotropy excluded on purpose)
_COMPARED = [
    "fov", "volume", "center_x", "center_y", "min_x", "min_y", "max_x", "max_y",
    "transcript_count", "perimeter_area_ratio", "solidity",
]

pytestmark = pytest.mark.skipif(
    not BOUNDARIES.exists(), reason="segmentation_utils fixture not present"
)


def _cell_by_gene_csv(tmp_path):
    """The fixture ships cell-by-gene as parquet; VPT/derive reads CSV."""
    csv = tmp_path / "cell_by_gene.csv"
    pd.read_parquet(CELL_BY_GENE).to_csv(csv)
    return csv


def test_metadata_matches_vpt(tmp_path):
    expected = pd.read_parquet(EXPECTED_META)
    got = derive_entity_metadata(BOUNDARIES, _cell_by_gene_csv(tmp_path))

    assert list(got.index) == list(expected.index)
    assert list(got.columns) == list(expected.columns)
    assert got.index.name == "EntityID"

    got_a = got.reindex(index=expected.index)
    for col in _COMPARED:
        np.testing.assert_allclose(
            got_a[col].to_numpy(float), expected[col].to_numpy(float),
            rtol=1e-6, atol=1e-6, equal_nan=True, err_msg=f"column {col}",
        )


def test_anisotropy_present_and_finite():
    """anisotropy is not asserted against VPT, but must be computed & sane."""
    got = derive_entity_metadata(BOUNDARIES)
    aniso = got["anisotropy"].to_numpy(float)
    finite = aniso[~np.isnan(aniso)]
    assert finite.size > 0
    assert (finite >= 1.0 - 1e-9).all()  # major/minor ratio is always >= 1


def test_transcript_count_optional():
    """Without cell_by_gene, transcript_count is NaN; other columns still computed."""
    got = derive_entity_metadata(BOUNDARIES)
    assert got["transcript_count"].isna().all()
    assert got["volume"].notna().any()


def test_output_file_written(tmp_path):
    out = tmp_path / "cell_metadata.csv"
    derive_entity_metadata(BOUNDARIES, _cell_by_gene_csv(tmp_path), output_metadata=out)
    assert out.exists()
    on_disk = pd.read_csv(out, index_col=0)
    assert on_disk.index.name == "EntityID"
    assert list(on_disk.columns) == [
        "fov", "volume", "center_x", "center_y", "min_x", "min_y", "max_x",
        "max_y", "anisotropy", "transcript_count", "perimeter_area_ratio",
        "solidity",
    ]


# --------------------------------------------------------------------------- #
# synthetic known-answer test                                                 #
# --------------------------------------------------------------------------- #
def test_synthetic_volume_and_center(tmp_path):
    """One cell, a unit square extruded over 2 z-planes of 1.5 um thickness."""
    polys = gpd.GeoDataFrame(
        {
            "EntityID": [7, 7],
            "ZIndex": [0, 1],
            "ZLevel": [1.5, 3.0],
            "geometry": [geom.box(0, 0, 2, 1), geom.box(0, 0, 2, 1)],
        },
        geometry="geometry",
    )
    bpath = tmp_path / "b.parquet"
    polys.to_parquet(bpath)
    meta = derive_entity_metadata(bpath)

    row = meta.loc[7]
    # area per plane = 2; thickness per plane = 1.5; two planes -> volume = 6
    assert row["volume"] == pytest.approx(2 * 1.5 + 2 * 1.5)
    assert row["center_x"] == pytest.approx(1.0)
    assert row["center_y"] == pytest.approx(0.5)
    assert row["min_x"] == pytest.approx(0.0)
    assert row["max_x"] == pytest.approx(2.0)
    assert row["anisotropy"] == pytest.approx(2.0)  # 2x1 rectangle
    assert row["solidity"] == pytest.approx(1.0)     # convex shape
