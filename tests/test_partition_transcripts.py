"""Tests for the pure-Python VPT `partition-transcripts` reimplementation.

The fixture under tests/data/partition_transcripts/ was produced by running the
real `vpt partition-transcripts` binary on a small spatially-trimmed region
(216 cells, ~13k transcripts). These tests assert bit-for-bit equivalence.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.segmentation_utils import partition_transcripts

DATA = Path(__file__).parent / "data" / "segmentation_utils"
BOUNDARIES = DATA / "cellpose_micron_space.parquet"
TRANSCRIPTS = DATA / "detected_transcripts.parquet"
EXPECTED_CBG = DATA / "expected_cell_by_gene.parquet"
EXPECTED_CELL_ID = DATA / "expected_cell_id.parquet"

pytestmark = pytest.mark.skipif(
    not BOUNDARIES.exists(), reason="partition-transcripts fixture not present"
)


# --------------------------------------------------------------------------- #
# equivalence vs the real VPT binary (committed fixture)                      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("chunk_size", [None, 2000, 5000])
def test_cell_by_gene_matches_vpt(chunk_size):
    expected = pd.read_parquet(EXPECTED_CBG)
    got = partition_transcripts(BOUNDARIES, TRANSCRIPTS, chunk_size=chunk_size)

    # identical values, row order, and column order
    assert list(got.index) == list(expected.index)
    assert list(got.columns) == list(expected.columns)
    assert got.index.name == "cell"
    pd.testing.assert_frame_equal(got, expected, check_dtype=True)


@pytest.mark.parametrize("chunk_size", [None, 2000])
def test_cell_id_matches_vpt(chunk_size, tmp_path):
    out_tx = tmp_path / "out.csv"
    partition_transcripts(
        BOUNDARIES, TRANSCRIPTS,
        output_transcripts=out_tx, chunk_size=chunk_size,
    )
    mine = pd.read_csv(out_tx)
    expected = pd.read_parquet(EXPECTED_CELL_ID)  # tx_uid, cell_id

    merged = mine[["tx_uid", "cell_id"]].merge(
        expected, on="tx_uid", suffixes=("_mine", "_vpt")
    )
    assert len(merged) == len(expected)
    assert (merged["cell_id_mine"] == merged["cell_id_vpt"]).all()


def test_cell_by_gene_row_sums_equal_assigned_counts():
    """Every assigned transcript is counted exactly once per (cell, gene)."""
    cbg = partition_transcripts(BOUNDARIES, TRANSCRIPTS)
    expected = pd.read_parquet(EXPECTED_CBG)
    assert int(cbg.values.sum()) == int(expected.values.sum())


def test_output_files_written(tmp_path):
    cbg_path = tmp_path / "cell_by_gene.csv"
    partition_transcripts(BOUNDARIES, TRANSCRIPTS, output_entity_by_gene=cbg_path)
    assert cbg_path.exists()
    on_disk = pd.read_csv(cbg_path, index_col=0)
    assert on_disk.index.name == "cell"


# --------------------------------------------------------------------------- #
# synthetic edge cases                                                        #
# --------------------------------------------------------------------------- #
def _write_synthetic(tmp_path):
    """Two unit-square cells at z=0; one cell only at z=1. Known-answer setup."""
    polys = gpd.GeoDataFrame(
        {
            "EntityID": [10, 20, 30],
            "ZIndex": [0, 0, 1],
            "geometry": [
                geom.box(0, 0, 1, 1),      # cell 10, z0
                geom.box(2, 2, 3, 3),      # cell 20, z0
                geom.box(0, 0, 1, 1),      # cell 30, z1
            ],
        },
        geometry="geometry",
    )
    bpath = tmp_path / "bounds.parquet"
    polys.to_parquet(bpath)

    tx = pd.DataFrame(
        {
            "tx_uid": [0, 1, 2, 3, 4, 5],
            "barcode_id": [1, 1, 2, 2, 1, 3],
            "gene": ["A", "A", "B", "B", "A", "C"],
            "global_x": [0.5, 2.5, 0.5, 0.5, 5.0, 0.5],
            "global_y": [0.5, 2.5, 0.5, 0.5, 5.0, 0.5],
            "global_z": [0, 0, 0, 1, 0, 9],  # last: out-of-range z -> unassigned
        }
    )
    tpath = tmp_path / "tx.parquet"
    tx.to_parquet(tpath, index=False)
    return bpath, tpath


def test_synthetic_assignment(tmp_path):
    bpath, tpath = _write_synthetic(tmp_path)
    out_tx = tmp_path / "out.csv"
    cbg = partition_transcripts(bpath, tpath, output_transcripts=out_tx)

    # rows = all three cells, sorted; columns ordered by barcode_id: A(1),B(2),C(3)
    assert list(cbg.index) == [10, 20, 30]
    assert list(cbg.columns) == ["A", "B", "C"]

    # tx0 (A,z0) -> cell 10 ; tx1 (A,z0) -> cell 20 ; tx2 (B,z0) inside cell 10's
    # box at z0 -> cell 10 ; tx3 (B,z1) inside cell 30's box at z1 -> cell 30
    expected_counts = {
        (10, "A"): 1, (10, "B"): 1, (10, "C"): 0,
        (20, "A"): 1, (20, "B"): 0,
        (30, "B"): 1,
    }
    got_counts = {k: int(cbg.loc[k[0], k[1]]) for k in expected_counts}
    assert got_counts == expected_counts

    # tx4 (5,5) is inside no cell; tx5 has out-of-range z -> both unassigned
    mine = pd.read_csv(out_tx).set_index("tx_uid")["cell_id"].to_dict()
    assert mine == {0: 10, 1: 20, 2: 10, 3: 30, 4: -1, 5: -1}


def test_boundary_point_excluded(tmp_path):
    """A point exactly on the polygon boundary is NOT contained (GEOS within)."""
    polys = gpd.GeoDataFrame(
        {"EntityID": [1], "ZIndex": [0], "geometry": [geom.box(0, 0, 1, 1)]},
        geometry="geometry",
    )
    bpath = tmp_path / "b.parquet"
    polys.to_parquet(bpath)
    tx = pd.DataFrame(
        {
            "tx_uid": [0, 1],
            "barcode_id": [1, 1],
            "gene": ["A", "A"],
            "global_x": [0.0, 0.5],  # on-edge vs interior
            "global_y": [0.5, 0.5],
            "global_z": [0, 0],
        }
    )
    tpath = tmp_path / "t.parquet"
    tx.to_parquet(tpath, index=False)
    cbg = partition_transcripts(bpath, tpath)
    assert cbg["A"].to_dict() == {1: 1}  # only the interior point counts
