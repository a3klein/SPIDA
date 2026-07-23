"""Tests for process_custom_segmentation (segmentation/main.py) + ingest_custom_polygons.

Covers the user-provided-segmentation path on synthetic inputs: 2D micron boundaries
expanded to cylinders across the transcript planes, pixel->micron transform, 3D input
with strict plane matching (valid + mismatch), user-metadata merge, and boundaries-only
single-plane output. Asserts the segmentation schema, EntityID assignment + caller-id
preservation, and step auto-selection. sum-signals (image) path is exercised by the
segmentation_utils tests; here we cover the boundary/transcript/metadata composition.
"""

import gzip
import shutil

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.S.segmentation.main import process_custom_segmentation

_SCHEMA = ["ID", "Geometry", "EntityID", "ParentID", "ParentType", "Type", "ZIndex", "ZLevel"]


def _boxes_2d(tmp_path, n=3):
    """n non-overlapping 15x15 boxes, one row per cell, keyed by 'my_cell'."""
    rows = [{"my_cell": f"c{i}", "Geometry": geom.box(i * 20, 0, i * 20 + 15, 15)}
            for i in range(n)]
    g = gpd.GeoDataFrame(rows, geometry="Geometry")
    p = tmp_path / "bounds.parquet"
    g.to_parquet(p)
    return p


def _transcripts(tmp_path, n=3, planes=(0, 1, 2)):
    """One transcript per (cell, plane) placed inside each box. barcode_id<->gene is
    1:1 (as in real MERFISH); cell i uses gene g{i%2} / barcode i%2."""
    rows = []
    for i in range(n):
        for z in planes:
            rows.append({"global_x": i * 20 + 7, "global_y": 7, "global_z": z,
                         "gene": f"g{i % 2}", "barcode_id": i % 2})
    p = tmp_path / "tx.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _boxes_3d(tmp_path, n=2, planes=(0, 1, 2)):
    rows = [{"my_cell": f"c{i}", "zp": z, "Geometry": geom.box(i * 20, 0, i * 20 + 15, 15)}
            for i in range(n) for z in planes]
    g = gpd.GeoDataFrame(rows, geometry="Geometry")
    p = tmp_path / "bounds3d.parquet"
    g.to_parquet(p)
    return p


def test_2d_micron_with_transcripts(tmp_path):
    out = tmp_path / "out"
    process_custom_segmentation(_boxes_2d(tmp_path), out, cell_id_col="my_cell",
                                transcripts_path=_transcripts(tmp_path))
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    assert all(c in bnd.columns for c in _SCHEMA)
    assert "my_cell" in bnd.columns                              # caller id preserved
    assert sorted(bnd["ZIndex"].unique()) == [0, 1, 2]           # cylinder across tx planes
    assert np.allclose(sorted(bnd["ZLevel"].unique()), [1.5, 3.0, 4.5])  # 1-based
    assert bnd["EntityID"].nunique() == 3
    assert bnd.groupby("my_cell")["EntityID"].nunique().eq(1).all()   # mapping inherent
    cbg = pd.read_csv(out / "cell_by_gene.csv", index_col=0)
    assert cbg.values.sum() == 9                                 # 3 cells x 3 planes
    assert (out / "detected_transcripts.csv").exists()
    meta = pd.read_csv(out / "cell_metadata.csv", index_col=0)
    assert len(meta) == 3 and "volume" in meta.columns


def test_pixel_space_applies_transform(tmp_path):
    m2m = tmp_path / "m2m.csv"
    m2m.write_text("2.0 0.0 0.0\n0.0 2.0 0.0\n0.0 0.0 1.0\n")   # scale x2 micron->pixel
    out = tmp_path / "out"
    process_custom_segmentation(_boxes_2d(tmp_path, n=1), out, cell_id_col="my_cell",
                                boundaries_space="pixel", micron_to_mosaic_path=m2m)
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    # inverse (halve) of the pixel box(0,0,15,15) -> micron area (15/2)^2
    assert np.isclose(bnd.geometry.area.iloc[0], (15 / 2) ** 2)


def test_3d_valid(tmp_path):
    out = tmp_path / "o3"
    process_custom_segmentation(_boxes_3d(tmp_path), out, cell_id_col="my_cell",
                                boundary_z_col="zp", n_z_planes=3)
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    assert sorted(bnd["ZIndex"].unique()) == [0, 1, 2]
    assert bnd["EntityID"].nunique() == 2
    assert bnd.groupby("my_cell")["EntityID"].nunique().eq(1).all()


def test_3d_boundary_plane_mismatch_raises(tmp_path):
    b = _boxes_3d(tmp_path, planes=(0, 1))                       # only 2 planes
    with pytest.raises(ValueError):
        process_custom_segmentation(b, tmp_path / "o", cell_id_col="my_cell",
                                    boundary_z_col="zp", n_z_planes=3)


def test_3d_transcript_plane_mismatch_raises(tmp_path):
    b = _boxes_3d(tmp_path, planes=(0, 1, 2))
    t = _transcripts(tmp_path, n=2, planes=(0, 1))               # tx only 2 planes
    with pytest.raises(ValueError):
        process_custom_segmentation(b, tmp_path / "o", cell_id_col="my_cell",
                                    boundary_z_col="zp", n_z_planes=3, transcripts_path=t)


def test_metadata_merge(tmp_path):
    md = tmp_path / "md.csv"
    pd.DataFrame({"my_cell": ["c0", "c1", "c2"], "cell_type": ["A", "B", "A"]}).to_csv(md, index=False)
    out = tmp_path / "out"
    process_custom_segmentation(_boxes_2d(tmp_path), out, cell_id_col="my_cell",
                                transcripts_path=_transcripts(tmp_path), metadata_path=md)
    meta = pd.read_csv(out / "cell_metadata.csv", index_col=0)
    assert "cell_type" in meta.columns and "volume" in meta.columns   # user + derived merged
    assert set(meta["cell_type"].dropna()) == {"A", "B"}


def test_2d_boundaries_only_single_plane(tmp_path):
    out = tmp_path / "o"
    process_custom_segmentation(_boxes_2d(tmp_path), out, cell_id_col="my_cell")
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    assert sorted(bnd["ZIndex"].unique()) == [0]                 # no tx/images -> single 2D plane
    assert not (out / "cell_by_gene.csv").exists()               # no transcripts -> no partition
    assert len(pd.read_csv(out / "cell_metadata.csv", index_col=0)) == 3


def _boxes_geojson(tmp_path, name="bounds.geojson", n=3):
    rows = [{"my_cell": f"c{i}", "geometry": geom.box(i * 20, 0, i * 20 + 15, 15)}
            for i in range(n)]
    p = tmp_path / name
    gpd.GeoDataFrame(rows, geometry="geometry").to_file(p, driver="GeoJSON")
    return p


def test_geojson_boundary_input(tmp_path):
    # the .geojson branch of _read_boundaries (real proseg format) -> same schema output
    out = tmp_path / "out"
    process_custom_segmentation(_boxes_geojson(tmp_path), out, cell_id_col="my_cell")
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    assert all(c in bnd.columns for c in _SCHEMA)
    assert "my_cell" in bnd.columns
    assert bnd["EntityID"].nunique() == 3
    assert bnd.groupby("my_cell")["EntityID"].nunique().eq(1).all()


def test_geojson_gz_boundary_input(tmp_path):
    # the .geojson.gz branch (gzip -> read_file)
    plain = _boxes_geojson(tmp_path, name="b.geojson")
    gzp = tmp_path / "b.geojson.gz"
    with open(plain, "rb") as fi, gzip.open(gzp, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    out = tmp_path / "out"
    process_custom_segmentation(gzp, out, cell_id_col="my_cell")
    assert gpd.read_parquet(out / "boundaries_micron.parquet")["EntityID"].nunique() == 3


def test_images_without_transform_raises(tmp_path):
    with pytest.raises(ValueError, match="micron_to_mosaic"):
        process_custom_segmentation(_boxes_2d(tmp_path), tmp_path / "o",
                                    cell_id_col="my_cell", images_dir=str(tmp_path))


def test_metadata_without_cell_id_raises(tmp_path):
    md = tmp_path / "md.csv"
    pd.DataFrame({"my_cell": ["c0"], "cell_type": ["A"]}).to_csv(md, index=False)
    with pytest.raises(ValueError, match="cell_id_col"):
        process_custom_segmentation(_boxes_2d(tmp_path), tmp_path / "o",
                                    transcripts_path=_transcripts(tmp_path), metadata_path=md)


def test_pixel_without_transform_raises(tmp_path):
    with pytest.raises(ValueError, match="pixel"):
        process_custom_segmentation(_boxes_2d(tmp_path), tmp_path / "o",
                                    cell_id_col="my_cell", boundaries_space="pixel")


def test_3d_without_boundary_z_col_raises(tmp_path):
    with pytest.raises(ValueError, match="boundary_z_col"):
        process_custom_segmentation(_boxes_3d(tmp_path), tmp_path / "o",
                                    cell_id_col="my_cell", n_z_planes=3)


def test_2d_with_boundary_z_col_raises(tmp_path):
    with pytest.raises(ValueError, match="boundary_z_col"):
        process_custom_segmentation(_boxes_3d(tmp_path), tmp_path / "o",
                                    cell_id_col="my_cell", boundary_z_col="zp", n_z_planes=1)


def test_transcript_z_in_microns_converted(tmp_path):
    # transcript z given as micron ZLevel = (plane+1)*z_spacing -> converts to planes 0,1,2
    b = _boxes_2d(tmp_path)
    rows = []
    for i in range(3):
        for zi in (0, 1, 2):
            rows.append({"global_x": i * 20 + 7, "global_y": 7,
                         "global_z": (zi + 1) * 1.5, "gene": "g0", "barcode_id": 0})
    t = tmp_path / "tx_um.csv"
    pd.DataFrame(rows).to_csv(t, index=False)
    out = tmp_path / "out"
    process_custom_segmentation(b, out, cell_id_col="my_cell", transcripts_path=t,
                                transcript_z_in_microns=True)
    bnd = gpd.read_parquet(out / "boundaries_micron.parquet")
    assert sorted(bnd["ZIndex"].unique()) == [0, 1, 2]
    assert pd.read_csv(out / "cell_by_gene.csv", index_col=0).values.sum() == 9
