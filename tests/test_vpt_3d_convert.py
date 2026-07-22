"""Stage 4A regression: after run_cellpose stops stamping ZIndex/ZLevel/EntityID,
the vpt-backend 3D path (_3d_convert_geometry) must derive them itself and
reproduce the release behavior (same ZIndex/ZLevel/geometry/cell-structure).
"""
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

from spida.S.segmentation.vpt_utils import _3d_convert_geometry

_SPACING = 1.5


def _sq(x, y, s=4):
    return Polygon([(x, y), (x, y + s), (x + s, y + s), (x + s, y)])


def _raw_3d(n_cells=5, n_z=7):
    """RAW cellpose 3D output: only ID, z, Geometry (post-Stage-4)."""
    rows = []
    for cid in range(1, n_cells + 1):
        for z in range(n_z):
            rows.append({"ID": cid, "z": z, "geometry": _sq(cid * 10, z * 10)})
    return gpd.GeoDataFrame(rows, geometry="geometry")


def _write(tmp_path, gdf, region="region_T"):
    (tmp_path / region).mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(tmp_path / region / "polygons.parquet", index=False)


def test_3d_convert_derives_schema_from_raw(tmp_path):
    _write(tmp_path, _raw_3d(n_cells=5, n_z=7))
    _3d_convert_geometry(str(tmp_path), str(tmp_path), "region_T",
                         convert_micron=False, spacing_z=_SPACING)
    out = gpd.read_parquet(tmp_path / "region_T" / "cellpose_micron_space.parquet")

    # ZIndex == z; ZLevel == (z+1)*spacing (1-based)
    assert (out["ZIndex"] == out["z"]).all()
    assert np.allclose(out["ZLevel"], (out["z"] + 1) * _SPACING)
    # one EntityID per cellpose ID, consistent across the cell's z-planes
    assert out.groupby("ID")["EntityID"].nunique().eq(1).all()
    assert out["EntityID"].nunique() == out["ID"].nunique() == 5
    # schema bookkeeping cols present
    assert (out["Type"] == "cell").all()
    for c in ("ParentID", "ParentType", "Name"):
        assert c in out.columns
    # geometry passes through unchanged when convert_micron=False
    assert out.geometry.area.round(3).eq(16.0).all()


def _apply_min_z_on_ID(g, min_z):
    id_vc = g.groupby("ID")["z"].count()
    return g[g["ID"].isin(id_vc[id_vc >= min_z].index)].copy()


def _apply_min_z_on_EntityID(g, min_z):     # the OLD (pre-4A) logic
    g = g.assign(EntityID=pd.factorize(g["ID"])[0] + 1)
    eid_vc = g.groupby("EntityID")["z"].count()
    return g[g["EntityID"].isin(eid_vc[eid_vc >= min_z].index)].copy()


def test_min_z_filter_on_ID_equals_EntityID(tmp_path):
    # run_cellpose now filters min-z on ID instead of EntityID=factorize(ID)+1;
    # both must keep the same cells (bijection).
    g = _raw_3d(n_cells=6, n_z=7)
    g = g[~((g.ID == 3) & (g.z >= 2))]            # cell 3 present in only 2 planes
    g = g[~((g.ID == 5) & (g.z >= 1))]            # cell 5 present in only 1 plane
    keep_id = set(_apply_min_z_on_ID(g, 3)["ID"].unique())
    keep_eid = set(_apply_min_z_on_EntityID(g, 3)["ID"].unique())
    assert keep_id == keep_eid == {1, 2, 4, 6}


def test_min_z_filter_equivalence_adversarial_ids():
    # non-contiguous, unsorted IDs -> factorize order (appearance) differs from
    # sorted ID order. Also uneven per-ID z-plane counts. The ID- and EntityID-
    # based filters must still keep exactly the same rows for every min_z.
    rows = []
    plane_counts = {22: 7, 3: 5, 10: 1, 7: 3, 99: 2, 4: 6}   # ID -> # z-planes present
    for cid, npl in plane_counts.items():
        for z in range(npl):
            rows.append({"ID": cid, "z": z, "geometry": _sq(cid, z)})
    g = gpd.GeoDataFrame(rows, geometry="geometry")
    assert list(pd.factorize(g["ID"])[0][:1]) == [0]          # sanity: appearance-ordered
    for min_z in (1, 2, 3, 5, 7):
        a = _apply_min_z_on_ID(g, min_z).sort_values(["ID", "z"]).reset_index(drop=True)
        b = _apply_min_z_on_EntityID(g, min_z).sort_values(["ID", "z"]).reset_index(drop=True)
        # same rows kept (ID+z), for every threshold
        assert a[["ID", "z"]].equals(b[["ID", "z"]]), f"min_z={min_z}"
    # spot-check the expected keep set at min_z=3
    assert set(_apply_min_z_on_ID(g, 3)["ID"].unique()) == {22, 3, 7, 4}


def test_3d_convert_reproduces_release_zlevels(tmp_path):
    # Feed the OLD-style stamped input; the new derive-from-(ID,z) must reproduce
    # the release's ZIndex/ZLevel exactly (EntityID values are timestamp-based and
    # non-deterministic in both, so only structure is compared).
    raw = _raw_3d(n_cells=4, n_z=7)
    old = raw.copy()
    old["ZIndex"] = old["z"]
    old["ZLevel"] = (old["z"] + 1) * _SPACING
    old["EntityID"] = pd.factorize(old["ID"])[0] + 1
    _write(tmp_path, raw)
    _3d_convert_geometry(str(tmp_path), str(tmp_path), "region_T",
                         convert_micron=False, spacing_z=_SPACING)
    out = gpd.read_parquet(tmp_path / "region_T" / "cellpose_micron_space.parquet")
    m = out.merge(old[["ID", "z", "ZIndex", "ZLevel"]], on=["ID", "z"], suffixes=("", "_old"))
    assert (m["ZIndex"] == m["ZIndex_old"]).all()
    assert np.allclose(m["ZLevel"], m["ZLevel_old"])
