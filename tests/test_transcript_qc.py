import numpy as np
import pandas as pd
import anndata as ad
import geopandas as gpd
import shapely.geometry as geom
import pytest

from spida.P import transcript_qc as tq
from spida.utilities.tiling import hex_area_from_size


class _FakeDaskDF:
    def __init__(self, df):
        self._df = df

    def compute(self):
        return self._df


class _FakeSData:
    def __init__(self, points):
        self._points = points
        self.tables = {}
        self.written = []

    def __getitem__(self, key):
        if key == "points":
            return self._points
        raise KeyError(key)

    def write_element(self, key):
        self.written.append(key)


@pytest.fixture
def points_df():
    # points strictly inside [0, 0, 10, 10]
    return pd.DataFrame(
        {
            "x": [2.0, 3.0, 4.0, 5.0, 6.0],
            "y": [2.0, 3.0, 4.0, 5.0, 6.0],
            "gene": ["A", "A", "B", "Blank-1", "C"],
        }
    )


def test_ensure_geodataframe_from_dataframe(points_df):
    gdf = tq._ensure_geodataframe(points_df, x_col="x", y_col="y")
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert "geometry" in gdf.columns
    assert len(gdf) == len(points_df)


def test_ensure_geodataframe_from_dask_like(points_df):
    dask_like = _FakeDaskDF(points_df)
    gdf = tq._ensure_geodataframe(dask_like, x_col="x", y_col="y")
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert len(gdf) == len(points_df)


def test_hex_table_key():
    assert tq._hex_table_key(30, 0) == "adata_hex_s30_o0"
    assert tq._hex_table_key(30.5, 0.1) == "adata_hex_s30.5_o0.1"


def test_add_hex_metadata():
    adata = ad.AnnData(X=np.array([[1, 2], [3, 4]]), obs=pd.DataFrame(index=["0", "1"]), var=pd.DataFrame(index=["g1", "g2"]))
    tq._add_hex_metadata(adata, hex_size=10, hex_overlap=0.2)
    assert "hexgrid" in adata.uns
    assert adata.uns["hexgrid"]["hex_size"] == 10
    assert adata.uns["hexgrid"]["hex_overlap"] == 0.2
    assert np.isclose(adata.uns["hexgrid"]["hex_area_um2"], hex_area_from_size(10))


def test_apply_hex_filter_density_supersedes_min_transcripts():
    obs = pd.DataFrame(
        {
            "tz_count": [100, 10],
            "density": [0.1, 10.0],
        },
        index=["h1", "h2"],
    )
    adata = ad.AnnData(X=np.array([[1, 0], [0, 1]]), obs=obs, var=pd.DataFrame(index=["g1", "g2"]))

    # If min_density is set, it should be used even if min_transcripts would keep h1
    out = tq._apply_hex_filter(adata, min_transcripts=50, min_density=1.0)
    assert list(out.obs.index) == ["h2"]


def test_overlay_hexes_basic(monkeypatch, points_df):
    # Patch grid creation to deterministic single polygon containing all points
    poly = geom.Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

    def _fake_grid(bounds, hex_size, overlap):
        return gpd.GeoDataFrame({"geometry": [poly]}, geometry="geometry")

    monkeypatch.setattr(tq, "create_hexagonal_grid", _fake_grid)

    adata_hex, grid = tq.overlay_hexes(
        points_df,
        x_col="x",
        y_col="y",
        hex_size=5,
        hex_overlap=0.0,
        gene_col="gene",
        return_grid=True,
    )

    # Blank-* gene should be removed
    assert "Blank-1" not in adata_hex.var_names
    assert "tz_count" in adata_hex.obs.columns
    assert "n_genes" in adata_hex.obs.columns
    assert "density" in adata_hex.obs.columns
    assert "hexgrid" in adata_hex.uns
    assert len(grid) == 1


def test_get_or_create_hex_adata_caches(monkeypatch, points_df):
    sdata = _FakeSData(points_df)
    calls = {"n": 0}

    def _fake_overlay(*args, **kwargs):
        calls["n"] += 1
        obs = pd.DataFrame(
            {
                "geometry": [geom.Point(1, 1)],
                "tz_count": [3],
                "n_genes": [2],
                "density": [0.5],
            },
            index=["0"],
        )
        return ad.AnnData(X=np.array([[1, 2]]), obs=obs, var=pd.DataFrame(index=["A", "B"]))

    monkeypatch.setattr(tq, "overlay_hexes", _fake_overlay)

    a1, key1 = tq.get_or_create_hex_adata(sdata, points_key="points", hex_size=30, hex_overlap=0)
    a2, key2 = tq.get_or_create_hex_adata(sdata, points_key="points", hex_size=30, hex_overlap=0)

    assert key1 == key2 == "adata_hex_s30_o0"
    assert calls["n"] == 1
    assert "adata_hex_s30_o0" in sdata.tables
    assert "adata_hex_s30_o0" in sdata.written
    assert a1 is a2


def test_transcript_qc_requires_hex_area_if_density_missing():
    obs = pd.DataFrame({"tz_count": [10]}, index=["h1"])
    adata = ad.AnnData(X=np.array([[1, 0]]), obs=obs, var=pd.DataFrame(index=["g1", "g2"]))

    with pytest.raises(ValueError, match="hex_area_um2"):
        tq.transcript_qc(adata, min_density=0.1)


def test_run_transcript_qc(monkeypatch):
    sdata = _FakeSData(points=pd.DataFrame({"x": [1], "y": [1], "gene": ["A"]}))

    obs = pd.DataFrame(
        {
            "geometry": [geom.Point(1, 1), geom.Point(2, 2)],
            "tz_count": [100, 5],
            "n_genes": [10, 2],
            "density": [2.0, 0.1],
        },
        index=["h1", "h2"],
    )
    adata_hex = ad.AnnData(X=np.array([[1, 2], [0, 1]]), obs=obs, var=pd.DataFrame(index=["A", "B"]))
    adata_hex.uns["hexgrid"] = {"hex_size": 30, "hex_overlap": 0, "hex_area_um2": hex_area_from_size(30)}

    def _fake_get_or_create(*args, **kwargs):
        return adata_hex, "adata_hex_s30_o0"

    monkeypatch.setattr(tq, "get_or_create_hex_adata", _fake_get_or_create)

    adata_qc, grid, key = tq.run_transcript_qc(
        sdata=sdata,
        points_key="points",
        min_transcripts=50,
    )

    assert key == "adata_hex_s30_o0"
    assert list(adata_qc.obs.index) == ["h1"]
    assert "filtered" in grid.columns
    assert bool(grid.loc["h1", "filtered"]) is False
    assert bool(grid.loc["h2", "filtered"]) is True


def test_run_cluster_hexes_saves_clustered_table(monkeypatch):
    sdata = _FakeSData(points=pd.DataFrame({"x": [1], "y": [1], "gene": ["A"]}))

    obs = pd.DataFrame(
        {
            "geometry": [geom.Point(1, 1), geom.Point(2, 2)],
            "tz_count": [100, 5],
            "n_genes": [10, 2],
            "density": [2.0, 0.1],
        },
        index=["h1", "h2"],
    )
    adata_hex = ad.AnnData(X=np.array([[1, 2], [0, 1]]), obs=obs, var=pd.DataFrame(index=["A", "B"]))
    adata_hex.uns["hexgrid"] = {"hex_size": 30, "hex_overlap": 0, "hex_area_um2": hex_area_from_size(30)}

    def _fake_get_or_create(*args, **kwargs):
        return adata_hex, "adata_hex_s30_o0"

    def _fake_cluster(adata_hex, **kwargs):
        out = adata_hex.copy()
        out.obs["leiden"] = pd.Categorical(["0"] * out.n_obs)
        return out

    monkeypatch.setattr(tq, "get_or_create_hex_adata", _fake_get_or_create)
    monkeypatch.setattr(tq, "cluster_hexes", _fake_cluster)

    clustered, cluster_key = tq.run_cluster_hexes(
        sdata=sdata,
        points_key="points",
        min_density=1.0,  # keeps only h1
    )

    assert cluster_key == "adata_hex_s30_o0_clustered"
    assert cluster_key in sdata.tables
    assert cluster_key in sdata.written
    assert "leiden" in clustered.obs.columns
    assert list(clustered.obs.index) == ["h1"]