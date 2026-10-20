"""transcript-qc / cluster-hexes on a real (synthetic) zarr store: guards, overwrite, and the shared hex table.

No fakes: real overlay_hexes, transcript_qc, cluster_hexes and spatialdata writes.
"""
import os
import json
import warnings

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")
sd = pytest.importorskip("spatialdata")

import dask.dataframe as dd
from click.testing import CliRunner
from spatialdata.models import PointsModel

from spida.io.layout import Dataset, element_keys
from spida._constants import POINTS_KEY

PTS = element_keys("default", None)[POINTS_KEY]
HEX = "adata_hex_s30.0_o0.0"
CLUSTERED = f"{HEX}_clustered"
CL_KW = dict(min_cells=1, min_genes=5, n_top_genes=10, pca_comps=5, leiden_resolution=0.5)


@pytest.fixture
def ds(tmp_path):
    # 600 x 300 um; density rises left to right, so the transcript threshold splits the hexes
    rng = np.random.default_rng(0)
    n = 60000
    x = 600 * np.sqrt(rng.uniform(0, 1, n))
    pts = pd.DataFrame({"x": x, "y": rng.uniform(0, 300, n),
                        "gene": pd.Categorical(rng.choice([f"g{i}" for i in range(20)], n))})
    s = sd.SpatialData(points={PTS: PointsModel.parse(dd.from_pandas(pts, npartitions=2), feature_key="gene")})
    s.attrs["spida"] = {"namespace": "", "exp": None, "reg": None}
    d = Dataset(zarr=tmp_path / "s.zarr", figures=tmp_path / "fig")
    s.write(d.zarr)
    return d


def _mtime(d, rel):
    return max(f.stat().st_mtime_ns for f in (d.zarr / rel).rglob("*") if f.is_file())


def _qc(d, **k):
    from spida.P.main import transcript_qc_region
    transcript_qc_region(d, hex_size=30.0, hex_overlap=0.0, **k)


def _cluster(d, **k):
    from spida.P.main import cluster_hexes_region
    cluster_hexes_region(d, hex_size=30.0, hex_overlap=0.0, **{**CL_KW, **k})


def _filtered(d):
    return int(sd.read_zarr(d.zarr)["transcript_qc_shapes"]["filtered"].sum())


def test_first_run_saves_qc_shapes_and_hex_table(ds):
    _qc(ds, min_transcripts=100)
    r = sd.read_zarr(ds.zarr)
    g = r["transcript_qc_shapes"]
    assert HEX in r.tables
    assert set(g.geom_type) == {"Polygon"} and len(g) == r[HEX].n_obs
    assert 0 < g["filtered"].sum() < len(g)  # the threshold splits the region


def test_rerun_without_overwrite_raises_before_compute(ds, monkeypatch):
    import spida.P.transcript_qc as tq
    _qc(ds, min_transcripts=100)
    before = {rel: _mtime(ds, rel) for rel in ("shapes/transcript_qc_shapes", f"tables/{HEX}")}
    monkeypatch.setattr(tq, "get_or_create_hex_adata", lambda *a, **k: pytest.fail("computed before the guard"))
    with pytest.raises(FileExistsError, match="shapes/transcript_qc_shapes"):
        _qc(ds, min_transcripts=300)
    assert {rel: _mtime(ds, rel) for rel in before} == before


def test_rerun_with_overwrite_replaces_shapes_reuses_hex_table(ds):
    _qc(ds, min_transcripts=100)
    n100, hex_t = _filtered(ds), _mtime(ds, f"tables/{HEX}")
    _qc(ds, min_transcripts=300, overwrite=True)
    assert _filtered(ds) > n100  # stricter threshold -> more hexes flagged
    assert _mtime(ds, f"tables/{HEX}") == hex_t  # shared hex table not rewritten
    assert not any("__tmp" in p for p in sd.read_zarr(ds.zarr).elements_paths_on_disk())


def test_overwrite_result_equals_fresh_run(ds, tmp_path):
    # overwrite path gives exactly what a first run with the same settings gives
    _qc(ds, min_transcripts=100)
    _qc(ds, min_transcripts=300, overwrite=True)
    swapped = sd.read_zarr(ds.zarr)["transcript_qc_shapes"]
    from spida.P.main import transcript_qc_region
    import shutil
    fresh = Dataset(zarr=tmp_path / "fresh.zarr", figures=tmp_path / "fig2")
    shutil.copytree(ds.zarr, fresh.zarr)
    s = sd.read_zarr(fresh.zarr)
    s.delete_element_from_disk("transcript_qc_shapes")
    transcript_qc_region(fresh, hex_size=30.0, hex_overlap=0.0, min_transcripts=300)
    first = sd.read_zarr(fresh.zarr)["transcript_qc_shapes"]
    pd.testing.assert_frame_equal(pd.DataFrame(swapped.drop(columns="geometry")),
                                  pd.DataFrame(first.drop(columns="geometry")))
    assert swapped.geometry.geom_equals(first.geometry).all()


def test_cluster_after_qc_not_blocked_and_reuses_hex_table(ds):
    _qc(ds, min_transcripts=100)
    hex_t, qc_t = _mtime(ds, f"tables/{HEX}"), _mtime(ds, "shapes/transcript_qc_shapes")
    _cluster(ds, min_transcripts=100)
    r = sd.read_zarr(ds.zarr)
    assert CLUSTERED in r.tables and "leiden" in r[CLUSTERED].obs
    assert _mtime(ds, f"tables/{HEX}") == hex_t and _mtime(ds, "shapes/transcript_qc_shapes") == qc_t


def test_qc_after_cluster_not_blocked(ds):
    _cluster(ds, min_transcripts=100)  # creates the hex table
    hex_t, cl_t = _mtime(ds, f"tables/{HEX}"), _mtime(ds, f"tables/{CLUSTERED}")
    _qc(ds, min_transcripts=100)
    assert "transcript_qc_shapes" in sd.read_zarr(ds.zarr).shapes
    assert _mtime(ds, f"tables/{HEX}") == hex_t and _mtime(ds, f"tables/{CLUSTERED}") == cl_t


def test_cluster_rerun_guard_and_overwrite(ds):
    _cluster(ds, min_transcripts=100)
    t0 = _mtime(ds, f"tables/{CLUSTERED}")
    with pytest.raises(FileExistsError, match=CLUSTERED):
        _cluster(ds, min_transcripts=100)
    assert _mtime(ds, f"tables/{CLUSTERED}") == t0
    for f in (ds.zarr / "tables" / CLUSTERED).rglob("*"):  # coarse clock: make "newer" unambiguous
        if f.is_file():
            os.utime(f, ns=(t0 - 10**10, t0 - 10**10))
    t0 = _mtime(ds, f"tables/{CLUSTERED}")
    _cluster(ds, min_transcripts=300, overwrite=True)
    r = sd.read_zarr(ds.zarr)
    assert _mtime(ds, f"tables/{CLUSTERED}") > t0
    assert r[CLUSTERED].n_obs < r[HEX].n_obs  # stricter threshold kept fewer hexes
    assert not any("__tmp" in p for p in r.elements_paths_on_disk())


def test_persist_off_writes_nothing(ds):
    import spida.P.transcript_qc as tq
    s = sd.read_zarr(ds.zarr)
    tq.run_transcript_qc(s, points_key=PTS, hex_size=30.0, hex_overlap=0.0, min_transcripts=100, persist_qc_shapes=False)
    assert "transcript_qc_shapes" not in sd.read_zarr(ds.zarr).shapes


def test_cli_flag(ds, tmp_path):
    from spida.P import cli as p_cli
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"IMAGE_STORE_PATH": str(tmp_path / "img")}))
    args = ["--config", str(cfg), "transcript-qc", "--zarr_path", str(ds.zarr), "--hex_size", "30",
            "--min_transcripts", "100"]
    assert CliRunner().invoke(p_cli.cli, args).exit_code == 0
    res = CliRunner().invoke(p_cli.cli, args)
    assert res.exit_code != 0 and isinstance(res.exception, FileExistsError)
    res = CliRunner().invoke(p_cli.cli, args + ["--overwrite"])
    assert res.exit_code == 0, res.output
