"""spida.io.anndata: SPIDA's AnnData files (h5ad or zarr-encoded directory), obs-only reads, temp-swap writes."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from spida.io.anndata import read_adata, read_obs, write_adata


def _adata(n=20, seed=0):
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame({"label": pd.Categorical(rng.choice(list("abc"), n))}, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(rng.poisson(2, (n, 5)).astype(np.float32), obs=obs)


def test_h5ad_round_trip(tmp_path):
    a = _adata()
    path = tmp_path / "sub" / "t.h5ad"  # parent created
    write_adata(path, a)
    r = read_adata(path)
    assert np.array_equal(r.X, a.X)
    pd.testing.assert_frame_equal(r.obs, a.obs)
    assert [p.name for p in path.parent.iterdir()] == ["t.h5ad"]  # no tmp left


def test_zarr_dir_read(tmp_path):
    a = _adata()
    a.write_zarr(tmp_path / "t.zarr")
    r = read_adata(tmp_path / "t.zarr")
    assert np.array_equal(r.X, a.X)
    pd.testing.assert_frame_equal(read_obs(tmp_path / "t.zarr"), a.obs)


def test_read_backed(tmp_path):
    write_adata(tmp_path / "t.h5ad", _adata())
    r = read_adata(tmp_path / "t.h5ad", backed="r")
    assert r.isbacked
    r.file.close()


def test_read_obs_closes_file(tmp_path):
    path = tmp_path / "t.h5ad"
    write_adata(path, _adata())
    obs = read_obs(path)
    assert list(obs["label"].cat.categories) == ["a", "b", "c"]
    write_adata(path, _adata(n=9))  # file not held open: replacing it works
    assert len(read_obs(path)) == 9


def test_write_replaces_existing(tmp_path):
    path = tmp_path / "t.h5ad"
    write_adata(path, _adata(n=20))
    write_adata(path, _adata(n=7))
    assert read_adata(path).n_obs == 7


def test_failed_write_keeps_original(tmp_path, monkeypatch):
    path = tmp_path / "t.h5ad"
    write_adata(path, _adata(n=20))

    def boom(self, filename, *a, **k):
        open(filename, "wb").write(b"partial")
        raise RuntimeError("disk full")

    monkeypatch.setattr(ad.AnnData, "write_h5ad", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        write_adata(path, _adata(n=7))
    monkeypatch.undo()
    assert read_adata(path).n_obs == 20
    assert [p.name for p in tmp_path.iterdir()] == ["t.h5ad"]  # partial tmp removed
