import os
import warnings

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")
sd = pytest.importorskip("spatialdata")

import dask.dataframe as dd
import tifffile
from spatialdata.models import PointsModel

from spida.io.layout import Dataset, element_keys
from spida._constants import POINTS_KEY

PTS = element_keys("default", None)[POINTS_KEY]
KW = dict(a_min_um2=100.0, gb=5.0)  # small synthetic tissue


@pytest.fixture
def ds(tmp_path):
    rng = np.random.default_rng(0)
    n = 20000
    pts = pd.DataFrame({"x": rng.uniform(0, 300, n), "y": rng.uniform(0, 300, n), "gene": rng.choice(list("AB"), n)})
    s = sd.SpatialData(points={PTS: PointsModel.parse(dd.from_pandas(pts, npartitions=1), feature_key="gene")})
    s.attrs["spida"] = {"namespace": "", "exp": None, "reg": None}
    d = Dataset(zarr=tmp_path / "s.zarr", raw=tmp_path / "raw")
    s.write(d.zarr)
    img = d.raw_images()
    img.mkdir(parents=True)
    tifffile.imwrite(img / "mosaic_DAPI_z0.tif", np.zeros((300, 300), np.uint16))
    (img / "micron_to_mosaic_pixel_transform.csv").write_text("1 0 0\n0 1 0\n0 0 1\n")
    return d


def _run(d, **k):
    from spida.P.main import tissue_mask_region
    tissue_mask_region(d, **KW, **k)


def _state(d):
    r = sd.read_zarr(d.zarr)
    tif = d.raw_images() / "tissue_mask.tif"
    return len(r["qc_mask"]), int(r["tissue_mask"].values.sum()), tif.stat().st_mtime


def test_first_run_writes_all_three(ds):
    _run(ds)
    r = sd.read_zarr(ds.zarr)
    assert "qc_mask" in r.shapes and "tissue_mask" in r.labels
    assert (ds.raw_images() / "tissue_mask.tif").is_file()


def test_rerun_refuses_and_lists_all(ds):
    _run(ds)
    before = _state(ds)
    with pytest.raises(FileExistsError) as e:
        _run(ds, rho_min=1.0)
    msg = str(e.value)
    assert "shapes/qc_mask" in msg and "labels/tissue_mask" in msg and "tissue_mask.tif" in msg
    assert _state(ds) == before


def test_rerun_overwrite_replaces_all_three(ds):
    _run(ds)
    tif = ds.raw_images() / "tissue_mask.tif"
    os.utime(tif, ns=(tif.stat().st_mtime_ns - 10**10,) * 2)  # coarse clock: make "newer" unambiguous
    t0 = tif.stat().st_mtime_ns
    _run(ds, erode_sigma=3.0, overwrite=True)  # more erosion -> a smaller mask
    r = sd.read_zarr(ds.zarr)
    assert not any("__tmp" in p for p in r.elements_paths_on_disk())
    assert (ds.raw_images() / "tissue_mask.tif").stat().st_mtime_ns > t0
    eroded = int(r["tissue_mask"].values.sum())  # read now: the element is lazy
    _run(ds, overwrite=True)
    assert int(sd.read_zarr(ds.zarr)["tissue_mask"].values.sum()) > eroded


def test_no_raster_ignores_existing_tif(ds):
    (ds.raw_images() / "tissue_mask.tif").write_bytes(b"x")
    _run(ds, write_raster=False)
    assert "qc_mask" in sd.read_zarr(ds.zarr).shapes
