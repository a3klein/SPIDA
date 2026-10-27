"""File names SPIDA reads and writes, kept in io/formats and io/schema (phase 5c-2)."""
from pathlib import Path

import pytest

from spida.io.formats.merscope import PLANE_RE, TRANSCRIPTS_CSV, TRANSCRIPTS_PARQUET, raw_transcripts
from spida.io.schema.decon import plane_name


def test_raw_transcripts_prefers_parquet(tmp_path):
    (tmp_path / TRANSCRIPTS_CSV).write_text("")
    assert raw_transcripts(tmp_path) == tmp_path / TRANSCRIPTS_CSV
    (tmp_path / TRANSCRIPTS_PARQUET).write_text("")
    assert raw_transcripts(tmp_path) == tmp_path / TRANSCRIPTS_PARQUET


@pytest.mark.parametrize("name,stain,z", [("mosaic_DAPI_z3.tif", "DAPI", "3"), ("mosaic_Cellbound1-2_z12.tif", "Cellbound1-2", "12")])
def test_plane_re_matches_raw_planes(name, stain, z):
    m = PLANE_RE.search(name)
    assert m and (m["stain"], m["z"]) == (stain, z)


@pytest.mark.parametrize("name", ["mosaic_DAPI_z3.decon.tif", "mosaic_DAPI_z_stack_z3.decon.tif",
                                  "mosaic_DAPI_z_stack.tif", "mosaic_DAPI_z3_masks.tif"])
def test_plane_re_skips_derived_images(name):
    assert PLANE_RE.search(name) is None


# The writer (S/image_prep/orchestrator.transform_image) names its stitched outputs from the input planes:
#   3D: stack = f"{first_plane.stem[:-1]}_stack.tif"; out = f"{Path(stack).stem}_z{z:01d}.decon.tif"
#   2D: out = f"{plane.stem}.decon.tif"
# These mirror those two expressions; if the writer changes, change them here and in io/schema/decon.py.
def _writer_3d(planes, z):
    stack = Path(planes[0]).parent / f"{Path(planes[0]).stem[:-1]}_stack.tif"
    return f"{stack.stem}_z{z:01d}.decon.tif"


def _writer_2d(plane):
    return f"{Path(plane).stem}.decon.tif"


@pytest.mark.parametrize("stain", ["DAPI", "PolyT"])
def test_decon_names_reader_matches_writer(stain):
    planes = sorted(f"mosaic_{stain}_z{i}.tif" for i in range(7))  # THM2: z0..z6
    for z in range(7):
        assert plane_name(stain, z, z_stack=True) == _writer_3d(planes, z)
    assert plane_name(stain, 3) == _writer_2d(f"mosaic_{stain}_z3.tif")  # PU 2D: z3 only


def test_decon_reader_names_real_thm2_files():
    """Names seen on disk (THM2 UWA7648 gold raw images, 2026-10-25)."""
    on_disk = {f"mosaic_{s}_z_stack_z{z}.decon.tif" for s in ("DAPI", "PolyT") for z in range(7)}
    assert {plane_name(s, z, z_stack=True) for s in ("DAPI", "PolyT") for z in range(7)} == on_disk
