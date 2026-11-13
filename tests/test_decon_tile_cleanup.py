"""Decon tile dir cleanup: removed after a successful run only when empty."""
from spida.S.image_prep.orchestrator import _remove_if_empty


def test_empty_tile_dir_removed(tmp_path):
    d = tmp_path / "decon-tiles-DAPI"
    d.mkdir()
    _remove_if_empty(d)
    assert not d.exists()


def test_tile_dir_with_files_kept(tmp_path):
    # e.g. tiles of another channel, or anything a --continue_stalled run needs
    d = tmp_path / "decon-tiles-DAPI"
    d.mkdir()
    (d / "tile_PolyT_0_0.tif").write_bytes(b"x")
    _remove_if_empty(d)
    assert (d / "tile_PolyT_0_0.tif").exists()


def test_missing_dir_is_fine(tmp_path):
    _remove_if_empty(tmp_path / "never_made")
