"""generate-load-qc-figs must survive a region processed without transcript-qc,
cluster-hexes, and/or tissue-mask, plotting whatever inputs are present."""

from types import SimpleNamespace
from unittest import mock

import click
import pytest

SHAPES_KEY = "transcript_qc_shapes"
HEX_KEY = "adata_hex_s30.0_o0.0_clustered"
MASK_KEY = "qc_mask"


class _FakeSData(dict):
    def __init__(self, shapes=(), tables=()):
        super().__init__()
        self.shapes = {k: SimpleNamespace(copy=lambda: "grid") for k in shapes}
        self.tables = {
            k: SimpleNamespace(copy=lambda: SimpleNamespace(obs="obs")) for k in tables
        }
        self.update(self.shapes)
        self.update(self.tables)


@pytest.mark.parametrize(
    "shapes, tables, want_tz, want_hex, want_mask",
    [
        ((SHAPES_KEY, MASK_KEY), (HEX_KEY,), True, True, True),
        ((SHAPES_KEY,), (HEX_KEY,), True, True, False),
        ((SHAPES_KEY,), (), True, False, False),
        ((), (HEX_KEY,), False, True, False),
        ((), (), False, False, False),
        ((MASK_KEY,), (), False, False, True),
    ],
    ids=["everything", "no_mask", "no_cluster_hexes", "no_transcript_qc",
         "neither", "mask_only"],
)
def test_optional_inputs_are_skipped_not_fatal(shapes, tables, want_tz, want_hex, want_mask):
    from spida.site import generate_load_qc_figs

    sdata = _FakeSData(shapes, tables)
    with mock.patch("spida.site._load_metadata",
                    return_value=(sdata, "/img", "/data",
                                  ("img_key", "shapes_key", "table_key", "points_key"), None)), \
         mock.patch("spida.P.transcript_qc._obs_to_grid_geodf", return_value="geodf"), \
         mock.patch("spida.pl.site_figs.plot_load_images") as load_images, \
         mock.patch("spida.pl.site_figs.plot_tz_qc") as tz_qc, \
         mock.patch("spida.pl.site_figs.plot_tz_hex_qc") as hex_qc, \
         mock.patch("spida.pl.site_figs.plot_tissue_mask") as tissue_mask, \
         mock.patch("spida.pl.site_figs.plot_decon_images") as decon:
        # pass_context injects ctx itself, so the callback needs a live context.
        with click.Context(generate_load_qc_figs):
            generate_load_qc_figs.callback(
                exp="EXP", reg="REG", zarr_path="/nonexistent.zarr", raw_path=None, seg_path=None,
                figures_path=None, lab=None, brain_region="WB",
                zarr_store=None, site_dir="/site", plot_decon=False, naming_map=None,
            )

    assert load_images.called
    assert tz_qc.called is want_tz
    assert hex_qc.called is want_hex
    assert tissue_mask.called is want_mask
    assert not decon.called
