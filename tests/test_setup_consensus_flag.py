"""Tests for the --consensus/--no-consensus wiring on setup-adata-region.

The flag exists because ALLCools consensus clustering fans out over ``n_jobs=24``
workers whose memory scales with cell count: measured 52 GB at ~133k cells, and
above ~150k it exhausts a 62 GB node. Critically the failure is an OOM that kills
the NODE, so SLURM blames the node and requeues the job -- the ``except`` fallback
inside ``_calc_embeddings`` never fires. That makes this a wiring bug class worth
pinning: a flag that is silently dropped somewhere in
cli -> P.main -> run_setup -> _calc_embeddings looks exactly like it worked.

GPU-free and data-free: the clustering itself is mocked, so what is asserted is
only which branch the flag selects and that the value survives each hand-off.
"""

from __future__ import annotations

import pytest

from spida.io.layout import RegionMeta
from spida.utilities import ad_utils


# ---------------------------------------------------------------------------
# _calc_embeddings: consensus_cluster -> which clustering call happens
# ---------------------------------------------------------------------------

@pytest.fixture
def stub_clustering(monkeypatch):
    """Neutralise everything in _calc_embeddings except the clustering branch."""
    calls = {"consensus": 0, "leiden": 0, "leiden_kwargs": None}

    def fake_consensus(adata, **kw):
        calls["consensus"] += 1
        return adata, None

    def fake_leiden(adata, **kw):
        calls["leiden"] += 1
        calls["leiden_kwargs"] = kw

    monkeypatch.setattr(ad_utils, "_consensus_cluster", fake_consensus)
    monkeypatch.setattr(ad_utils.sc.tl, "leiden", fake_leiden)
    return calls


def _run_calc(monkeypatch, stub, consensus_cluster):
    """Call _calc_embeddings with the embedding/PCA machinery stubbed out."""
    import anndata as ad
    import numpy as np

    adata = ad.AnnData(np.zeros((4, 3), dtype="float32"))
    adata.obsm["base_pca"] = np.zeros((4, 2), dtype="float32")

    # Passing use_rep skips the PCA / significant_pc_test block entirely, and
    # do_umap/do_tsne=False skip the embeddings, so the graph build is the only
    # thing left between the entry point and the clustering branch.
    monkeypatch.setattr(ad_utils.sc.pp, "neighbors", lambda *a, **k: None)

    ad_utils._calc_embeddings(
        adata,
        use_rep="base_pca",
        key_added="base_",
        consensus_cluster=consensus_cluster,
        do_umap=False,
        do_tsne=False,
    )
    return adata


def test_consensus_true_uses_consensus_clustering(monkeypatch, stub_clustering):
    _run_calc(monkeypatch, stub_clustering, consensus_cluster=True)
    assert stub_clustering["consensus"] == 1
    assert stub_clustering["leiden"] == 0


def test_consensus_false_uses_plain_leiden(monkeypatch, stub_clustering):
    _run_calc(monkeypatch, stub_clustering, consensus_cluster=False)
    assert stub_clustering["consensus"] == 0
    assert stub_clustering["leiden"] == 1


def test_plain_leiden_writes_the_same_obs_key(monkeypatch, stub_clustering):
    """Downstream (rank_genes_groups, the seg_clust figures) reads base_leiden.

    If the no-consensus path wrote a different key the plots would silently lose
    their cluster labels rather than fail, so pin the key explicitly.
    """
    _run_calc(monkeypatch, stub_clustering, consensus_cluster=False)
    assert stub_clustering["leiden_kwargs"]["key_added"] == "base_leiden"


# ---------------------------------------------------------------------------
# The hand-offs: cli -> P.main -> run_setup -> _calc_embeddings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value", [True, False])
def test_run_setup_forwards_consensus_to_calc_embeddings(monkeypatch, value):
    import anndata as ad
    import numpy as np
    import scipy.sparse as scp
    from spida.P import setup_adata

    seen = {}

    def fake_calc(adata, **kw):
        seen.update(kw)

    monkeypatch.setattr(setup_adata, "_calc_embeddings", fake_calc)
    monkeypatch.setattr(setup_adata, "normalize_adata", lambda *a, **k: None)
    monkeypatch.setattr(setup_adata.sc.tl, "rank_genes_groups", lambda *a, **k: None)

    adata = ad.AnnData(scp.csr_matrix(np.ones((3, 2), dtype="float32")))
    adata.obs["pass_qc"] = True
    adata.obs["volume"] = 1.0
    adata.obs["CENTER_X"] = 0.0
    adata.obs["CENTER_Y"] = 0.0

    setup_adata.run_setup(adata, RegionMeta("exp", "reg", "seg"), consensus_cluster=value)
    assert seen["consensus_cluster"] is value


@pytest.mark.parametrize("value", [True, False])
def test_main_forwards_consensus_to_run_setup(monkeypatch, value):
    from spida.P import main as p_main
    import spida.P.setup_adata as setup_adata

    seen = {}

    def fake_run_setup(adata, *a, **kw):
        seen.update(kw)
        return adata

    from spida._constants import TABLE_KEY

    from pathlib import Path
    from spida.io.layout import Dataset

    monkeypatch.setattr(setup_adata, "run_setup", fake_run_setup)
    monkeypatch.setattr(p_main, "read_table", lambda *a, **k: object())
    monkeypatch.setattr(p_main, "write_table", lambda *a, **k: None)  # setup writes via write_table
    monkeypatch.setattr(p_main, "require", lambda *a, **k: True)        # phase 3e: store is fake here

    p_main.setup_adata_region(Dataset(exp="exp", reg="reg", zarr=Path("/nonexistent.zarr")), "pfx",
                              consensus_cluster=value)
    assert seen["consensus_cluster"] is value


@pytest.mark.parametrize(
    "argv_flag,expected",
    [("--consensus", True), ("--no-consensus", False), (None, True)],
)
def test_cli_flag_maps_to_consensus_cluster(monkeypatch, tmp_path, argv_flag, expected):
    """Default (no flag) must stay True so existing invocations are unchanged."""
    from click.testing import CliRunner
    from spida.P import cli as p_cli
    from spida.P import main as p_main

    seen = {}
    monkeypatch.setattr(
        p_main, "setup_adata_region", lambda *a, **kw: seen.update(kw)
    )

    args = ["setup-adata-region", "pfx", "--exp", "exp", "--reg", "reg"]
    if argv_flag:
        args.append(argv_flag)

    result = CliRunner().invoke(p_cli.cli, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output
    assert seen["consensus_cluster"] is expected
