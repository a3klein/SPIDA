"""annotate_adata with --plot, with and without --filter_annotations
(HTH step 5 failed on the unfiltered path, 2026-11-11: plot_annot asked for total_filt regardless)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from spida.I import mmc

LEVELS = {"class": ["A", "B"], "subclass": ["A1", "A2", "B1"]}


def _query(n: int = 60) -> ad.AnnData:
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=sparse.csr_matrix(rng.poisson(3, (n, 10)).astype(float)),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(10)]),
    )
    adata.obsm["base_umap"] = rng.normal(size=(n, 2))
    adata.obsm["spatial"] = rng.uniform(0, 100, (n, 2))
    return adata


def _mmc_csv(path, obs_names):
    # the two columns per level that MMCAnnotator.transfer_labels reads
    rng = np.random.default_rng(1)
    df = pd.DataFrame(index=obs_names)
    for level, labels in LEVELS.items():
        df[f"{level}_name"] = rng.choice(labels, len(obs_names))
        df[f"{level}_bootstrapping_probability"] = rng.uniform(0, 1, len(obs_names))
    df.to_csv(path)
    return path


@pytest.mark.parametrize("filter_annot", [False, True])
def test_annotate_adata_plots_with_and_without_filter(tmp_path, monkeypatch, filter_annot):
    query = _query()
    csv = _mmc_csv(tmp_path / "csv_results.csv", query.obs_names)
    monkeypatch.setattr(mmc.MMCConfig, "validate", lambda self: None)
    monkeypatch.setattr(mmc.MMCAnnotator, "annotate", lambda self, *a, **k: (None, csv))

    pdf = tmp_path / "plots" / "annot.pdf"
    out = mmc.annotate_adata(
        query, "ref", mmc.MMCConfig(mmc_store_path=tmp_path / "store"),
        tmp_path / "query.h5ad", tmp_path / "out",
        filter_annot=filter_annot, plot_path=pdf,
    )
    assert pdf.is_file() and pdf.stat().st_size > 0
    assert {f"mmc_ref_{lvl}" for lvl in LEVELS} <= set(out.obs.columns)
    assert ("mmc_ref_total_filt" in out.obs.columns) == filter_annot
