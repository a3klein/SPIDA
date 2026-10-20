import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from spida.config import SpidaConfig
from spida.io.layout import dataset_from


@pytest.fixture
def cfg_file(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({k: str(tmp_path / v) for k, v in {
        "ZARR_STORAGE_PATH": "z", "ANNDATA_STORE_PATH": "ad", "ANNOTATIONS_STORE_PATH": "an",
        "IMAGE_STORE_PATH": "img", "MMC_DIR": "mmc"}.items()}))
    (tmp_path / "mmc").mkdir()
    return p


def _table(n=4):
    return ad.AnnData(np.ones((n, 2), dtype=np.float32), obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))


def _zarr_with_table(zarr, key):
    (zarr / "tables").mkdir(parents=True)
    _table().write_zarr(zarr / "tables" / key)


def test_dataset_has_table_locations(tmp_path):
    cfg = SpidaConfig(anndata_store=tmp_path / "ad", annotation_store=tmp_path / "an")
    d = dataset_from(cfg, "E", "region_1")
    assert d.table_h5ad("p_table_filt") == tmp_path / "ad" / "E" / "region_1" / "p_table_filt.h5ad"
    assert d.annotation_dir("p_table_filt", "mmc", "ref") == tmp_path / "an" / "E" / "region_1" / "p_table_filt" / "mmc" / "ref"
    d2 = dataset_from(cfg, anndata=tmp_path / "x", annotations=tmp_path / "y")
    assert d2.table_h5ad("k") == tmp_path / "x" / "k.h5ad" and d2.annotation_dir("k", "m", "r") == tmp_path / "y" / "k" / "m" / "r"


def test_load_table_copies_then_reuses(tmp_path):
    from spida.I._io import load_table

    d = dataset_from(SpidaConfig(), zarr=tmp_path / "s.zarr", anndata=tmp_path / "ad")
    _zarr_with_table(tmp_path / "s.zarr", "p_table_filt")
    a, path = load_table(d, "p", "_filt")
    assert path == tmp_path / "ad" / "p_table_filt.h5ad" and path.exists()
    a.obs["mmc_r_Class"] = "x"
    a.write_h5ad(path)
    assert "mmc_r_Class" in load_table(d, "p", "_filt")[0].obs          # copy reused
    assert "mmc_r_Class" not in load_table(d, "p", "_filt", refresh=True)[0].obs


def test_refresh_table_logs_dropped_annotations(tmp_path, caplog):
    from spida.I._io import load_table

    d = dataset_from(SpidaConfig(), zarr=tmp_path / "s.zarr", anndata=tmp_path / "ad")
    _zarr_with_table(tmp_path / "s.zarr", "p_table_filt")
    a, path = load_table(d, "p", "_filt")
    with caplog.at_level("WARNING"):
        load_table(d, "p", "_filt", refresh=True)  # nothing to drop yet
    assert "drops" not in caplog.text
    a.obs["mmc_r_Class"] = "x"
    a.obs["allcools_AIT_Subclass"] = "y"
    a.obs["base_leiden"] = "0"
    a.write_h5ad(path)
    with caplog.at_level("WARNING"):
        load_table(d, "p", "_filt", refresh=True)
    assert "drops 2 annotation columns: mmc_r_Class, allcools_AIT_Subclass" in caplog.text


def test_load_table_adata_file_needs_no_store(tmp_path):
    from spida.I._io import load_table

    f = tmp_path / "collab.h5ad"
    _table().write_h5ad(f)
    a, path = load_table(dataset_from(SpidaConfig(), annotations=tmp_path / "an"), "p", adata_file=f)
    assert path == f and a.n_obs == 4


def test_tag_columns():
    from spida.I._io import tag_columns

    a = _table()
    a.obs["allcools_Subclass"] = "x"
    a.obs["c2c_allcools_label_Subclass"] = "y"
    a.obs["all_annot"] = "z"
    a.obs["allcools_refA_Subclass"] = "stale"
    names = tag_columns(a, ["allcools_Subclass", "c2c_allcools_label_Subclass", "all_annot"], "allcools", "refA")
    assert names == {"allcools_Subclass": "allcools_refA_Subclass",
                     "c2c_allcools_label_Subclass": "allcools_refA_c2c_label_Subclass",
                     "all_annot": "allcools_refA_all_annot"}
    assert list(a.obs.columns) == list(names.values()) and (a.obs["allcools_refA_Subclass"] == "x").all()


def test_annotate_region_writes_schema_paths(tmp_path, monkeypatch):
    from spida.I import mmc

    d = dataset_from(SpidaConfig(zarr_store=tmp_path / "z", anndata_store=tmp_path / "ad",
                                 annotation_store=tmp_path / "an"), "E", "region_1")
    _zarr_with_table(d.zarr, "p_table_filt")
    seen = {}

    def fake_annotate(self, query_path, identifier, output_dir, **kw):
        seen.update(query=query_path, out=output_dir)
        csv = tmp_path / "res.csv"
        idx = ad.read_h5ad(query_path).obs_names
        pd.DataFrame({"Class_name": "A", "Class_bootstrapping_probability": 0.9}, index=idx).to_csv(csv)
        return None, csv

    monkeypatch.setattr(mmc.MMCAnnotator, "annotate", fake_annotate)
    config = mmc.MMCConfig(mmc_store_path=tmp_path)
    for ref in ("refA", "refB"):
        mmc.annotate_region(d, "p", "_filt", ref, config)

    copy = d.table_h5ad("p_table_filt")
    assert seen["query"] == copy
    assert seen["out"] == tmp_path / "an" / "E" / "region_1" / "p_table_filt" / "mmc" / "refB"
    obs = ad.read_h5ad(copy).obs
    assert {"mmc_refA_Class", "mmc_refB_Class", "mmc_refA_Class_transfer_score"} <= set(obs.columns)


def test_mmc_region_cli(monkeypatch, cfg_file, tmp_path):
    from spida.I import cli as i_cli

    calls = []
    monkeypatch.setattr("spida.I.mmc.annotate_region", lambda *a, **k: calls.append((a, k)))
    res = CliRunner().invoke(i_cli.cli, ["--config", str(cfg_file), "mmc-annotation-region", "p", "refA",
                                         "--exp", "E", "--reg", "region_1", "--l1_thr=0.7"])
    assert res.exit_code == 0, res.output
    (d, prefix, suffix, ident, config), k = calls[0]
    assert (prefix, suffix, ident) == ("p", "_filt", "refA")
    assert d.anndata == tmp_path / "ad" / "E" / "region_1" and d.annotations == tmp_path / "an" / "E" / "region_1"
    assert config.mmc_store_path == tmp_path / "mmc" and k["l1_thr"] == 0.7


def test_mmc_region_cli_adata_file_unnamed(monkeypatch, cfg_file, tmp_path):
    from spida.I import cli as i_cli

    f = tmp_path / "collab.h5ad"
    _table().write_h5ad(f)
    calls = []
    monkeypatch.setattr("spida.I.mmc.annotate_region", lambda *a, **k: calls.append((a, k)))
    res = CliRunner().invoke(i_cli.cli, ["--config", str(cfg_file), "mmc-annotation-region", "collab", "refA",
                                         "--adata_file", str(f), "--annotations_path", str(tmp_path / "out")])
    assert res.exit_code == 0, res.output
    (d, *_), k = calls[0]
    assert not d.named and d.annotations == tmp_path / "out" and k["adata_file"] == str(f)


def test_i_store_options_are_config_fields():
    import click
    from spida.config import CONFIG_KEY_TO_FIELD
    from spida.I.cli import cli

    mapped = set(CONFIG_KEY_TO_FIELD.values())
    assert not {p.name for c in cli.commands.values() for p in c.params
                if isinstance(p, click.Option) and "store" in p.name} - mapped


def test_integrate_adata_tags_query_columns(tmp_path, monkeypatch):
    pytest.importorskip("ALLCools")
    from spida.I import allcools

    ref = _table(6)
    ref.write_h5ad(tmp_path / "ref.h5ad")

    def fake_run(ref_adata, qry_adata, output_dir=None, **kw):
        qry_adata.obs["allcools_Subclass"] = "A"
        qry_adata.obs["c2c_allcools_label_Subclass"] = "A"
        qry_adata.obsm["allcools_Subclass_score_matrix"] = np.zeros((qry_adata.n_obs, 1))
        return qry_adata, None

    monkeypatch.setattr(allcools, "run_allcools_seurat", fake_run)
    qry, _, _ = allcools.integrate_adata(_table(), tmp_path / "ref.h5ad", "refA", tmp_path / "q.h5ad", None,
                                         rna_cell_type_column="Subclass")
    obs = ad.read_h5ad(tmp_path / "q.h5ad").obs
    assert {"allcools_refA_Subclass", "allcools_refA_c2c_label_Subclass"} == set(obs.columns)
    assert "allcools_refA_Subclass_score_matrix" in qry.obsm
