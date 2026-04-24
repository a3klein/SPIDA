from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import re

from spida.utilities.sd_utils import _gen_keys
from spida._constants import IMAGE_KEY, SHAPES_KEY, POINTS_KEY, TABLE_KEY


def load_naming_map(path: str | Path | None) -> pd.DataFrame | None:
    """Load the naming_map.csv, indexed by Name column. Returns None if path is None."""
    if path is None:
        return None
    return pd.read_csv(path, index_col="Name")


def lookup_naming_entry(
    naming_map: pd.DataFrame | None,
    exp_name: str,
    reg_name: str,
    lab: str | None = None,
) -> dict | None:
    """
    Return {brain_region, exp_alias, donor, lab} for the given experiment/region/lab,
    or None if the naming map is absent or has no matching row.
    """
    if naming_map is None:
        return None
    key = f"{exp_name}_{reg_name}_{lab}" if lab is not None else f"{exp_name}_{reg_name}"
    if key not in naming_map.index:
        return None
    row = naming_map.loc[key]
    return {
        "brain_region": row["Brain Region"],
        "exp_alias": row["EXP name"],
        "donor": row["Donor name"],
        "lab": row["Lab"],
    }


def _load_metadata(
    exp_name: str,
    reg_name: str,
    prefix_name: str | None = None,
    lab: str | None = None,
    brain_region: str = "WB",
    zarr_store: str | None = None,
    site_dir: str | Path | None = None,
    naming_map: str | Path | None = None,
):
    import spatialdata as sd

    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    if site_dir is None:
        site_dir = os.getenv("SPIDA_SITE_DIR", None)
    if isinstance(site_dir, str):
        site_dir = Path(site_dir)

    naming_map_df = load_naming_map(naming_map)
    entry = lookup_naming_entry(naming_map_df, exp_name, reg_name, lab)
    br = (entry["brain_region"] if entry else brain_region) or "WB"

    dir_name = f"{exp_name}_{reg_name}" if lab is None else f"{exp_name}_{reg_name}_{lab}"
    img_dir = site_dir / "images" / br / dir_name
    img_dir.mkdir(parents=True, exist_ok=True)
    data_dir = site_dir / "data" / br / dir_name
    data_dir.mkdir(parents=True, exist_ok=True)
    sdata_path = Path(zarr_store) / exp_name / reg_name
    sdata = sd.read_zarr(sdata_path)

    keys = _gen_keys("default", exp_name, reg_name)
    image_key = keys[IMAGE_KEY]
    if prefix_name is not None:
        keys = _gen_keys(prefix_name, exp_name, reg_name)
        shapes_key = keys[SHAPES_KEY]
        table_key = keys[TABLE_KEY]
        points_key = keys[POINTS_KEY]
    else:
        shapes_key = None
        table_key = None
        points_key = None

    return sdata, img_dir, data_dir, (image_key, shapes_key, table_key, points_key), naming_map_df


def generate_soma_gene_proportions(points, segmentation_key: str, out_dir: str | Path):
    if "cell_id" in points.columns:
        points["in_cell"] = (points["cell_id"] > 0).map({True: "soma", False: "outside"}).astype("category")
        use_col = "cell_id"
    else:
        points["in_cell"] = (~points["assignment"].isna()).map({True: "soma", False: "outside"}).astype("category")
        use_col = "x"
    gass = points.groupby(["in_cell", "gene"], observed=False)[use_col].count().reset_index()
    gass = gass.pivot(index="gene", columns="in_cell", values=use_col).fillna(0)
    gass_norm = gass.apply(lambda x: x / x.sum(), axis=1).sort_values(by="soma")
    gass_norm["total counts"] = gass.sum(axis=1)
    gass_norm["log total counts"] = np.log10(gass_norm["total counts"] + 1)
    out_path = Path(out_dir) / f"{segmentation_key}_soma_gene_proportions.csv"
    gass_norm.to_csv(out_path)
    return gass_norm


def append_brain_region_qc_metrics(
    adata,
    brain_region: str,
    exp_name: str,
    reg_name: str,
    segmentation: str,
    site_dir: str | Path | None,
    lab: str | None = None,
    neuron_type_col: str | None = None,
    naming_map: pd.DataFrame | None = None,
) -> None:
    if site_dir is None:
        site_dir = os.getenv("SPIDA_SITE_DIR", None)
    if isinstance(site_dir, str):
        site_dir = Path(site_dir)
    if site_dir is None:
        return

    entry = lookup_naming_entry(naming_map, exp_name, reg_name, lab)
    if entry is not None:
        brain_region_alias = entry["exp_alias"]
        donor = entry["donor"]
        replicate = entry["lab"]
        br = entry["brain_region"] or "WB"
    else:
        alias_match = re.search(r"4x1-([^/]+)-(?:E|Q)", exp_name)
        brain_region_alias = alias_match.group(1) if alias_match else exp_name
        donor_match = re.search(r"region_([^/_]+)", reg_name)
        donor = donor_match.group(1) if donor_match else reg_name
        replicate = lab or "unknown"
        br = brain_region if brain_region is not None else "WB"

    data_root = site_dir / "data" / br
    data_root.mkdir(parents=True, exist_ok=True)
    output_file = data_root / "qc_metrics.csv"

    df_obs = adata.obs.copy()
    if neuron_type_col and neuron_type_col in df_obs.columns:
        df_obs["neuron_type"] = df_obs[neuron_type_col].astype(str)
    else:
        df_obs["neuron_type"] = "all"

    required_cols = ["nCount_RNA", "nFeature_RNA", "nCount_RNA_per_Volume", "volume"]
    missing = [col for col in required_cols if col not in df_obs.columns]
    if missing:
        return

    df_obs = df_obs[required_cols + ["neuron_type"]].copy()
    df_obs["brain_region"] = brain_region_alias
    df_obs["brain_region_alias"] = brain_region_alias
    df_obs["donor"] = donor
    df_obs["replicate"] = replicate
    df_obs["experiment"] = exp_name
    df_obs["region"] = reg_name
    df_obs["segmentation"] = segmentation

    group_cols = [
        "brain_region",
        "donor",
        "replicate",
        "experiment",
        "region",
        "segmentation",
        "neuron_type",
    ]

    df_count = df_obs.groupby(group_cols)["nCount_RNA"].median().reset_index(name="nCount_RNA_median")
    df_ft = df_obs.groupby(group_cols)["nFeature_RNA"].median().reset_index(name="nFeature_RNA_median")
    df_v = df_obs.groupby(group_cols)["volume"].median().reset_index(name="volume_median")
    df_tot = df_obs.groupby(group_cols)["nCount_RNA"].count().reset_index(name="total_cells")

    df_out = df_count.merge(df_ft, on=group_cols, how="inner")
    df_out = df_out.merge(df_v, on=group_cols, how="inner")
    df_out = df_out.merge(df_tot, on=group_cols, how="inner")

    if output_file.exists():
        existing = pd.read_csv(output_file)
        existing = existing[
            ~(
                (existing["experiment"] == exp_name)
                & (existing["region"] == reg_name)
                & (existing["segmentation"] == segmentation)
            )
        ]
        df_out = pd.concat([existing, df_out], ignore_index=True)

    df_out.to_csv(output_file, index=False)
