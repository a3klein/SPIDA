"""Archived 2026-09-27 (IO redesign phase 1, I module): legacy MMC wrappers from I/mmc.py (lines 894-1136), I/mmc_cli.py, I/_utils.py. Not importable as-is."""

# ==== I/mmc.py legacy block ====

# ============================================================================
# Legacy Functions (Backward Compatibility)
# ============================================================================

def setup_mmc(
    ref_path: str,
    brain_region: str,
    codebook: str,
    hierarchy_list: list,
    codebook_path: str = None,
    mmc_store_path: str = None,
    ref_norm: str = "log2CPM",
    **kwargs,
):
    """
    Legacy setup function for MapMyCells integration.
    
    DEPRECATED: Use setup_annotation_pipeline() instead.
    
    Parameters:
    ref_path (str): Path to the reference data.
    brain_region (str): Brain region for the annotation.
    codebook (str): Codebook for the annotation.
    hierarchy_list (list): List of hierarchy levels for the annotation.
    codebook_path (str, optional): Path to the codebook file. Defaults to None.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    ref_norm (str, optional): Normalization method for the reference AnnData.X. Defaults to "log2CPM".
    **kwargs: Additional keyword arguments.
    """
    logger.warning(
        "setup_mmc() is deprecated. Use setup_annotation_pipeline() instead."
    )
    
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
        if not mmc_store_path:
            raise ValueError(
                "Please provide a mmc_store_path or set the MMC_DIR environment variable."
            )

    if not codebook_path:
        gene_panel_path = os.getenv("GENE_PANEL_PATH")
        if not gene_panel_path:
            raise ValueError(
                "Please provide a codebook path or set the GENE_PANEL_PATH environment variable."
            )
        codebook_path = glob.glob(f"{gene_panel_path}/*{codebook}*.csv")[0]

    logger.info(f"ref_path: {ref_path}")
    logger.info(f"hierarchy_list: {hierarchy_list}")
    logger.info(f"brain_region: {brain_region}")
    logger.info(f"codebook: {codebook}")
    logger.info(f"codebook_path: {codebook_path}")
    logger.info(f"mmc_store_path: {mmc_store_path}")
    logger.info(f"ref_norm: {ref_norm}")

    config = MMCConfig(
        mmc_store_path=Path(mmc_store_path),
        anndata_store_path=Path(os.getenv("ANNDATA_STORE_PATH", "")),
        annotations_store_path=Path(os.getenv("ANNOTATIONS_STORE_PATH", "")),
        zarr_storage_path=Path(os.getenv("ZARR_STORAGE_PATH", "")),
        gene_panel_path=Path(os.getenv("GENE_PANEL_PATH", "")),
        ref_norm=ref_norm,
        **kwargs,
    )
    
    preprocessor = MMCPreprocessor(config)
    preprocessor.setup_annotation_pipeline(
        ref_path, brain_region, codebook, hierarchy_list, codebook_path, **kwargs
    )

    logger.info("DONE")
    return 0


def run_mmc(
    query_adata: ad.AnnData,
    brain_region: str,
    codebook: str,
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Legacy annotation function for MapMyCells.
    
    DEPRECATED: Use annotate_region() or annotate_experiment() instead.
    
    Parameters:
    query_adata (ad.AnnData): Query AnnData object to annotate.
    brain_region (str): Brain region for the annotation.
    codebook (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    logger.warning(
        "run_mmc() is deprecated. Use annotate_region() or annotate_experiment() instead."
    )

    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
        if not mmc_store_path:
            raise ValueError(
                "Please provide a mmc_store_path or set the MMC_DIR environment variable."
            )

    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")
        if not anndata_store_path:
            raise ValueError(
                "Please provide an anndata_store_path or set the ANNDATA_STORE_PATH environment variable."
            )

    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")
        if not annotations_store_path:
            raise ValueError(
                "Please provide an annotations_store_path or set the ANNOTATIONS_STORE_PATH environment variable."
            )

    exp = query_adata.uns.get("experiment")
    seg_name = query_adata.uns.get("segmentation")
    donor = query_adata.uns.get("donor")
    query_path = f"{anndata_store_path}/{exp}/{seg_name}/adata_{donor}.h5ad"
    output_path = f"{annotations_store_path}/{exp}/{seg_name}/mmc/{donor}"
    Path(output_path).mkdir(parents=True, exist_ok=True)
    logger.info(f"Writing query adata to temporary path: {query_path}")
    Path(query_path).parent.mkdir(parents=True, exist_ok=True)
    query_adata.write_h5ad(query_path)

    logger.info("Running MMC: ")
    logger.info(f"mmc_store_path: {mmc_store_path}")
    logger.info(f"brain_region: {brain_region}")
    logger.info(f"codebook: {codebook}")
    logger.info(f"query path: {query_path}")
    logger.info(f"output_path: {output_path}")

    config = MMCConfig(
        mmc_store_path=Path(mmc_store_path),
        anndata_store_path=Path(anndata_store_path),
        annotations_store_path=Path(annotations_store_path),
        zarr_storage_path=Path(os.getenv("ZARR_STORAGE_PATH", "")),
        gene_panel_path=Path(os.getenv("GENE_PANEL_PATH", "")),
        **kwargs,
    )
    
    annotator = MMCAnnotator(config)
    extended_path, csv_path = annotator.annotate(
        query_path, brain_region, codebook, output_path, **kwargs
    )

    logger.info("Transferring labels to the adata object")
    pre_cols = query_adata.obs.columns
    query_adata = annotator.transfer_labels(query_adata, csv_path)
    post_cols = query_adata.obs.columns
    logger.info(f"Added columns to adata.obs: {set(post_cols) - set(pre_cols)}")

    query_adata.write_h5ad(query_path)

    return 0


def mmc_annotation_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str,
    brain_region: str,
    codebook: str,
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Legacy region-level annotation function.
    
    DEPRECATED: Use annotate_region() instead.
    """
    logger.warning(
        "mmc_annotation_region() is deprecated. Use annotate_region() instead."
    )

    # Getting the sdata object
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    adata = ad.read_zarr(f"{zarr_path}/tables/{prefix_name}_table{suffix}")

    return run_mmc(
        adata,
        brain_region,
        codebook,
        mmc_store_path,
        anndata_store_path,
        annotations_store_path,
        **kwargs,
    )


def mmc_annotation_experiment(
    exp_name: str,
    prefix_name: str,
    suffix: str,
    brain_region: str,
    codebook: str,
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Legacy experiment-level annotation function.
    
    DEPRECATED: Use annotate_experiment() instead.
    """
    logger.warning(
        "mmc_annotation_experiment() is deprecated. Use annotate_experiment() instead."
    )

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    for reg in regions:
        reg_name = reg.split("/")[-1]
        mmc_annotation_region(
            exp_name,
            reg_name,
            prefix_name,
            suffix=suffix,
            brain_region=brain_region,
            codebook=codebook,
            mmc_store_path=mmc_store_path,
            anndata_store_path=anndata_store_path,
            annotations_store_path=annotations_store_path,
            **kwargs,
        )

    return 0

# ==== I/mmc_cli.py ====
"""
CLI wrapper functions for MapMyCells integration.

These functions adapt the new MMC architecture to work with Click CLI commands.
"""

import logging
from pathlib import Path
from typing import List, Optional

import anndata as ad

from spida.I.mmc import (
    MMCConfig,
    MMCPreprocessor,
    MMCAnnotator,
    setup_annotation_pipeline,
    annotate_region,
    annotate_experiment,
)

logger = logging.getLogger(__name__)


def mmc_setup_cli(
    ref_path: str,
    identifier: str,
    hierarchy_list: List[str],
    gene_names_path: Optional[str] = None,
    gene_name_mapping_path: Optional[str] = None,
    mmc_store_path: Optional[str] = None,
    ref_norm: str = "log2CPM",
    n_cpu: int = 8,
    n_valid: int = 10,
    n_per_utility: int = 10,
    **kwargs,
) -> dict:
    """
    CLI wrapper for MMCPreprocessor setup.
    
    Parameters match the CLI command interface while delegating to the new architecture.
    
    Returns:
    --------
    dict with paths to precomputed files
    """
    logger.info("Starting MMC setup via CLI")
    
    # Build config
    import os
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
    
    config = MMCConfig(
        mmc_store_path=Path(mmc_store_path),
        anndata_store_path=Path(os.getenv("ANNDATA_STORE_PATH", "")),
        annotations_store_path=Path(os.getenv("ANNOTATIONS_STORE_PATH", "")),
        zarr_storage_path=Path(os.getenv("ZARR_STORAGE_PATH", "")),
        gene_panel_path=Path(os.getenv("GENE_PANEL_PATH", "")),
        ref_norm=ref_norm,
        n_cpu=n_cpu,
    )
    
    result = setup_annotation_pipeline(
        ref_path=ref_path,
        identifier=identifier,
        hierarchy_list=hierarchy_list,
        gene_names_path=gene_names_path,
        gene_name_mapping_path=gene_name_mapping_path,
        config=config,
        n_valid=n_valid,
        n_per_utility=n_per_utility,
        **kwargs,
    )
    
    logger.info(f"MMC setup complete. Results: {result}")
    return result


def mmc_annotation_region_cli(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str,
    identifier: str,
    mmc_store_path: Optional[str] = None,
    anndata_store_path: Optional[str] = None,
    annotations_store_path: Optional[str] = None,
    zarr_store_path: Optional[str] = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    filter_annot: bool = False,
    plot: bool = False,
    palette_path: Optional[str] = None,
    image_store_path: Optional[str] = None,
    **kwargs,
) -> int:
    """
    CLI wrapper for region-level annotation.
    
    Parameters match the CLI command interface while delegating to the new architecture.
    
    Returns:
    --------
    0 on success
    """
    logger.info(f"Starting MMC annotation for region: {exp_name}/{reg_name}")
    
    import os
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")
    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")
    if not zarr_store_path:
        zarr_store_path = os.getenv("ZARR_STORAGE_PATH")
    
    config = MMCConfig(
        mmc_store_path=Path(mmc_store_path),
        anndata_store_path=Path(anndata_store_path),
        annotations_store_path=Path(annotations_store_path),
        zarr_storage_path=Path(zarr_store_path),
        gene_panel_path=Path(os.getenv("GENE_PANEL_PATH", "")),
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
    )
    
    result = annotate_region(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        suffix=suffix,
        identifier=identifier,
        config=config,
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
        filter_annot=filter_annot,
        plot=plot,
        palette_path=palette_path,
        plot_path=image_store_path,
        **kwargs,
    )
    
    logger.info(f"MMC annotation for {exp_name}/{reg_name} complete")
    return result


def mmc_annotation_experiment_cli(
    exp_name: str,
    prefix_name: str,
    suffix: str,
    identifier: str,
    mmc_store_path: Optional[str] = None,
    anndata_store_path: Optional[str] = None,
    annotations_store_path: Optional[str] = None,
    zarr_store_path: Optional[str] = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    filter_annot: bool = False,
    plot: bool = False,
    palette_path: Optional[str] = None,
    image_store_path: Optional[str] = None,
    **kwargs,
) -> int:
    """
    CLI wrapper for experiment-level annotation.
    
    Parameters match the CLI command interface while delegating to the new architecture.
    
    Returns:
    --------
    0 on success
    """
    logger.info(f"Starting MMC annotation for experiment: {exp_name}")
    
    import os
    if not mmc_store_path:
        mmc_store_path = os.getenv("MMC_DIR")
    if not anndata_store_path:
        anndata_store_path = os.getenv("ANNDATA_STORE_PATH")
    if not annotations_store_path:
        annotations_store_path = os.getenv("ANNOTATIONS_STORE_PATH")
    if not zarr_store_path:
        zarr_store_path = os.getenv("ZARR_STORAGE_PATH")
    
    config = MMCConfig(
        mmc_store_path=Path(mmc_store_path),
        anndata_store_path=Path(anndata_store_path),
        annotations_store_path=Path(annotations_store_path),
        zarr_storage_path=Path(zarr_store_path),
        gene_panel_path=Path(os.getenv("GENE_PANEL_PATH", "")),
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
    )
    
    result = annotate_experiment(
        exp_name=exp_name,
        prefix_name=prefix_name,
        suffix=suffix,
        identifier=identifier,
        config=config,
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
        filter_annot=filter_annot,
        plot=plot,
        palette_path=palette_path,
        plot_path=image_store_path,
        **kwargs,
    )
    
    logger.info(f"MMC annotation for experiment {exp_name} complete")
    return result

# ==== I/_utils.py ====
import os
from pathlib import Path
import glob
import warnings
import anndata as ad


def backup_adata_region(
    exp_name: str, reg_name: str, prefix_name: str, adata_path: str = None
):
    """
    Backup function for AnnData objects. Due to env incompatibilities with
    SpatialData and some of the annotation / integration libraries.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        import spatialdata as sd
        from spida._utilities import _backup_adata, _gen_keys, _region_to_donor
        from spida._constants import TABLE_KEY

    print(
        "BACKING UP ADATA, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    _donor = _region_to_donor(reg_name)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # Getting the sdata object (right now from a constant zarr store path)
        zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
        zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
        sdata = sd.read_zarr(zarr_path)

    if not adata_path:
        adata_path = os.getenv("ANNDATA_STORE_PATH")
        if not adata_path:
            raise ValueError(
                "Please provide an adata_path or set the ANNDATA_STORE_PATH environment variable."
            )

    # getting the adata object
    adata = ad.read_h5ad(f"{adata_path}/{exp_name}/{prefix_name}/adata_{_donor}.h5ad")

    _backup_adata(sdata, adata, KEYS[TABLE_KEY])


def backup_adata_experiment(exp_name: str, prefix_name: str, adata_path: str = None):
    """
    Backup function for AnnData objects for an entire experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    for reg in regions:
        reg_name = reg.split("/")[-1]
        backup_adata_region(exp_name, reg_name, prefix_name, adata_path=adata_path)
