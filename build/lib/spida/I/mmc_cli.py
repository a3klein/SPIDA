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
