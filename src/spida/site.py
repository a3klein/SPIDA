### This notebook is meant to generate the necessary files for the SPIDA website. It should be run from the root of the SPIDA repository.
from __future__ import annotations

import os
import sys
from pathlib import Path
import json
import re

import click 
from rich_click import RichCommand, RichGroup, echo as rich_echo # type: ignore
from dotenv import load_dotenv  # type: ignore

from spida.settings import configure_logging_for_runtime
from spida.config import ConfigDefaultGroup
from spida.utilities.site_utils import (
    _load_metadata,
    append_brain_region_qc_metrics,
    generate_soma_gene_proportions,
)

import logging
logger = logging.getLogger(__name__)
# Might need to do this in every root script until spatialdata_plot merges their fix to the main branch
sys.ps1 = ">>> "  # Set the primary prompt string to indicate interactive mode

@click.group(cls=ConfigDefaultGroup)
@click.option(
    "config",
    "--config",
    default=".env",
    type=click.Path(exists=True),
    help="Path to the configuration file (default: .env)",
)
@click.pass_context
def cli(ctx, config):
    ctx.ensure_object(dict)
    pass


@cli.command(cls=RichCommand, help="Generate QC figures for a given sample for the SPIDA website.")
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("--lab", type=str, default=None, help="Name of the lab (optional)")
@click.option("--brain-region", type=str, default="WB", help="Name of the brain region (default: WB)")
@click.option("--plot-decon", is_flag=True, default=False, help="Whether to plot deconvolved images (default: False)")
@click.option("--zarr-store", type=str, default=None, help="Path to the zarr store (optional, defaults to ZARR_STORAGE_PATH env variable)")
@click.option("--site-dir", type=str, default=None, help="Path to the SPIDA site directory (optional, defaults to SPIDA_SITE_DIR env variable)")
@click.pass_context
def generate_load_qc_figs(
    ctx,
    exp_name : str, 
    reg_name : str,
    lab : str = None,
    brain_region : str = "WB",
    zarr_store : str = None,
    site_dir : str | Path = None,
    plot_decon : bool = False
): 
    """
    Generate QC figures for a given experiment and region, and save them to the appropriate directory for the SPIDA website.
    """
    from spida.P.transcript_qc import _obs_to_grid_geodf
    from spida.pl.site_figs import plot_load_images, plot_decon_images, plot_tz_qc, plot_tz_hex_qc

    logger.info(f"Generating load QC figures for experiment {exp_name} region {reg_name}")

    sdata, img_dir, data_dir, keys = _load_metadata(
        exp_name=exp_name,
        reg_name=reg_name,
        lab=lab,
        brain_region=brain_region,
        zarr_store=zarr_store,
        site_dir=site_dir
    )
    
    image_key = keys[0]
    decon_image_key = "decon_image"
    tz_qc_key = "transcript_qc_shapes"
    tz_qc_table_key = "adata_hex_s30.0_o0.0_clustered"
    
    try: 
        grid=sdata[tz_qc_key].copy()
    except KeyError:
        logger.info(f"Key {tz_qc_key} not found in sdata for experiment {exp_name} region {reg_name}. Skipping QC plots.")
        tz_qc_grid_key = "adata_hex_s30.0_o0.0"
        adata_hex = sdata[tz_qc_grid_key].copy()
        adata_hex.obs = _obs_to_grid_geodf(adata_hex.obs)    

    cmap="RdYlGn_r"
    adata_hex = sdata[tz_qc_table_key].copy()
    adata_hex.obs = _obs_to_grid_geodf(adata_hex.obs)    
    
    plot_load_images(sdata, image_key, img_dir)
    if plot_decon:
        plot_decon_images(sdata, decon_image_key, img_dir)
    plot_tz_qc(grid, img_dir, cmap=cmap)
    plot_tz_hex_qc(adata_hex, img_dir)

@cli.command(cls=RichCommand, help="Generate QC figures for a given segmentation for the SPIDA website.")
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--lab", type=str, default=None, help="Name of the lab (optional)")
@click.option("--brain-region", type=str, default="WB", help="Name of the brain region (default: WB)")
@click.option("--zarr-store", type=str, default=None, help="Path to the zarr store (optional, defaults to ZARR_STORAGE_PATH env variable)")
@click.option("--site-dir", type=str, default=None, help="Path to the SPIDA site directory (optional, defaults to SPIDA_SITE_DIR env variable)")
@click.option(
    "--neuron-type-col",
    type=str,
    default=None,
    help="Column in adata.obs to use for neuron type grouping (optional).",
)
@click.pass_context
def generate_seg_qc_figs(
    ctx,
    exp_name : str, 
    reg_name : str,
    prefix_name : str,
    lab : str = None,
    brain_region : str = "WB",
    zarr_store : str = None,
    site_dir : str | Path = None,
    neuron_type_col: str | None = None,
): 
    """
    Generate QC figures for a given segmentation, and save them to the appropriate directory for the SPIDA website.
    """
    from spida.pl.site_figs import plot_segload, plot_seg_qc, plot_cell_cluster, plot_seg_clust_dotplot

    sdata, img_dir, data_dir, keys = _load_metadata(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        lab=lab,
        brain_region=brain_region,
        zarr_store=zarr_store,
        site_dir=site_dir
    )
    logger.info(f"Generating Seg QC figures for experiment {exp_name} region {reg_name}, segmentation {prefix_name}")
    image_key = keys[0]
    shapes_key = keys[1]
    table_key = keys[2]
    points_key = keys[3]

    if shapes_key not in sdata:
        logger.warning(f"Shapes key {shapes_key} not found in sdata for experiment {exp_name} region {reg_name}. Skipping segmentation QC plots for {prefix_name}.")
    else: 
        plot_segload(sdata, image_key, shapes_key, prefix_name, img_dir)
    if table_key not in sdata:
        logger.warning(f"Table key {table_key} not found in sdata for experiment {exp_name} region {reg_name}. Skipping segmentation QC plots for {prefix_name}.")
    else:
        adata_seg = sdata[table_key].copy()
        logger.info("plot_seg_qc")
        plot_seg_qc(adata_seg, prefix_name, img_dir)
        logger.info("append_brain_region_qc_metrics")
        append_brain_region_qc_metrics(
            adata_seg,
            brain_region=brain_region,
            exp_name=exp_name,
            reg_name=reg_name,
            segmentation=prefix_name,
            site_dir=site_dir,
            lab=lab,
            neuron_type_col=neuron_type_col,
        )
    
    filt_table_key = table_key + "_filt"
    if filt_table_key not in sdata:
        logger.warning(f"Table key {filt_table_key} not found in sdata for experiment {exp_name} region {reg_name}. Skipping segmentation QC plots for {prefix_name}.")
    else: 
        adata_seg = sdata[filt_table_key].copy()
        logger.info("plot_cell_cluster")
        plot_cell_cluster(adata_seg, prefix_name, img_dir)
        try: 
            logger.info("plot_seg_clust_dotplot")
            plot_seg_clust_dotplot(adata_seg, prefix_name, img_dir)
        except ValueError as e:
            logger.warning(f"Could not generate dotplot for {prefix_name} in experiment {exp_name} region {reg_name} due to error: {e}")

    if points_key not in sdata: 
        logger.warning(f"Points key {points_key} not found in sdata for experiment {exp_name} region {reg_name}. Skipping soma gene proportion calculation for {prefix_name}.")
    else:
        points = sdata[points_key].compute()
        logger.info("generate_soma_gene_proportions")
        generate_soma_gene_proportions(points, prefix_name, data_dir)


@cli.command(cls=RichCommand, help="Generate data.json manifest for the SPIDA website from the images and data directories.")
@click.option("--site-dir", type=str, default=None, help="Path to the SPIDA site directory (optional, defaults to SPIDA_SITE_DIR env variable)")
@click.option("--experiment-alias-re", type=str, default=r"4x1-([^/]+)-(?:E|Q)", help="Regular expression to extract experiment alias from folder names (default: r'4x1-([^/]+)-(?:E|Q)')")
@click.option("--region-alias-re", type=str, default=r"region_([^/_]+)", help="Regular expression to extract region alias from folder names (default: r'region_([^/_]+)')")
@click.pass_context
def generate_data_manifest(
    ctx,
    site_dir : str | Path = None,
    experiment_alias_re : str = r"4x1-([^/]+)-(?:E|Q)",
    region_alias_re : str = r"region_([^/_]+)",
): 
    """Generate data.json manifest from images and data directories.

    Assumes folder structure:
        images/{brain_region}/{experiment}_{region}
        data/{brain_region}/{experiment}_{region}

    Alias rules: (these may be subject to change as this expands)
    - experiment alias: text between "4x1-" and "-E" or "-Q"
    - region alias (donor): text after "region_"
    """
    if site_dir is None:
        site_dir = os.getenv("SPIDA_SITE_DIR", None)
    if isinstance(site_dir, str):
        site_dir = Path(site_dir)
    IMAGES_DIR = site_dir / "images" # absolutes
    DATA_DIR = site_dir / "data" # absolutes
    OUTPUT_FILE = site_dir / "data.json" # absolutes

    EXPERIMENT_ALIAS_RE = re.compile(experiment_alias_re)
    REGION_ALIAS_RE = re.compile(region_alias_re)

    def experiment_alias(name: str) -> str:
        match = EXPERIMENT_ALIAS_RE.search(name)
        return match.group(1) if match else name

    def region_alias(name: str) -> str:
        match = REGION_ALIAS_RE.search(name)
        return match.group(1) if match else name

    def collect_files(folder: Path, suffix: str | None = None) -> list[str]:
        if not folder.exists():
            return []
        files = [
            p.name
            for p in folder.iterdir()
            if p.is_file() and not p.name.startswith(".") and (suffix is None or p.suffix == suffix)
        ]
        return sorted(files)

    if not IMAGES_DIR.exists():
        logger.warning(f"Images directory not found: {IMAGES_DIR}. Cannot generate data manifest.")
        raise ValueError(f"Images directory not found: {IMAGES_DIR}. Cannot generate data manifest.")

    brain_regions = []
    for region_dir in sorted(p for p in IMAGES_DIR.iterdir() if p.is_dir()):
        region_entry = {
            "name": region_dir.name,
            "alias": region_dir.name,
            "experiments": [],
        }

        for experiment_dir in sorted(p for p in region_dir.iterdir() if p.is_dir()):
            data_dir = DATA_DIR / region_dir.name / experiment_dir.name
            region_entry["experiments"].append(
                {
                    "id": experiment_dir.name,
                    "label": experiment_dir.name,
                    "alias": experiment_alias(experiment_dir.name),
                    "region_alias": region_alias(experiment_dir.name),
                    "image_files": collect_files(experiment_dir),
                    "data_files": collect_files(data_dir, suffix=".csv"),
                }
            )

        brain_regions.append(region_entry)

    output = {"brain_regions": brain_regions}
    OUTPUT_FILE.write_text(json.dumps(output, indent=2), encoding="utf-8")
    logger.info(f"Wrote {OUTPUT_FILE}")

#TODO: for the generate manifest also add a fodler for .csv files in the same layout as images (?) --> This gets annoying for integrating the stuff from multiple labs? 

#TODO: For each experiment extract metadata such as genes/cell ... for each segmentation and save it to one csv shared across all experiment!

cli.add_command(generate_load_qc_figs)
cli.add_command(generate_seg_qc_figs)
cli.add_command(generate_data_manifest)

if __name__ == "__main__":
    
    # Configure root logger with INFO level handlers (to allow INFO messages through)
    env = configure_logging_for_runtime(
        level=logging.INFO,  # Handlers need to accept INFO level
    )
    # Set root logger level to WARNING to suppress other modules
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)
    
    # Set the entire spida package to INFO level as well
    spida_logger = logging.getLogger('spida')
    spida_logger.setLevel(logging.INFO)
    # Set the level for the current logger
    logger.setLevel(logging.INFO) # Set the level for the current logger
    logger.info(f"Logging configured for environment: {env}")
            
    cli()