#!/usr/bin/env python3
"""
CLI interface for SPIDA P module (PreProcessing) using click
"""

from pathlib import Path
import click
from rich_click import RichCommand # type: ignore
import logging

# from spida.utilities.script_utils import parse_click_kwargs, JSONParam
from spida.settings import configure_logging_for_runtime
from spida.config import ConfigDefaultGroup
from spida.io.options import dataset_options, dataset_overrides, dataset_from_options, run_config
from spida.utilities.script_utils import parse_click_kwargs


def _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, **roots):
    return dataset_from_options(run_config(ctx, **roots), exp, reg, zarr_path, raw_path, seg_path, figures_path)


import logging
import warnings

from spida.utilities.script_utils import parse_list

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

@click.group(cls=ConfigDefaultGroup)
@click.option(
    "config",
    "--config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="JSON config file (optional; unset settings come from the environment).",
)
@click.pass_context
def cli(ctx, config):
    ctx.ensure_object(dict)
    pass

# Show how a call's configuration resolves
@click.command(cls=RichCommand, help="Print the resolved configuration: --config JSON, then the environment, then these options.")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Directory holding the VPT binary, --backend vpt only (default: VPT_BIN_PATH in the config)")
@click.option("rust_bin_path", "--rust_bin_path", default=None, type=click.Path(), help="Directory holding the proseg binary (default: RUST_BIN_PATH in the config)")
@click.option("deconwolf_config_file", "--deconwolf_config_file", default=None, type=click.Path(), help="deconwolf .config.ini: binary, PSFs, [Optics] (default: DECONWOLF_CONFIG in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("anndata_store", "--anndata_store", default=None, type=click.Path(), help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.option("annotation_store", "--annotation_store", default=None, type=click.Path(), help="Root of the annotation intermediates (default: ANNOTATION_STORE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.pass_context
def test_command(
    ctx,
    vpt_bin_path,
    rust_bin_path,
    deconwolf_config_file,
    zarr_store,
    root_path,
    segmentation_store,
    anndata_store,
    annotation_store,
    image_store,
):
    """Print every SpidaConfig field as this call resolves it."""
    import dataclasses
    cfg = run_config(ctx, vpt_bin_path=vpt_bin_path, rust_bin_path=rust_bin_path,
                     deconwolf_config_file=deconwolf_config_file, zarr_store=zarr_store, root_path=root_path,
                     segmentation_store=segmentation_store, anndata_store=anndata_store,
                     annotation_store=annotation_store, image_store=image_store)
    for f in dataclasses.fields(cfg):
        click.echo(click.style(f.name, fg='magenta', bold=True) + f": {click.style(str(getattr(cfg, f.name)), fg='green')}")



@cli.command(
    name="filter-cells-region",
    cls=RichCommand,
    aliases=["filter_cells_region"],
    help="Filter cells of the PREFIX_NAME segmentation of one region (--exp/--reg, or the path options).",
    )
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--seg_fam", type=str, default=None, help="Which preset of column names to use when running filtering (default is None which uses vpt outputs)")
@click.option("--cutoffs_path", type=click.Path(exists=True, file_okay=True), default=None, help="QC cutoffs JSON (default: CUTOFFS_PATH in the config)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("--overwrite", is_flag=True, default=False, help="Re-filter a table that was already filtered (default: False)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@dataset_overrides
@click.pass_context
def filter_cells_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    seg_fam:str = None,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Filter cells of the PREFIX_NAME segmentation of one region.
    """
    for key, value in ctx.params.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    cfg = run_config(ctx, zarr_store=zarr_store, image_store=image_store, cutoffs_path=cutoffs_path)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from spida.P.main import filter_cells_region as func
    func(dataset, prefix_name, seg_fam=seg_fam, cutoffs_path=cfg.cutoffs_path,
         plot=plot, image_path=image_path, overwrite=overwrite)

@cli.command(
    name="filter-cells-all",
    cls=RichCommand,
    aliases=["filter_cells_all"],
    help="Filter cells of the PREFIX_NAME segmentation for every region of --exp.",
)
@click.argument("prefix_name", type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("--seg_fam", type=str, default=None, help="Which preset of column names to use when running filtering (default is None which uses vpt outputs)")
@click.option("--cutoffs_path", type=click.Path(exists=True), default=None, help="QC cutoffs JSON (default: CUTOFFS_PATH in the config)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--overwrite", is_flag=True, default=False, help="Re-filter regions already filtered; without it they are skipped (default: False)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.pass_context
def filter_cells_all(
    ctx,
    prefix_name: str,
    exp_name: str,
    cutoffs_path: Path = None,
    seg_fam:str = None,
    plot: bool = False,
    overwrite: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Filter cells of the PREFIX_NAME segmentation for every region of an experiment.
    """
    cfg = run_config(ctx, zarr_store=zarr_store, image_store=image_store, cutoffs_path=cutoffs_path)
    from spida.P.main import filter_cells_all as func
    func(cfg, exp_name, prefix_name, seg_fam=seg_fam, cutoffs_path=cfg.cutoffs_path, plot=plot, overwrite=overwrite)

@cli.command(
    name="plot-filtering-region",
    cls=RichCommand,
    aliases=["plot_filtering_region"],
    help="Plot the filtering results of the PREFIX_NAME segmentation of one region.",
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@dataset_overrides
@click.pass_context
def plot_filtering_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    image_store: Path = None,
    image_path: Path = None,
    suffix: str = "",
    zarr_store: Path = None,
):
    """
    Plot the filtering results of the PREFIX_NAME segmentation of one region.
    """
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import plot_filtering_region as func
    func(dataset, prefix_name, image_path=image_path, suffix=suffix)

@cli.command(
    name="write-adata",
    cls=RichCommand,
    aliases=["write_adata"],
    help="Export segmentation tables of --exp (one --reg, or all regions) to h5ad.",
)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("--reg", "reg_name", type=str, default=None, help="Region name (default: every region of --exp)")
@click.option("--prefix_names", type=parse_list, required=True, help="List of prefix names")
@click.option("--suffix", type=str, default="", help="Table suffix, e.g. _filt (default: none)")
@click.option("anndata_store", "--anndata_store", default=None, type=click.Path(), help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Replace existing h5ad copies; without it, one region raises and all regions skip existing files (default: False)")
@click.pass_context
def write_adata(
    ctx,
    exp_name,
    reg_name: str = None,
    prefix_names: list[str] = None,
    suffix: str = "",
    anndata_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Export segmentation tables to AnnData files, ANNDATA_STORE_PATH/<exp>/<reg>/<prefix>_table<suffix>.h5ad.
    """
    cfg = run_config(ctx, zarr_store=zarr_store, anndata_store=anndata_store)
    from spida.P.main import write_adata as func
    func(cfg, exp_name, reg_name=reg_name, prefix_names=prefix_names, suffix=suffix, overwrite=overwrite)

@cli.command(
    name="setup-adata-region",
    cls=RichCommand,
    aliases=["setup_adata_region"],
    help="Setup the PREFIX_NAME table of one region for downstream analysis (--exp/--reg, or the path options).",
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--suffix", type=str, default="", help="Suffix for the table key (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--do_tsne", is_flag=True, default=False, help="Also compute the tSNE embedding (needed by the setup plots; default: False = UMAP only)")
@click.option("--consensus/--no-consensus","consensus_cluster",default=True,help="Consensus clustering (200 Leiden repeats) vs a single Leiden run. ")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Rebuild PREFIX_NAME's table+suffix if it exists; drops doublet / RESOLVI results on it (default: False)")
@dataset_overrides
@click.pass_context
def setup_adata_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    suffix: str = "",
    plot=False,
    do_tsne: bool = False,
    consensus_cluster: bool = True,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Setup the AnnData object for downstream analysis.
    This involves normalizing data, calculating PCA, umap, tsne, and leiden clusters
    """
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import setup_adata_region as func
    func(dataset, prefix_name, suffix=suffix, plot=plot, do_tsne=do_tsne,
         consensus_cluster=consensus_cluster, image_path=image_path, overwrite=overwrite)


@cli.command(
    name="setup-adata-all",
    cls=RichCommand,
    aliases=["setup_adata_all"],
    help="Setup the PREFIX_NAME tables for every region of --exp.",
)
@click.argument("prefix_name", type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("--suffix", type=str, default="", help="Suffix for the output table key in the spatialdata object (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Rebuild existing tables; without it, those regions are skipped (default: False)")
@click.pass_context
def setup_adata_all(
    ctx,
    prefix_name: str,
    exp_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Setup the AnnData objects for every region of an experiment.
    """
    cfg = run_config(ctx, zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import setup_adata_all as func
    func(cfg, exp_name, prefix_name, suffix=suffix, plot=plot, image_path=image_path, overwrite=overwrite)

@cli.command(
    name="plot-setup-region",
    cls=RichCommand,
    aliases=["plot_setup_region"],
    help="Plot the setup results of the PREFIX_NAME table of one region.",
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@dataset_overrides
@click.pass_context
def plot_setup_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    image_path: Path = None,
    image_store: Path = None,
    suffix: str = "",
    zarr_store: Path = None,
):
    """
    Plot the setup results of the PREFIX_NAME table of one region.
    """
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import plot_setup_region as func
    func(dataset, prefix_name, image_path=image_path, suffix=suffix)

@cli.command(
    name="remove-doublets-region",
    cls=RichCommand,
    aliases=["remove_doublets_region"],
    help="Remove doublets from the AnnData object for a specific region",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("--threshold", type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Re-run on a table that already has doublet calls (default: False)")
@dataset_overrides
@click.pass_context
def remove_doublets_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Remove doublets from the AnnData object for a specific region in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import remove_doublets_region as func
    func(dataset, prefix_name, threshold=threshold, suffix=suffix, plot=plot, model_kwargs=kwargs, overwrite=overwrite)

@cli.command(
    name="remove-doublets-all",
    cls=RichCommand,
    aliases=["remove_doublets_all"],
    help="Remove doublets from all AnnData objects in an experiment, for a given segmentation given by PREFIX_NAME",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("prefix_name", type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("--threshold", type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Re-run regions that already have doublet calls; without it, they are skipped (default: False)")
@click.pass_context
def remove_doublets_all(
    ctx,
    prefix_name: str,
    exp_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Remove doublets from all regions in an experiment
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    cfg = run_config(ctx, zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import remove_doublets_all as func
    func(cfg, exp_name, prefix_name, threshold=threshold, suffix=suffix, plot=plot, model_kwargs=kwargs, overwrite=overwrite)

@cli.command(
    name="resolvi-region",
    cls=RichCommand,
    aliases=["resolvi_region"],
    help="Perform RESOLVI clustering on the AnnData object for a specific region in an experiment",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--suffix", type=str, default="", help="Suffix for the AnnData object (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--model_save_path", type=click.Path(exists=True), default=None, help="Path to save the ResolVI model")
@click.option("--trained", type=bool, default=False, help="Whether the model is already trained (default: False)")
@click.option("--max_epochs", type=int, default=100, help="Maximum number of epochs for training (default: 100)")
@click.option("--layer", type=str, default="raw", help="Layer of the AnnData object to use (default: 'raw')")
@click.option("--batch_key", type=str, default=None, help="Key for batch information in the AnnData object (default: None)")
@click.option("--categorical_covariates", type=parse_list, default=None, help="List of categorical covariates to use (default: [])")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Re-run on a table that already has RESOLVI results (default: False)")
@dataset_overrides
@click.pass_context
def resolvi_region(
    ctx,
    prefix_name: str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    model_save_path: Path = None,
    trained: bool = False,
    max_epochs: int = 100,
    layer: str = "raw",
    batch_key: str = None,
    categorical_covariates: list = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Perform RESOLVI clustering on the AnnData object for a specific region in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import resolvi_cluster_region as func
    func(
        dataset,
        prefix_name,
        suffix=suffix,
        plot=plot,
        model_save_path=model_save_path,
        trained=trained,
        max_epochs=max_epochs,
        layer=layer,
        batch_key=batch_key,
        categorical_covariates=categorical_covariates,
        model_kwargs=kwargs,
        overwrite=overwrite,
    )

@cli.command(
    name="resolvi-all",
    cls=RichCommand,
    aliases=["resolvi_all"],
    help="Run the Resolvi algorithm for all regions in an experiment",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("prefix_name", type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("--suffix", type=str, default="_filt", help="Suffix for the table key in the spatialdata object (default: '_filt')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--model_save_path", type=click.Path(), default=None, help="Path to save the ResolVI model")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Re-run regions that already have RESOLVI results; without it, they are skipped (default: False)")
@click.pass_context
def resolvi_all(
    ctx,
    prefix_name: str,
    exp_name: str,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    model_save_path: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    """
    Run the Resolvi algorithm for all regions in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")
    
    cfg = run_config(ctx, zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import resolvi_cluster_all as func
    func(cfg, exp_name, prefix_name, suffix=suffix, plot=plot,
         model_save_path=model_save_path, model_kwargs=kwargs, overwrite=overwrite)

@cli.command(
    name="plot-resolvi-region",
    cls=RichCommand,
    aliases=["plot_resolvi_region"],
    help="Plot the resolvi results of the PREFIX_NAME table of one region.",
)
@click.argument("prefix_name", type=str)
@dataset_options
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@dataset_overrides
@click.pass_context
def plot_resolvi_region(
    ctx,
    prefix_name:str,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    image_store:Path = None,
    suffix:str = "",
    zarr_store: Path = None,
):
    """
    Plot the resolvi results of the PREFIX_NAME table of one region.
    """
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import plot_resolvi_region as func
    func(dataset, prefix_name, suffix=suffix)

@cli.command(
    name="combine-datasets",
    cls=RichCommand,
    aliases=["combine_datasets"],
    help="Combine multiple datasets into a single SpatialData object + AnnData object",
)
@click.argument("experiment_list", type=parse_list)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="_filt", help="Suffix to append to the table keys (default: '_filt')")
@click.option("--lab_name", type=str, default="salk", help="Name of the lab (default: 'salk')")
@click.option("--project_name", type=str, default="BG_salk", help="Name of the project (default: 'BG_salk')")
@click.option("--zarr_store", type=click.Path(exists=True), default=None, help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.option("--table_keys", type=parse_list, default=None, help="List of table keys to include (default: None, which uses all keys with the specified suffix)")
@click.pass_context
def combine_datasets(
    ctx,
    experiment_list: list[str],
    prefix_name: str,
    suffix: str = "_filt",
    lab_name: str = "salk",
    project_name: str = "BG_salk",
    zarr_store: Path = None,
    anndata_store: Path = None,
    table_keys: list[str] = None,
):
    """
    Aggregate spatial data from multiple experiments into a single SpatialData object.
    """
    cfg = run_config(ctx, zarr_store=zarr_store, anndata_store=anndata_store)
    from spida.P.main import combine_datasets as func
    func(
        experiment_list=experiment_list,
        prefix_name=prefix_name,
        suffix=suffix,
        lab_name=lab_name,
        project_name=project_name,
        config=cfg,
        table_keys=table_keys,
    )

@cli.command(
    name="setup-dataset",
    cls=RichCommand,
    aliases=["setup_dataset"],
    help="Perform standard PCA + Leiden clustering on a joint dataset object",
)
@click.argument("dataset_name", type=str)
@click.option("--scale", type=bool, default=False, help="Whether to scale the data (default: False)")
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.pass_context
def setup_dataset(
    ctx,
    dataset_name: str,
    scale: bool = False,
    anndata_store: Path = None,
):
    """
    Perform standard PCA + Leiden clustering on a joint dataset object.
    """
    cfg = run_config(ctx, anndata_store=anndata_store)
    from spida.P.main import setup_dataset as func
    func(dataset_name=dataset_name, scale=scale, anndata_store=cfg.anndata_store)

@cli.command(
    name="resolvi-dataset",
    cls=RichCommand,
    aliases=["resolvi_dataset"],
    help="Perform ResolVI workflow on a joint dataset object",
)
@click.argument("dataset_name", type=str)
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.option("--model_save_path", type=click.Path(), default=None, help="Path to save the ResolVI model")
@click.option("--trained", type=bool, default=False, help="Whether the model is already trained (default: False)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--max_epochs", type=int, default=100, help="Maximum number of epochs for training (default: 100)")
@click.option("--categorical_covariates", type=parse_list, default=['donor', 'brain_region'], help="List of categorical covariates to use (default: ['donor', 'brain_region'])")
@click.option("--batch_key", type=str, default="dataset_id", help="Key for batch information in the AnnData object (default: 'dataset_id')")
@click.option("--layer", type=str, default="counts", help="Layer of the AnnData object to use (default: 'counts')")
@click.pass_context
def resolvi_dataset(
    ctx,
    dataset_name:str,
    anndata_store:Path = None,
    model_save_path:Path = None, 
    trained:bool = False,
    image_store:Path = None,
    max_epochs:int = 100,
    layer:str = "counts", 
    batch_key:str = "dataset_id",
    categorical_covariates:list = ['donor', 'brain_region'],
    ): 
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    cfg = run_config(ctx, anndata_store=anndata_store, image_store=image_store)
    from spida.P.main import resolvi_dataset as func
    func(
        dataset_name=dataset_name,
        anndata_store=cfg.anndata_store,
        image_store=cfg.image_store,
        model_save_path=model_save_path,
        trained=trained,
        max_epochs=max_epochs,
        layer=layer,
        batch_key=batch_key,
        categorical_covariates=categorical_covariates,
        model_kwargs=kwargs,
    )


@cli.command(
    name="plot-dataset",
    cls=RichCommand,
    aliases=["plot_dataset"],
    help="Plot the dataset resolvi + setup results ",
)
@click.argument("dataset_name", type=str)
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("--show", is_flag=True, help="Whether to show the plot (default: False)")
@click.pass_context
def plot_dataset(
    ctx,
    dataset_name : str,
    anndata_store: str | Path = None,
    image_store: str | Path | None = None,
    show: bool = False
):
    """
    Plot the dataset results for a specific region in an experiment.
    """
    cfg = run_config(ctx, anndata_store=anndata_store, image_store=image_store)
    from spida.P.main import plot_dataset_setup as func
    func(dataset_name=dataset_name, anndata_store=cfg.anndata_store, image_store=cfg.image_store, show=show)

@cli.command(
    name="call-region-tz",
    cls=RichCommand,
    aliases=["call_region_tz"],
    help="Call the transcript-based spatial regions for a specific region in an experiment",
)
@click.argument("prefix_name", type=str, default="default")
@click.argument("geoms_name", type=str, default="wm_region")
@dataset_options
@click.option("--use_genes", type=parse_list, default="BCAS1", help="List of genes to use for calling regions (default: ['BCAS1'])")
@click.option("--save_geoms_path", type=click.Path(), default=None, help="Path to save the geometries (default: None)")
@click.option("--dsc_comp_min_size", type=int, default=5, help="Minimum size of the DSC components (default: 5)")
@click.option("--hex_size", type=int, default=50, help="Size of the hexagons (default: 50)")
@click.option("--hex_overlap", type=int, default=0, help="Overlap between hexagons (default: 0)")
@click.option("--gmm_ncomp", type=str, default="auto", help="Number of GMM components, 'auto calculates the best number of components based on the BIC (default: 'auto')")
@click.option("--gmm_cov_type", type=str, default="full", help="Covariance type for GMM (default: 'full')")
@click.option("--gene_agreement_thr", type=float, default=0.25, help="Threshold for gene agreement (default: 0.25)")
@click.option("--top_n_comp", type=int, default=1, help="Top N components to consider (default: 1)")
@click.option("--plot", is_flag=True, help="Whether to generate plots (default: False)")
@click.option("--image_path", type=click.Path(), default=None, help="Path to save the plots (default: None)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Replace existing GEOMS_NAME shapes (and --save_geoms_path) (default: False)")
@dataset_overrides
@click.pass_context
def call_region_tz(
    ctx,
    prefix_name : str = "default",
    geoms_name : str = "wm_region",
    exp=None, reg=None, zarr_path=None, raw_path=None, seg_path=None, figures_path=None,
    use_genes : str | list[str] = "BCAS1",
    save_geoms_path : str | None = None,
    dsc_comp_min_size : int = 5,
    hex_size : int = 50,
    hex_overlap : int = 0,
    gmm_ncomp : int | str = "auto",
    gmm_cov_type : str = "full",
    gene_agreement_thr : float = 0.25,
    top_n_comp : int = 1,
    plot : bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None,
    overwrite: bool = False,
): 
    """
    Call the transcript-based spatial regions for a specific region in an experiment.
    """
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import call_region_tz as func
    func(
        dataset,
        prefix_name=prefix_name,
        geoms_name=geoms_name,
        use_genes=use_genes,
        save_geoms_path=save_geoms_path,
        dsc_comp_min_size=dsc_comp_min_size,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gmm_ncomp=gmm_ncomp,
        gmm_cov_type=gmm_cov_type,
        gene_agreement_thr=gene_agreement_thr,
        top_n_comp=top_n_comp,
        gen_plots=plot,
        image_path=image_path,
        overwrite=overwrite,
    )

@cli.command(
    name="transcript-qc",
    cls=RichCommand,
    aliases=["transcript_qc"],
    help="Run transcript hex QC for one region (--exp/--reg, or --zarr_path).",
)
@dataset_options
@click.option("--prefix_name", type=str, default="default")
@click.option("--qc_shapes_key", type=str, default="transcript_qc_shapes", help="Key to store the QC shapes in the spatialdata object (default: 'transcript_qc_shapes')")
@click.option("--hex_size", type=float, default=30, help="Hex size (distance center->vertex) (default: 30)")
@click.option("--hex_overlap", type=float, default=0, help="Hex overlap (default: 0)")
@click.option("--gene_col", type=str, default="gene", help="Gene column name (default: 'gene')")
@click.option("--x_col", type=str, default="x", help="X coordinate column (default: 'x')")
@click.option("--y_col", type=str, default="y", help="Y coordinate column (default: 'y')")
@click.option("--min_transcripts", type=int, default=100, help="Minimum transcripts per hex (default: 100)")
@click.option("--min_density", type=float, default=None, help="Minimum density (per µm²) (default: None)")
@click.option("--plot", is_flag=True, help="Whether to plot the results (default: False)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Replace existing QC shapes; the shared hex table is reused either way (default: False)")
@dataset_overrides
@click.pass_context
def transcript_qc(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    prefix_name: str = "default",
    qc_shapes_key: str = "transcript_qc_shapes",
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int = 100,
    min_density: float = None,
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import transcript_qc_region as func
    func(
        dataset,
        prefix_name=prefix_name,
        qc_shapes_key=qc_shapes_key,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        min_transcripts=min_transcripts,
        min_density=min_density,
        plot=plot,
        overwrite=overwrite,
    )

@cli.command(
    name="tissue-mask",
    cls=RichCommand,
    aliases=["tissue_mask"],
    help="Compute the KDE tissue mask for one region: a polygon and a label layer in the zarr store, and tissue_mask.tif.",
)
@dataset_options
@click.option("--prefix_name", type=str, default="default")
@click.option("--gene_col", type=str, default="gene", help="Gene column name (default: 'gene')")
@click.option("--x_col", type=str, default="x", help="X coordinate column (default: 'x')")
@click.option("--y_col", type=str, default="y", help="Y coordinate column (default: 'y')")
@click.option("--rho_min", type=float, default=5.0, help="Min transcript density tx/100µm² for tissue (default: 5)")
@click.option("--ref_um", type=float, default=25.0, help="Fixed reference scale (µm) for the density estimate (default: 25)")
@click.option("--erode_sigma", type=float, default=1.0, help="Boundary erosion in reference-sigma units (default: 1)")
@click.option("--a_min_um2", type=float, default=50000.0, help="Drop connected components smaller than this µm² (default: 50000)")
@click.option("--gb", type=float, default=5.0, help="Micron grid spacing (default: 5)")
@click.option("--shapes_key", type=str, default="qc_mask", help="Key for the pass-region polygon shapes (default: 'qc_mask')")
@click.option("--labels_key", type=str, default="tissue_mask", help="Key for the boolean labels layer (default: 'tissue_mask')")
@click.option("--no_raster", is_flag=True, default=False, help="Skip writing the pixel tissue_mask.tif")
@click.option("--overwrite", is_flag=True, default=False, help="Replace an existing qc_mask, tissue_mask and tissue_mask.tif (default: False)")
@click.option("--tissue_mask_output_path", "tissue_mask_output_path", type=click.Path(file_okay=False), default=None, help="Directory to write tissue_mask.tif into (default: the region's raw images dir, where segment-region looks for it)")
@click.option("--plot", is_flag=True, help="Whether to plot the mask overview (default: False)")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@dataset_overrides
@click.pass_context
def tissue_mask(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    prefix_name: str = "default",
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    rho_min: float = 5.0,
    ref_um: float = 25.0,
    erode_sigma: float = 1.0,
    a_min_um2: float = 50000.0,
    gb: float = 5.0,
    shapes_key: str = "qc_mask",
    labels_key: str = "tissue_mask",
    no_raster: bool = False,
    overwrite: bool = False,
    tissue_mask_output_path: Path = None,
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
    root_path: Path = None,
):
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store, root_path=root_path)
    from spida.P.main import tissue_mask_region as func
    func(
        dataset,
        prefix_name=prefix_name,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        rho_min=rho_min,
        ref_um=ref_um,
        erode_sigma=erode_sigma,
        a_min_um2=a_min_um2,
        gb=gb,
        shapes_key=shapes_key,
        labels_key=labels_key,
        write_raster=not no_raster,
        plot=plot,
        mask_dir=tissue_mask_output_path,
        overwrite=overwrite,
    )

@cli.command(
    name="cluster-hexes",
    cls=RichCommand,
    aliases=["cluster_hexes"],
    help="Run hex clustering for one region (--exp/--reg, or --zarr_path).",
)
@dataset_options
@click.option("--prefix_name", type=str, default="default")
@click.option("--hex_size", type=float, default=30, help="Hex size (distance center->vertex)")
@click.option("--hex_overlap", type=float, default=0, help="Hex overlap")
@click.option("--gene_col", type=str, default="gene", help="Gene column name")
@click.option("--x_col", type=str, default="x", help="X coordinate column")
@click.option("--y_col", type=str, default="y", help="Y coordinate column")
@click.option("--min_transcripts", type=int, default=100, help="Minimum transcripts per hex")
@click.option("--min_density", type=float, default=None, help="Minimum density (per µm²)")
@click.option("--leiden_resolution", type=float, default=1.0, help="Leiden resolution")
@click.option("--min_cells", type=int, default=10, help="Min cells per gene")
@click.option("--min_genes", type=int, default=100, help="Min genes per hex")
@click.option("--n_top_genes", type=int, default=200, help="Number of HVGs")
@click.option("--pca_comps", type=int, default=50, help="PCA components")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(), default=None, help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("--overwrite", is_flag=True, default=False, help="Replace the existing clustered hex table; the shared hex table is reused either way (default: False)")
@dataset_overrides
@click.pass_context
def cluster_hexes(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    prefix_name: str = "default",
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int = 100,
    min_density: float = None,
    leiden_resolution: float = 1.0,
    min_cells: int = 10,
    min_genes: int = 100,
    n_top_genes: int = 200,
    pca_comps: int = 50,
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
    overwrite: bool = False,
):
    dataset = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                       zarr_store=zarr_store, image_store=image_store)
    from spida.P.main import cluster_hexes_region as func
    func(
        dataset,
        prefix_name=prefix_name,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        min_transcripts=min_transcripts,
        min_density=min_density,
        leiden_resolution=leiden_resolution,
        min_cells=min_cells,
        min_genes=min_genes,
        n_top_genes=n_top_genes,
        pca_comps=pca_comps,
        plot=plot,
        overwrite=overwrite,
    )

cli.add_command(test_command)
cli.add_command(filter_cells_region)
cli.add_command(filter_cells_all)
cli.add_command(setup_adata_region)
cli.add_command(setup_adata_all)
cli.add_command(remove_doublets_region)
cli.add_command(remove_doublets_all)
cli.add_command(resolvi_region)
cli.add_command(resolvi_all)
cli.add_command(combine_datasets)
cli.add_command(setup_dataset)
cli.add_command(resolvi_dataset)
cli.add_command(write_adata)
cli.add_command(plot_filtering_region)
cli.add_command(plot_setup_region)
cli.add_command(plot_resolvi_region)
cli.add_command(plot_dataset)
cli.add_command(call_region_tz)
cli.add_command(transcript_qc)
cli.add_command(cluster_hexes)

def main(): 
    """Main entry point for the CLI."""
    
    # TODO: Functionalize this configuration setup? 
    # - Technically once SPIDA is fully built this will move into the main entry point of the package.
    # - and this won't be done once in each module! 

     #LOGGING
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
            
    # CALL
    cli()

if __name__ == "__main__":
    main()
