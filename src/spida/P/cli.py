#!/usr/bin/env python3
"""
CLI interface for SPIDA P module (PreProcessing) using click
"""

import os
from pathlib import Path
import click
from rich_click import RichCommand, RichGroup # type: ignore
import logging

# from spida.utilities.script_utils import parse_click_kwargs, JSONParam
from spida.settings import configure_logging_for_runtime
from spida.config import ConfigDefaultGroup
from spida.utilities.script_utils import parse_click_kwargs, JSONParam


import sys
import argparse
import logging
import inspect
import warnings

from spida.utilities.script_utils import ParseKwargs, parse_path, parse_list

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

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

# For Testing the CLI Group default_map functionality
@click.command(cls=RichCommand, help="test_command")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Path to the VPT binary")
@click.option("rust_bin_path", "--rust_bin_path", default=None, type=click.Path(), help="Path to the Rust binary")
@click.option("deconwolf_config_file", "--deconwolf_config_file", default=None, type=click.Path(), help="Path to the Deconwolf config file")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Path to the processed root directory")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Path to the segmentation output directory")
@click.option("anndata_store", "--anndata_store", default=None, type=click.Path(), help="Path to the AnnData store directory")
@click.option("annotation_store", "--annotation_store", default=None, type=click.Path(), help="Path to the annotation store directory")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to the image store directory")
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
    """
    Printing all of the default values for the configuration parameters
    """
    click.echo(click.style("Default configuration values:", bold=True))
    for key, value in ctx.default_map.items():
        click.echo(click.style(key, fg='magenta', bold=True) + f": {click.style(value, fg='green')}")

    click.echo(click.style("\nPassed In Values:\n", bold=True))
    click.echo(click.style(f"vpt_bin_path: {vpt_bin_path}", fg='green'))
    click.echo(click.style(f"rust_bin_path: {rust_bin_path}", fg='green'))
    click.echo(click.style(f"deconwolf_config_file: {deconwolf_config_file}", fg='green'))
    click.echo(click.style(f"zarr_store: {zarr_store}", fg='green'))
    click.echo(click.style(f"root_path: {root_path}", fg='green'))
    click.echo(click.style(f"segmentation_store: {segmentation_store}", fg='green'))
    click.echo(click.style(f"anndata_store: {anndata_store}", fg='green'))
    click.echo(click.style(f"annotation_store: {annotation_store}", fg='green'))
    click.echo(click.style(f"image_store: {image_store}", fg='green'))  



@cli.command(
    name="filter-cells-region",
    cls=RichCommand,
    aliases=["filter_cells_region"],
    help="Filter cells for a specific PREFIX_NAME segmentation for REG_NAME in an experiment EXP_NAME.",
    )
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--seg_fam", type=str, default=None, help="Which preset of column names to use when running filtering (default is None which uses vpt outputs)")
@click.option("--cutoffs_path", type=click.Path(exists=True, file_okay=True), default=None, help="Path to the cutoffs JSON file")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to plot store")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def filter_cells_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    seg_fam:str = None,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Filter cells for a specific PREFIX_NAME segmentation for REG_NAME in an experiment EXP_NAME.
    """
    for key, value in ctx.params.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from spida.P.main import filter_cells_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        seg_fam=seg_fam,
        cutoffs_path=cutoffs_path,
        plot=plot,
        image_path=image_path,
        image_store=image_store,
        zarr_store=zarr_store,
    )

@cli.command(
    name="filter-cells-all",
    cls=RichCommand,
    aliases=["filter_cells_all"],
    help="Filter cells for all regions in an experiment EXP_NAME for a specific PREFIX_NAME segmentation.",
)
@click.argument("exp_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--seg_fam", type=str, default=None, help="Which preset of column names to use when running filtering (default is None which uses vpt outputs)")
@click.option("--cutoffs_path", type=click.Path(exists=True), default=None, help="Path to the cutoffs JSON file")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def filter_cells_all(
    ctx,
    exp_name: str,
    prefix_name: str,
    cutoffs_path: Path = None,
    seg_fam:str = None,
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Filter cells for all regions in an experiment EXP_NAME for a specific PREFIX_NAME segmentation.
    """
    from spida.P.main import filter_cells_all as func
    func(
        exp_name=exp_name,
        prefix_name=prefix_name,
        seg_fam=seg_fam,
        cutoffs_path=cutoffs_path,
        plot=plot,
        image_path=image_store,
        zarr_store=zarr_store,
    )

@cli.command(
    name="plot-filtering-region",
    cls=RichCommand,
    aliases=["plot_filtering_region"],
    help="Plot the filtering results for a specific region in an experiment.",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def plot_filtering_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_store: Path = None,
    image_path: Path = None,
    suffix: str = "",
    zarr_store: Path = None,
):
    """
    Plot the filtering results for a specific region in an experiment.
    """
    from spida.P.main import plot_filtering_region as func  
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        image_path=image_path,
        image_store=image_store,
        suffix=suffix,
        zarr_store=zarr_store,
    )

@cli.command(
    name="write-adata",
    cls=RichCommand,
    aliases=["write_adata"],
    help="Write AnnData objects to disk for a specific experiment and region",
)
@click.argument("exp_name", type=str)
@click.option("--reg_name", type=str, default=None, help="Name of the region (optional, if not provided writes all regions)")
@click.option("--prefix_names", type=parse_list, required=True, help="List of prefix names")
@click.option("--output_path", type=click.Path(exists=True), default=None, help="Output directory path (default is None)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def write_adata(
    ctx,
    exp_name,
    reg_name: str = None,
    prefix_names: list[str] = None,
    output_path: Path = "/ceph/cephatlas/aklein/bican/data/anndatas/",
    zarr_store: Path = None,
):
    """
    Write AnnData objects to disk for a specific experiment and region
    """
    from spida.P.main import write_adata as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_names=prefix_names,
        output_path=output_path,
        zarr_store=zarr_store,
    )

@cli.command(
    name="setup-adata-region",
    cls=RichCommand,
    aliases=["setup_adata_region"],
    help="Setup the anndata object for downstream analysis for a specific region",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="", help="Suffix for the table key (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def setup_adata_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str = "",
    plot=False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Setup the AnnData object for downstream analysis.
    This involves normalizing data, calculating PCA, umap, tsne, and leiden clusters

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """
    from spida.P.main import setup_adata_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        suffix=suffix,
        plot=plot,
        image_path=image_path,
        image_store=image_store,
        zarr_store=zarr_store,
    )


@cli.command(
    name="setup-adata-all",
    cls=RichCommand,
    aliases=["setup_adata_all"],
    help="Setup the AnnData objects for all regions in an experiment",
)
@click.argument("exp_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="", help="Suffix for the output table key in the spatialdata object (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def setup_adata_all(
    ctx,
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Setup the AnnData objects for all regions in an experiment.
    This function iterates over all regions and calls setup_adata for each region.
    """
    from spida.P.main import setup_adata_all as func
    func(
        exp_name=exp_name,
        prefix_name=prefix_name,
        suffix=suffix,
        plot=plot,
        image_path=image_path,
        image_store=image_store,
        zarr_store=zarr_store,
    )

@cli.command(
    name="plot-setup-region",
    cls=RichCommand,
    aliases=["plot_setup_region"],
    help="Plot the setup results for a specific region",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--image_path", type=click.Path(exists=True), default=None, help="Path to store the image")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def plot_setup_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    image_store: Path = None,
    suffix: str = "",
    zarr_store: Path = None,
):
    """
    Plot the setup results for a specific region in an experiment.
    This function generates a PDF with the setup plots
    """
    from spida.P.main import plot_setup_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        image_path=image_path,
        image_store=image_store,
        suffix=suffix,
        zarr_store=zarr_store,
    )

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
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("--threshold", type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def remove_doublets_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Remove doublets from the AnnData object for a specific region in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from spida.P.main import remove_doublets_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        threshold=threshold,
        suffix=suffix,
        plot=plot,
        image_path=image_store,
        model_kwargs=kwargs,
        zarr_store=zarr_store,
    )

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
@click.argument("exp_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("--threshold", type=float, default=0.5, help="Threshold for doublet detection (default: 0.5)")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def remove_doublets_all(
    ctx,
    exp_name: str,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    zarr_store: Path = None,
):
    """
    Remove doublets from all regions in an experiment
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from spida.P.main import remove_doublets_all as func
    func(
        exp_name=exp_name,
        prefix_name=prefix_name,
        threshold=threshold,
        suffix=suffix,
        plot=plot,
        image_path=image_store,
        zarr_store=zarr_store,
        model_kwargs=kwargs,
    )

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
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="", help="Suffix for the AnnData object (default: '')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to save the plot")
@click.option("--model_save_path", type=click.Path(exists=True), default=None, help="Path to save the ResolVI model")
@click.option("--trained", type=bool, default=False, help="Whether the model is already trained (default: False)")
@click.option("--max_epochs", type=int, default=100, help="Maximum number of epochs for training (default: 100)")
@click.option("--layer", type=str, default="raw", help="Layer of the AnnData object to use (default: 'raw')")
@click.option("--batch_key", type=str, default=None, help="Key for batch information in the AnnData object (default: None)")
@click.option("--categorical_covariates", type=parse_list, default=None, help="List of categorical covariates to use (default: [])")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def resolvi_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
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
):
    """
    Perform RESOLVI clustering on the AnnData object for a specific region in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from spida.P.main import resolvi_cluster_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        suffix=suffix,
        plot=plot,
        image_path=image_store,
        model_save_path=model_save_path,
        trained=trained,
        max_epochs=max_epochs,
        layer=layer,
        batch_key=batch_key,
        categorical_covariates=categorical_covariates,
        zarr_store=zarr_store,
        model_kwargs=kwargs,
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
@click.argument("exp_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--suffix", type=str, default="_filt", help="Suffix for the table key in the spatialdata object (default: '_filt')")
@click.option("--plot", is_flag=True, help="Whether to plot the results")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("--model_save_path", type=click.Path(), default=None, help="Path to save the ResolVI model")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def resolvi_all(
    ctx,
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_store: Path = None,
    model_save_path: Path = None,
    zarr_store: Path = None,
):
    """
    Run the Resolvi algorithm for all regions in an experiment.
    """
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")
    
    from spida.P.main import resolvi_cluster_all as func
    func(
        exp_name=exp_name,
        prefix_name=prefix_name,
        suffix=suffix,
        plot=plot,
        image_store=image_store,
        model_save_path=model_save_path,
        zarr_store=zarr_store,
        **kwargs,
    )

@cli.command(
    name="plot-resolvi-region",
    cls=RichCommand,
    aliases=["plot_resolvi_region"],
    help="Plot the resolvi results for a specific region",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("prefix_name", type=str)
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to the image store")
@click.option("--suffix", type=str, default="", help="Suffix for the keys in the spatialdata object (default: '')")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage")
@click.pass_context
def plot_resolvi_region(
    ctx,
    exp_name:str,
    reg_name:str,
    prefix_name:str,
    image_store:Path = None,
    suffix:str = "",
    zarr_store: Path = None,
):
    """
    Plot the setup results for a specific region in an experiment.
    This function generates a PDF with the setup plots.
    """
    from spida.P.main import plot_resolvi_region as func

    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        image_path=image_store,
        suffix=suffix,
        zarr_store=zarr_store,
    )

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
@click.option("--zarr_store", type=click.Path(exists=True), default=None, help="Path to the Zarr directory")
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Path to the output directory")
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
    from spida.P.main import combine_datasets as func
    func(
        experiment_list=experiment_list,
        prefix_name=prefix_name,
        suffix=suffix,
        lab_name=lab_name,
        project_name=project_name,
        zarr_store=zarr_store,
        anndata_store=anndata_store,
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
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Path to the AnnData directory")
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
    from spida.P.main import setup_dataset as func
    func(
        dataset_name=dataset_name,
        scale=scale,
        anndata_store=anndata_store,
    )

@cli.command(
    name="resolvi-dataset",
    cls=RichCommand,
    aliases=["resolvi_dataset"],
    help="Perform ResolVI workflow on a joint dataset object",
)
@click.argument("dataset_name", type=str)
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Path to the AnnData directory")
@click.option("--model_save_path", type=click.Path(), default=None, help="Path to save the ResolVI model")
@click.option("--trained", type=bool, default=False, help="Whether the model is already trained (default: False)")
@click.option("--image_store", type=click.Path(exists=True), default=None, help="Path to save the plot")
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

    from spida.P.main import resolvi_dataset as func
    func(
        dataset_name=dataset_name,
        anndata_store=anndata_store,
        model_save_path=model_save_path,
        trained=trained,
        image_path=image_store,
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
@click.option("--anndata_store", type=click.Path(exists=True), default=None, help="Path to the anndata file")
@click.option("--image_store", type=click.Path(), default=None, help="Path to save the plot (default: None)")
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
    from spida.P.main import plot_dataset as func
    func(
        dataset_name=dataset_name,
        anndata_store=anndata_store,
        image_path=image_store,
        show=show,
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
