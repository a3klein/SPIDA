"""
click CLI interface for SPIDA segmentation pipeline using argparse.
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

logger = logging.getLogger(__name__)
# warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


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


@click.command(cls=RichCommand, help="Ingest spatial data for a given EXP_NAME and REG_NAME.")
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("type", "--type", default="merscope", type=str, help="Type of the data to ingest (default: merscope)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("source", "--source", default="machine", type=str, help="Source of the data (default: machine)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to store the zarr files (default: None, uses environment variable)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context
def ingest_region(
    ctx,
    exp_name: str,
    reg_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    root_path : str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Ingest spatial data for a specific EXP_NAME and REG_NAME.
    """
    from .io.main import ingest_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        type=type,
        prefix_name=prefix_name,
        source=source,
        plot=plot,
        root_path=root_path,
        zarr_store=zarr_store,
        image_store=image_store,
        **kwargs,
    )

@click.command(cls=RichCommand, help="Ingest spatial data for all regions in a given EXPERIMENT.")
@click.argument("exp_name", type=str)
@click.option("type", "--type", default="merscope", type=str, help="Type of the data to ingest (default: merscope)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("source", "--source", default="machine", type=str, help="Source of the data (default: machine)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to store the zarr files (default: None, uses environment variable)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context
def ingest_all(
    ctx,
    exp_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    root_path : str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Ingest spatial data for all regions in a specific EXPERIMENT.
    """
    from .io.main import ingest_all as func
    func(
        exp_name=exp_name,
        type=type,
        prefix_name=prefix_name,
        source=source,
        plot=plot,
        root_path=root_path,
        zarr_store=zarr_store,
        image_store=image_store,
        **kwargs,
    )

@click.command(cls=RichCommand, help="Load segmentation data for a given EXP_NAME and REG_NAME.")
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("seg_dir", type=click.Path(exists=True, dir_okay=True))
@click.option("type", "--type", default="vpt", type=str, help="Type of the segmentation data (default: vpt)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to store the zarr files (default: None, uses environment variable)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context 
def load_segmentation_region(
    ctx,
    exp_name: str,
    reg_name: str,
    seg_dir: str,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **load_kwargs,
):
    """
    Load segmentation data for a specific EXP_NAME and REG_NAME.
    """
    # kwargs
    extra_args = ctx.args
    load_kwargs = parse_click_kwargs(extra_args)
    for key, value in load_kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from .io.main import load_segmentation_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        seg_dir=seg_dir,
        type=type,
        prefix_name=prefix_name,
        plot=plot,
        zarr_store=zarr_store,
        image_store=image_store,
        **load_kwargs,
    )

@click.command(cls=RichCommand, help="Load segmentation data for all regions in a given EXPERIMENT.")
@click.argument("exp_name", type=str)
@click.argument("seg_dir", type=click.Path(exists=True, dir_okay=True))
@click.option("type", "--type", default="vpt", type=str, help="Type of the segmentation data (default: vpt)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to store the zarr files (default: None, uses environment variable)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context 
def load_segmentation_all(
    ctx,
    exp_name: str,
    seg_dir: str,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **load_kwargs,
):
    """
    Load segmentation data for all regions in a specific EXPERIMENT.
    """
    # kwargs
    extra_args = ctx.args
    load_kwargs = parse_click_kwargs(extra_args)
    for key, value in load_kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from .io.main import load_segmentation_all as func
    func(
        exp_name=exp_name,
        seg_dir=seg_dir,
        type=type,
        prefix_name=prefix_name,
        plot=plot,
        zarr_store=zarr_store,
        image_store=image_store,
        **load_kwargs,
    )

@click.command(
    name="load-decon-images",
    cls=RichCommand,
    help="Load deconvolution images for a given EXP_NAME and REG_NAME.",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.argument("image_dir", type=click.Path(exists=True, dir_okay=True))
@click.option("image_name", "--image_name", default="decon_image", type=str, help="Name of the image key in the spatialdata object (default: decon_image)")
@click.option("suffix", "--suffix", default=".decon.tif", type=str, help="Suffix of the deconvolution image files (default: .decon.tif)")
@click.option("z_layer", "--z_layer", default=3, type=str, help="Z-layer to read from the deconvolution images (default: 3)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to store the zarr files (default: None, uses environment variable)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context
def load_deconvolution_region(
    ctx,
    exp_name: str,
    reg_name: str,
    image_dir: str,
    image_name: str = "decon_image",
    suffix: str = ".decon.tif",
    z_layer: int | str = 3,
    plot: bool = False,
    root_path: str | Path | None = None,
    zarr_store: str | Path | None = None,
    image_store: str | Path | None = None,
    **load_kwargs,
):
    """
    Load deconvolution images into spatialdata objects.
    """

    extra_args = ctx.args
    load_kwargs = parse_click_kwargs(extra_args)
    for key, value in load_kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from .io.main import load_deconvolution_region as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        image_dir=image_dir,
        image_name=image_name,
        suffix=suffix,
        z_layer=z_layer,
        plot=plot,
        root_path=root_path,
        zarr_store=zarr_store,
        image_store=image_store,
        **load_kwargs,
    )


@cli.command(
    name="test-vpt", 
    cls=RichCommand,
    help="Test VPT binary",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.pass_context
def test_vpt(ctx):
    """
    Test the VPT binary by running a simple command.
    """
    extra_args = ctx.args
    click.echo(extra_args)
    kwargs = parse_click_kwargs(extra_args)
    click.echo(kwargs)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    from .segmentation.vpt import _add_vpt_binary as func
    func()


@cli.command(
    name="run-segmentation-region",
    aliases=["run", "run_segmentation_region"],
    cls=RichCommand,
    help="['run', 'run_segmentation_region'], Run segmentation of kind TYPE for a given EXP_NAME and REG_NAME.",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("seg_type", type=str)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("input_dir", "--input_dir", default=None, type=click.Path(), help="Directory containing the input data (default: None, uses environment variable)")
@click.option("output_dir", "--output_dir", default=None, type=click.Path(), help="Directory to save the output data (default: None, uses environment variable)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Path to store the segmentation results (default: None, uses environment variable)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Path to the VPT binary (default: None, uses environment variable)")
@click.option("rust_bin_path", "--rust_bin_path", default=None, type=click.Path(), help="Path to the Rust binary (default: None, uses environment variable)")
@click.pass_context
def run_segmentation(
    ctx,
    seg_type: str,
    exp_name: str,
    reg_name: str,
    input_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    Run an implemented segmentation algorithm on a given region
    """
    # logger.info(ctx.obj.get('_raw_extra_args'))
    # logger.info(ctx.obj.get('extra_args'))
    # kwargs = ctx.obj['extra_kwargs']
    # for key, value in kwargs.items():
    #     logger.info(f"Parsing extra argument {key}={value}")

    extra_args = ctx.args
    click.echo(extra_args)
    kwargs = parse_click_kwargs(extra_args)
    click.echo(kwargs)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    

    from .segmentation.main import run_segmentation as func
    func(
        type=seg_type,
        exp_name=exp_name,
        reg_name=reg_name,
        input_dir=input_dir,
        output_dir=output_dir,
        root_path=root_path,
        segmentation_store=segmentation_store,
        vpt_bin_path=vpt_bin_path,
        kwargs=kwargs,
    )

@cli.command(
    name="segment-experiment",
    aliases=["experiment", "segment_experiment"],
    cls=RichCommand,
    help="['experiment', 'segment_experiment'], Run segmentation of kind TYPE for all regions in a given EXPERIMENT.",
)
@click.argument("type", type=str)
@click.argument("exp_name", type=str)
@click.option("input_dir", "--input_dir", default=None, type=click.Path(), help="Directory containing the input data (default: None, uses environment variable)")
@click.option("output_dir", "--output_dir", default=None, type=click.Path(), help="Directory to save the output data (default: None, uses environment variable)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Path to store the segmentation results (default: None, uses environment variable)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Path to the Zarr storage (default: None, uses environment variable)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Path to the VPT binary (default: None, uses environment variable)")
@click.option("rust_bin_path", "--rust_bin_path", default=None, type=click.Path(), help="Path to the Rust binary (default: None, uses environment variable)")
@click.option("kwargs", "--kwargs", default="{}", type=str, help="Additional keyword arguments as a JSON string (default: {})")
@click.pass_context
def segment_experiment(
    ctx,
    type: str,
    exp_name: str,
    input_dir: Path = None,
    output_dir: Path = None,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    zarr_store: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs
):
    """
    Run segmentation for all regions in an experiment.
    """
    from .segmentation.main import segment_experiment as func
    func(
        type=type,
        exp_name=exp_name,
        input_dir=input_dir,
        output_dir=output_dir,
        root_path=root_path,
        segmentation_store=segmentation_store,
        zarr_store=zarr_store,
        vpt_bin_path=vpt_bin_path,
        rust_bin_path=rust_bin_path,
        **kwargs,
    )

    
@cli.command(
    name="align-proseg",
    aliases=["align", "align_proseg"],
    cls=RichCommand,
    help="['align', 'align_proseg'], Align Proseg transcripts to seed transcripts for a given EXP_NAME and REG_NAME.",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("seed_prefix_name", "--seed_prefix_name", default="default", type=str, help="Prefix for the seed transcripts (default: default)")
@click.option("prefix_name", "--prefix_name", default="proseg", type=str, help="Prefix for the Proseg transcripts (default: proseg)")
@click.option("out_prefix_name", "--out_prefix_name", default="proseg_aligned", type=str, help="Prefix for the output aligned transcripts (default: proseg_aligned)")   
@click.option("input_dir", "--input_dir", default=None, type=click.Path(), help="Directory containing the input data (default: None, uses environment variable)")
@click.option("seg_dir", "--seg_dir", default=None, type=click.Path(), help="Directory containing the segmentation data (default: None, uses environment variable)")
@click.option("x", "--x", default="x", type=str, help="Column name for the x coordinates (default: x)")
@click.option("y", "--y", default="y", type=str, help="Column name for the y coordinates (default: y)")
@click.option("z", "--z", default="global_z", type=str, help="Column name for the z coordinates (default: global_z)")
@click.option("cell_column", "--cell_column", default="cell_id", type=str, help="Column name for the cell IDs (default: cell_id)")
@click.option("barcode_column", "--barcode_column", default="barcode_id", type=str, help="Column name for the barcode IDs (default: barcode_id)")
@click.option("gene_column", "--gene_column", default="gene", type=str, help="Column name for the gene names (default: gene)")
@click.option("fov_column", "--fov_column", default="fov", type=str, help="Column name for the field of view (default: fov)")
@click.option("cell_missing", "--cell_missing", default="-1", type=int, help="Value indicating missing cell IDs (default: -1)")
@click.option("min_jaccard", "--min_jaccard", default=0.4, type=float, help="Minimum Jaccard index for alignment (default: 0.4)")
@click.option("min_prob", "--min_prob", default=0.5, type=float, help="Minimum probability for alignment (default: 0.5)")
@click.option("filter_blank", "--filter_blank", is_flag=True, default=False, help="Whether to filter out blank cells (default: False)")
@click.option("cell_metadata_fname", "--cell_metadata_fname", default="merged_cell_metadata.csv", type=str, help="Filename for the cell metadata (default: merged_cell_metadata.csv)")
@click.option("cell_by_gene_fname", "--cell_by_gene_fname", default="merged_cell_by_gene.csv", type=str, help="Filename for the cell by gene data (default: merged_cell_by_gene.csv)")
@click.option("detected_transcripts_fname", "--detected_transcripts_fname", default="merged_transcript_metadata.csv", type=str, help="Filename for the detected transcripts data (default: merged_transcript_metadata.csv)")
@click.option("cell_polygons_fname", "--cell_polygons_fname", default="merged_cell_polygons.geojson", type=str, help="Filename for the cell polygons data (default: merged_cell_polygons.geojson)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Path to the VPT binary (default: None, uses environment variable)")
@click.option("kwargs", "--kwargs", default="{}", type=str, help="Additional keyword arguments as a JSON string (default: {})")
@click.pass_context
def align_proseg(
    ctx,
    exp_name: str,
    reg_name: str,
    seed_prefix_name: str = "default",
    prefix_name: str = "proseg",
    out_prefix_name: str = "proseg_aligned",
    input_dir: str | Path = None,
    seg_dir: str | Path = None,
    x: str = "x",
    y: str = "y",
    z: str = "global_z",
    cell_column: str = "cell_id",
    barcode_column: str = "barcode_id",
    gene_column: str = "gene",
    fov_column: str = "fov",
    cell_missing: int = "-1",
    min_jaccard: float = 0.4,
    min_prob: float = 0.5,
    filter_blank: bool = False,
    cell_metadata_fname: str = "merged_cell_metadata.csv",
    cell_by_gene_fname: str = "merged_cell_by_gene.csv",
    detected_transcripts_fname: str = "merged_transcript_metadata.csv",
    cell_polygons_fname: str = "merged_cell_polygons.geojson",
    vpt_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    Align Proseg transcripts to seed transcripts.
    """
    from .segmentation.main import align_proseg as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        seed_prefix_name=seed_prefix_name,
        prefix_name=prefix_name,
        out_prefix_name=out_prefix_name,
        input_dir=input_dir,
        seg_dir=seg_dir,
        x=x,
        y=y,
        z=z,
        cell_column=cell_column,
        barcode_column=barcode_column,
        gene_column=gene_column,
        fov_column=fov_column,
        cell_missing=cell_missing,
        min_jaccard=min_jaccard,
        min_prob=min_prob,
        filter_blank=filter_blank,
        cell_metadata_fname=cell_metadata_fname,
        cell_by_gene_fname=cell_by_gene_fname,
        detected_transcripts_fname=detected_transcripts_fname,
        cell_polygons_fname=cell_polygons_fname,
        vpt_bin_path=vpt_bin_path,
        **kwargs,
    )

@cli.command(
    name="filt-to-ids",
    aliases=["filt_to_ids", "filter_to_ids", "filter-to-ids"],
    cls=RichCommand,
    help="['filt_to_ids', 'filter_to_ids', 'filter-to-ids'] - Filter metadata and geometry to only include cells with IDs in the provided list.",
)
@click.argument("meta_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("geom_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("tz_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("cbg_path", type=click.Path(exists=True, dir_okay=False))
@click.pass_context
def filt_to_ids(
    ctx,
    meta_path: str | Path,
    geom_path: str | Path,
    tz_path: str | Path,
    cbg_path: str | Path,
):
    """
    Filter metadata and geometry to only include cells with IDs in the provided list.
    Additionally rename the index of the merged_cell_by_gene
    """
    from .segmentation.proseg import filt_to_ids as func
    func(
        meta_path=meta_path,
        geom_path=geom_path,
        tz_path=tz_path,
        cbg_path=cbg_path,
    )


@cli.command(
    name="align-segmentation",
    aliases=["align_segmentation", "align-segmentation"],
    cls=RichCommand,
    help="['align_segmentation', 'align-segmentation'] - Align segmentation results between different methods.",
)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("prefix1", "--prefix1", default="cellpose_nuc", type=str, help="Prefix for the first segmentation method (default: cellpose_nuc)")
@click.option("prefix2", "--prefix2", default="cellpose_cell", type=str, help="Prefix for the second segmentation method (default: cellpose_cell)")
@click.option("output_dir", "--output_dir", default=None, type=click.Path(), help="Directory to save the aligned segmentation results (default: None, uses environment variable)")
@click.option("geometry_mode", "--geometry_mode", default="larger", type=str, help="Geometry mode for alignment (default: larger)")
@click.option("cell_id", "--cell_id", default="EntityID", type=str, help="Column name for the cell IDs (default: EntityID for vpt)")
@click.option("min_intersection_area", "--min_intersection_area", default=0.0, type=float, help="Minimum intersection area for alignment (default: 0.0)")
@click.option("coordinate_system", "--coordinate_system", default="global", type=str, help="Coordinate system for alignment (default: global)")
@click.option("out_dir_name", "--out_dir_name", default="align", type=str, help="The name of the output directory to store the segmentation under (default: align)")
@click.option("zarr_store", "--zarr_store", default=None, type=str, help="Path to the zarr store (default: None)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=str, help="Path to the segmentation store (default: None)")
@click.option("root_path", "--root_path", default=None, type=str, help="Root path for the data (default: None)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=str, help="Path to the VPT binary (default: None)")
@click.pass_context
def align_geometries(
    ctx,
    exp_name,
    reg_name,
    prefix1 : str = "cellpose_nuc",
    prefix2 : str = "cellpose_cell",
    output_dir : str | None = None,
    geometry_mode : str = "larger", 
    cell_id : str = "EntityID",
    min_intersection_area : float = 0.0,
    coordinate_system : str = "global",
    out_dir_name : str = "align",
    zarr_store : str | None = None,
    segmentation_store: str | None = None,
    root_path : str | None = None,
    vpt_bin_path : str | None = None
): 
    """
    Align segmentation results between different methods.
    """
    from .segmentation.main import align_geometries as func
    func(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir=output_dir,
        geometry_mode=geometry_mode,
        cell_id=cell_id,
        min_intersection_area=min_intersection_area,
        coordinate_system=coordinate_system,
        out_dir_name=out_dir_name,
        zarr_store=zarr_store,
        segmentation_store=segmentation_store,
        root_path=root_path,
        vpt_bin_path=vpt_bin_path,
    )


@cli.command(
    name="decon-image",
    aliases=["decon_image"],
    cls=RichCommand,
    help="['decon_image', 'decon-image'] - Deconvolve an image  using Deconwolf.",
)
@click.option("image_path", "--image_path", "-i", default=None, type=click.Path(exists=True, dir_okay=True), help="Path to the image file to be deconvolved (default: None, throws an error if not provided)")
@click.option("data_org_path", "--data_org_path", default="{input}/dataorganization.csv", type=click.Path(), help="Path to the data organization CSV file (default: {input}/dataorganization.csv)")
@click.option("channels", "--channels", default="DAPI", type=str, help="Channels to use for deconvolution, comma-separated if multiple (default: DAPI)")
@click.option("tile_size", "--tile_size", default=2960, type=int, help="Tile size for processing (default: 2960)")
@click.option("overlap", "--overlap", default=100, type=int, help="Overlap size between tiles (default: 100)")
@click.option("output_dir", "--output_dir", '-o', default="tiles_output", type=click.Path(), help="Directory to save the output tiles (default: tiles_output)")
@click.option("filter_args", "--filter_args", default=None, type=(str,int), multiple=True, help="Additional filter arguments as a JSON string (default: None)")
@click.option("filter", "--filter", default="deconwolf", type=str, help="Filter to use for deconvolution (default: deconwolf)")
@click.option("gpu", "--gpu", is_flag=True, default=True, help="Whether to use GPU for deconvolution (default: True)")
@click.option("z_step", "--z_step", default=1.5, type=float, help="Z-step size for 3D images (default: 1.5)")
@click.option("continue_stalled", "--continue_stalled", is_flag=True, default=False, help="Whether to continue stalled processes (default: False)")
@click.option("thr_tiles", "--thr_tiles/--no_thr_tiles", is_flag=True, default=True, help="Whether to threshold tiles after deconvolution (default: True)")
@click.option("plot_thr", "--plot_thr/--no_plot_thr", is_flag=True, default=False, help="Whether to plot the thresholding intermediary (default: False)")
@click.option("match_pre", "--match_pre", is_flag=True, default=False, help="Whether to match pre-existing tiles histogram distribution (default: False)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Path to store the images (default: None, uses environment variable)")
@click.pass_context
def decon_image(
    ctx,
    image_path: str | Path | None = None,
    data_org_path: str | Path = "{input}/dataorganization.csv",
    channels: str | list[str] = "DAPI",
    tile_size: int | tuple = 2960,
    overlap: int | tuple = 100,
    output_dir: str | Path = "tiles_output",
    filter_args: dict = None,
    filter: str = "deconwolf",
    gpu: bool = True,
    z_step: float = 1.5,
    continue_stalled: bool = False,
    thr_tiles: bool = True,
    plot_thr: bool = False,
    match_pre: bool = False,
    image_store: str | Path | None = None,
    **kwargs,
):
    """
    Deconvolve an image using Deconwolf.
    """
    if image_path is None:
        raise ValueError("image_path must be provided.")
        
    filter_args = dict(filter_args) if filter_args is not None else None
    logger.info(f"Running decon_image")
    for key, value in ctx.params.items():
        logger.info(f"{key}: {value} ({type(value)})")


    from .decon_script import decon_image as func
    func(
        image_path=image_path,
        data_org_path=data_org_path,
        channels=channels,
        tile_size=tile_size,
        overlap=overlap,
        output_dir=output_dir,
        filter_args=filter_args,
        filter=filter,
        gpu=gpu,
        z_step=z_step,
        continue_stalled=continue_stalled,
        thr_tiles=thr_tiles,
        plot_thr=plot_thr,
        match_pre=match_pre,
        image_store=image_store,
        **kwargs,
    ) 

# cli.add_command(test_command)
cli.add_command(ingest_region)
cli.add_command(ingest_all)
cli.add_command(load_segmentation_region)
cli.add_command(load_segmentation_all)
cli.add_command(load_deconvolution_region)
cli.add_command(run_segmentation)
cli.add_command(segment_experiment)
cli.add_command(align_proseg)
cli.add_command(filt_to_ids)
cli.add_command(align_geometries)
cli.add_command(decon_image)

if __name__ == "__main__":
    
    # TODO: Functionalize this configuration setup? 
    # - Technically once SPIDA is fully built this will move into the main entry point of the package.
    # - and this won't be done once in each module! 

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