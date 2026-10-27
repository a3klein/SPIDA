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
from spida.io.options import dataset_options, dataset_overrides, dataset_from_options, run_config
from spida.utilities.script_utils import parse_click_kwargs, parse_literal, JSONParam

logger = logging.getLogger(__name__)

# warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


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


@click.command(cls=RichCommand, help="Ingest one region's raw data into a new zarr store (--exp/--reg, or the path options).")
@dataset_options
@click.option("type", "--type", default="merscope", type=str, help="Type of the data to ingest (default: merscope)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("source", "--source", default="machine", type=str, help="Source of the data (default: machine)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace this region's existing zarr store, including everything added after ingest (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@dataset_overrides
@click.pass_context
def ingest_region(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    overwrite: bool = False,
    root_path : str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Ingest one region's raw data into a new zarr store.
    """
    cfg = run_config(ctx, root_path=root_path, zarr_store=zarr_store, image_store=image_store)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .load.main import ingest_region as func
    func(dataset, type=type, prefix_name=prefix_name, source=source, plot=plot, overwrite=overwrite, **kwargs)

@click.command(cls=RichCommand, help="Ingest every region of --exp found under PROCESSED_ROOT_PATH.")
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("type", "--type", default="merscope", type=str, help="Type of the data to ingest (default: merscope)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Prefix for the keys in the spatialdata object (default: default)")
@click.option("source", "--source", default="machine", type=str, help="Source of the data (default: machine)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace existing zarr stores; without it, regions that already have one are skipped (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.pass_context
def ingest_all(
    ctx,
    exp_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    overwrite: bool = False,
    root_path : str | Path | None = None,
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **kwargs,
):
    """
    Ingest every region of an experiment found under the raw root.
    """
    cfg = run_config(ctx, root_path=root_path, zarr_store=zarr_store, image_store=image_store)
    from .load.main import ingest_all as func
    func(cfg, exp_name, type=type, prefix_name=prefix_name, source=source, plot=plot, overwrite=overwrite, **kwargs)

@click.command(cls=RichCommand, help="Load one segmentation into a region's zarr store (--exp/--reg, or the path options).")
@dataset_options
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("type", "--type", default="vpt", type=str, help="Type of the segmentation data (default: vpt)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Segmentation name, used for its output folder and its names in the zarr store (default: default)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("use_transcript_qc", "--transcript-qc", is_flag=True, default=False, help="Apply QC filtering at load-time using transcript QC regions. (default: False)")
@click.option("qc_shapes_key", "--qc_shapes_key", default="qc_mask", type=str,help="Shapes key holding the QC pass-region mask to apply (default: qc_mask; e.g. transcript_qc_shapes for the legacy hex mask)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace this segmentation's points, shapes and table if already loaded (default: False)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@dataset_overrides
@click.pass_context
def load_segmentation_region(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    segmentation_store: str | Path | None = None,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    use_transcript_qc: bool = False,
    overwrite: bool = False,
    qc_shapes_key: str = "qc_mask",
    zarr_store : str | Path | None = None,
    image_store : str | Path | None = None,
    **load_kwargs,
):
    """
    Load one segmentation into a region's zarr store.
    """
    # kwargs
    extra_args = ctx.args
    load_kwargs = parse_click_kwargs(extra_args)
    if use_transcript_qc:
        load_kwargs.setdefault("qc_shapes_key", qc_shapes_key)
        logger.info(
            "Transcript QC filtering enabled for load-segmentation-region with qc_shapes_key=%s",
            load_kwargs["qc_shapes_key"],
        )

    cfg = run_config(ctx, segmentation_store=segmentation_store, zarr_store=zarr_store, image_store=image_store)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .load.main import load_segmentation_region as func
    func(dataset, type=type, prefix_name=prefix_name, plot=plot, overwrite=overwrite, **load_kwargs)

@click.command(cls=RichCommand, help="Load one segmentation for every region of --exp that has it.")
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("type", "--type", default="vpt", type=str, help="Type of the segmentation data (default: vpt)")
@click.option("prefix_name", "--prefix_name", default="default", type=str, help="Segmentation name, used for its output folder and its names in the zarr store (default: default)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("use_transcript_qc", "--transcript-qc", is_flag=True, default=False, help="Apply QC filtering at load-time using transcript QC regions. (default: False)")
@click.option("qc_shapes_key", "--qc_shapes_key", default="qc_mask", type=str,help="Shapes key holding the QC pass-region mask to apply (default: qc_mask; e.g. transcript_qc_shapes for the legacy hex mask)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace segmentations already loaded; without it, those regions are skipped (default: False)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.pass_context
def load_segmentation_all(
    ctx,
    exp_name: str,
    segmentation_store: str | Path | None = None,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    use_transcript_qc: bool = False,
    overwrite: bool = False,
    qc_shapes_key: str = "qc_mask",
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
    if use_transcript_qc:
        load_kwargs.setdefault("qc_shapes_key", qc_shapes_key)
        logger.info(
            "Transcript QC filtering enabled for load-segmentation-all with qc_shapes_key=%s",
            load_kwargs["qc_shapes_key"],
        )
    cfg = run_config(ctx, segmentation_store=segmentation_store, zarr_store=zarr_store, image_store=image_store)
    from .load.main import load_segmentation_all as func
    func(cfg, exp_name, type=type, prefix_name=prefix_name, plot=plot, overwrite=overwrite, **load_kwargs)

@click.command(
    name="load-decon-images",
    cls=RichCommand,
    help="Load deconvolution images into a region's zarr store (--exp/--reg, or the path options).",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@dataset_options
@click.option("image_dir", "--image_dir", default=None, type=click.Path(exists=True, file_okay=False), help="Directory holding the deconvolved images (default: the region's raw images dir)")
@click.option("image_name", "--image_name", default="decon_image", type=str, help="Name of the image key in the spatialdata object (default: decon_image)")
@click.option("suffix", "--suffix", default=".decon.tif", type=str, help="Suffix of the deconvolution image files (default: .decon.tif)")
@click.option("z_layer", "--z_layer", default=3, type=str, help="Z-layer to read from the deconvolution images (default: 3)")
@click.option("plot", "--plot", is_flag=True, default=False, help="Whether to generate plots (default: False)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace the image if already loaded (default: False)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@dataset_overrides
@click.pass_context
def load_deconvolution_region(
    ctx,
    exp, reg, zarr_path, raw_path, seg_path, figures_path,
    image_dir: str | None = None,
    image_name: str = "decon_image",
    suffix: str = ".decon.tif",
    z_layer: int | str = 3,
    plot: bool = False,
    overwrite: bool = False,
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

    cfg = run_config(ctx, root_path=root_path, zarr_store=zarr_store, image_store=image_store)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .load.main import load_deconvolution_region as func
    func(dataset, image_dir=image_dir, image_name=image_name, suffix=suffix,
         z_layer=z_layer, plot=plot, overwrite=overwrite, **load_kwargs)


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
    name="segment-region",
    aliases=["segment"],
    cls=RichCommand,
    help="Run the segmentation backend for METHOD on EXP_NAME/REG_NAME (raw output only). "
         "Run in the backend's env (cellpose -> cellpose; proseg -> preprocessing; mesmer -> deepcell).",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
@click.argument("method", type=str)
@dataset_options
@click.option("--list-params", "--list_params", "list_params", is_flag=True, default=False, help="Print every run parameter METHOD accepts (with types, defaults and deprecated aliases) and exit. No dataset needed.")
@click.option("--prefix_name", "prefix_name", default=None, type=str, help="Segmentation name / output-dir leaf (default: method; e.g. cellpose_nuc / cellpose_cell). Must match process/load.")
@click.option("--root_path", "root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("--segmentation_store", "segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("--rust_bin_path", "rust_bin_path", default=None, type=click.Path(), help="Directory holding the proseg binary (default: RUST_BIN_PATH in the config)")
@click.option("--tissue_mask_path", "tissue_mask_path", default=None, type=click.Path(), help="Tissue mask restricting where segmentation runs (default: auto-discover <images_dir>/tissue_mask.tif, else ungated).")
@click.option("--overwrite", "overwrite", is_flag=True, default=False, help="Replace this segmentation's existing raw output (default: False)")
@dataset_overrides
@click.pass_context
def segment_region(ctx, method, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                   list_params, prefix_name, root_path, segmentation_store, rust_bin_path,
                   tissue_mask_path, overwrite):
    # Method-specific run parameters are NOT exposed as Click options.
    # `--list-params` makes them discoverable instead, from the method's own config class.
    if list_params:
        from .segmentation.backends import describe_method_params
        click.echo(describe_method_params(method))
        return

    kwargs = parse_click_kwargs(ctx.args)
    cfg = run_config(ctx, root_path=root_path, segmentation_store=segmentation_store, rust_bin_path=rust_bin_path)
    if method == "cellpose" and "model_dir" not in kwargs and cfg.cellpose_model_path is not None:
        kwargs["model_dir"] = str(cfg.cellpose_model_path)  # CELLPOSE_MODEL_PATH; --model_dir wins
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .segmentation.main import segment_region as func
    func(method, dataset, prefix_name=prefix_name,
         rust_bin_path=cfg.rust_bin_path, tissue_mask_path=tissue_mask_path, overwrite=overwrite, **kwargs)


@cli.command(
    name="process-segmentation-region",
    aliases=["process-segmentation"],
    cls=RichCommand,
    help="Normalize + post-process a segmenter's raw output into the segmentation schema for "
         "METHOD on EXP_NAME/REG_NAME. Run in the preprocessing env.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
@click.argument("method", type=str)
@dataset_options
@click.option("--version", "version", default=None, type=str, help="Method version (proseg: 2 or 3)")
@click.option("--prefix_name", "prefix_name", default=None, type=str, help="Segmentation name / output-dir leaf (default: method). Must match segment-region.")
@click.option("--backend", "backend", default="native", type=click.Choice(["native", "vpt"]), help="native (default) or vpt fallback")
@click.option("--root_path", "root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("--segmentation_store", "segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("--micron_per_z", "micron_per_z", default=1.5, type=float, help="Micron thickness per z-plane")
@click.option("--n_z_planes", "n_z_planes", default=7, type=int, help="Number of z-planes for 2D->3D replication")
@click.option("--n_jobs", "n_jobs", default=7, type=int, help="sum-signals parallel workers (default 7; IO-bound step, more is slower + costs more SUs)")
@click.option("--vpt_bin_path", "vpt_bin_path", default=None, type=click.Path(), help="Directory holding the VPT binary, --backend vpt only (default: VPT_BIN_PATH in the config)")
@click.option("--overwrite", "overwrite", is_flag=True, default=False, help="Replace existing segmentation-schema files (default: False)")
@dataset_overrides
@click.pass_context
def process_segmentation_region(ctx, method, exp, reg, zarr_path, raw_path, seg_path, figures_path,
                                version, prefix_name, backend, root_path, segmentation_store,
                                micron_per_z, n_z_planes, n_jobs, vpt_bin_path, overwrite):
    kwargs = parse_click_kwargs(ctx.args)
    cfg = run_config(ctx, root_path=root_path, segmentation_store=segmentation_store, vpt_bin_path=vpt_bin_path)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .segmentation.main import process_segmentation_region as func
    func(method, dataset, version=version, prefix_name=prefix_name, backend=backend,
         micron_per_z=micron_per_z, n_z_planes=n_z_planes, n_jobs=n_jobs,
         vpt_bin_path=cfg.vpt_bin_path, overwrite=overwrite, **kwargs)


@cli.command(
    name="process-custom-segmentation",
    aliases=["process-custom"],
    cls=RichCommand,
    help="Process a USER-PROVIDED segmentation (polygons + optional transcripts / stain images) into the segmentation schema.",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
@click.argument("boundaries_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("output_dir", type=click.Path())
@click.option("--boundaries_space", "boundaries_space", default="micron", type=click.Choice(["micron", "pixel"]), help="Coordinate space of the input polygons (default micron)")
@click.option("--micron_to_mosaic_path", "micron_to_mosaic_path", default=None, type=click.Path(), help="Path to micron_to_mosaic_pixel_transform.csv (required for pixel space and/or images)")
@click.option("--z_spacing", "z_spacing", default=1.5, type=float, help="Micron thickness per z-plane (default 1.5)")
@click.option("--n_z_planes", "n_z_planes", default=1, type=int, help="Number of z-planes; 1 => 2D cylinder expansion (default 1)")
@click.option("--cell_id_col", "cell_id_col", default=None, type=str, help="Caller's cell-id column in the boundaries (required for 3D / metadata merge)")
@click.option("--boundary_z_col", "boundary_z_col", default=None, type=str, help="Integer z-plane column in the boundaries for 3D input")
@click.option("--transcripts_path", "transcripts_path", default=None, type=click.Path(), help="Detected transcripts (.csv/.parquet); enables partition")
@click.option("--transcript_z_col", "transcript_z_col", default="global_z", type=str, help="Transcript z column (default global_z)")
@click.option("--transcript_z_in_microns", "transcript_z_in_microns", is_flag=True, default=False, help="Transcript z is micron ZLevel (convert to integer plane)")
@click.option("--images_dir", "images_dir", default=None, type=click.Path(), help="Dir of mosaic_{stain}_z{N}.tif images; enables sum-signals")
@click.option("--segmentation_z_index", "segmentation_z_index", default=None, type=int, help="For 2D, the single ZIndex plane for sum-signals (default middle image plane)")
@click.option("--metadata_path", "metadata_path", default=None, type=click.Path(), help="Optional user cell-metadata CSV merged onto derived metadata")
@click.option("--metadata_cell_id_col", "metadata_cell_id_col", default=None, type=str, help="Cell-id column in metadata_path if it differs from cell_id_col")
@click.option("--n_jobs", "n_jobs", default=7, type=int, help="sum-signals parallel workers (default 7)")
@click.option("--overwrite", "overwrite", is_flag=True, default=False, help="Replace existing segmentation-schema files in OUTPUT_DIR (default: False)")
@click.pass_context
def process_custom_segmentation(ctx, boundaries_path, output_dir, boundaries_space,
                                micron_to_mosaic_path, z_spacing, n_z_planes, cell_id_col,
                                boundary_z_col, transcripts_path, transcript_z_col,
                                transcript_z_in_microns, images_dir, segmentation_z_index,
                                metadata_path, metadata_cell_id_col, n_jobs, overwrite):
    kwargs = parse_click_kwargs(ctx.args)
    from .segmentation.main import process_custom_segmentation as func
    func(boundaries_path, output_dir, boundaries_space=boundaries_space,
         micron_to_mosaic_path=micron_to_mosaic_path, z_spacing=z_spacing,
         n_z_planes=n_z_planes, cell_id_col=cell_id_col, boundary_z_col=boundary_z_col,
         transcripts_path=transcripts_path, transcript_z_col=transcript_z_col,
         transcript_z_in_microns=transcript_z_in_microns, images_dir=images_dir,
         segmentation_z_index=segmentation_z_index, metadata_path=metadata_path,
         metadata_cell_id_col=metadata_cell_id_col, n_jobs=n_jobs, overwrite=overwrite, **kwargs)




@cli.command(
    name="align-segmentation",
    aliases=["align_segmentation"],
    cls=RichCommand,
    help="Align two segmentations of one region by overlap, and build the aligned segmentation's output files in SEGMENTATION_OUT_PATH/<exp>/<reg>/OUT_DIR_NAME.",
)
@dataset_options
@click.option("prefix1", "--prefix1", default="cellpose_nuc", type=str, help="First segmentation (default: cellpose_nuc)")
@click.option("prefix2", "--prefix2", default="cellpose_cell", type=str, help="Second segmentation (default: cellpose_cell)")
@click.option("geometry_mode", "--geometry_mode", default="larger", type=click.Choice(["larger", "prefix1", "prefix2", "intersection"]), help="Polygon kept per matched pair (default: larger)")
@click.option("cell_id", "--cell_id", default="EntityID", type=str, help="Column name for the cell IDs (default: EntityID)")
@click.option("min_intersection_area", "--min_intersection_area", default=0.0, type=float, help="Minimum intersection area for a match, pixels (default: 0.0)")
@click.option("coordinate_system", "--coordinate_system", default="global", type=str, help="Coordinate system to align in (default: global = mosaic pixels)")
@click.option("out_dir_name", "--out_dir_name", default="align", type=str, help="Name of the aligned segmentation, i.e. its output folder (default: align)")
@click.option("micron_per_z", "--micron_per_z", default=1.5, type=float, help="Micron thickness per z-plane (default: 1.5)")
@click.option("n_jobs", "--n_jobs", default=7, type=int, help="sum-signals workers (default: 7)")
@click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Root of the segmentation outputs (default: SEGMENTATION_OUT_PATH in the config)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root of the raw MERSCOPE output (default: PROCESSED_ROOT_PATH in the config)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, help="Replace an existing aligned segmentation (default: False)")
@dataset_overrides
@click.pass_context
def align_geometries(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, prefix1, prefix2, geometry_mode,
                     cell_id, min_intersection_area, coordinate_system, out_dir_name, micron_per_z, n_jobs,
                     zarr_store, segmentation_store, root_path, overwrite):
    """Align two segmentations of one region and build the aligned segmentation's output files."""
    cfg = run_config(ctx, zarr_store=zarr_store, segmentation_store=segmentation_store, root_path=root_path)
    dataset = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path)
    from .segmentation.main import align_geometries as func
    func(dataset, prefix1=prefix1, prefix2=prefix2, geometry_mode=geometry_mode, cell_id=cell_id,
         min_intersection_area=min_intersection_area, coordinate_system=coordinate_system,
         out_dir_name=out_dir_name, micron_per_z=micron_per_z, n_jobs=n_jobs, overwrite=overwrite)


@cli.command(
    name="transform-image",
    aliases=["transform_image", "decon-image", "decon_image"],
    cls=RichCommand,
    help="['transform_image', 'decon_image'] - Apply a per-tile transform (deconwolf or dark_sectioning) to a large image.",
)
@click.option("image_path", "--image_path", "-i", default=None, type=click.Path(exists=True, dir_okay=True), help="Path to the image file to be deconvolved (default: None, throws an error if not provided)")
@click.option("data_org_path", "--data_org_path", default="{input}/dataorganization.csv", type=click.Path(), help="Path to the data organization CSV file (default: {input}/dataorganization.csv)")
@click.option("channels", "--channels", default="DAPI", type=str, help="Channels to use for deconvolution, comma-separated if multiple (default: DAPI)")
@click.option("tile_size", "--tile_size", default=2960, type=int, help="Tile size for processing (default: 2960)")
@click.option("overlap", "--overlap", default=100, type=int, help="Overlap size between tiles (default: 100)")
@click.option("output_dir", "--tile_dir", "--output_dir", '-o', default="tiles_output", type=click.Path(), help="WORKING directory for per-tile intermediates: the split tiles, each tile's .decon.tif, and each tile's .decon.3d.tif projection. All of it is deleted after stitching unless --no_clear_tiles. Effectively scratch space. --output_dir/-o are kept as aliases and will be depreceated soon; --tile_dir is the accurate name. The stitched full-mosaic planes go to --stitched_dir.")
@click.option("stitched_dir", "--stitched_dir", default=None, type=click.Path(), help="Directory for the stitched *_z{n}.decon.tif outputs. Default writes within the directory of image_path and overwrites any existing deconvolution.")
@click.option("transform_args", "--transform_args", "--filter_args", default=None, type=(str,str), multiple=True, help="Arguments for the chosen transform, as repeated `--transform_args KEY VALUE` pairs. Values are passed through literal-eval (`iter 50`, `psigma 1.5`, `gpu True`); anything unparseable is kept as a string.")
@click.option("transform", "--transform", "--filter", default="deconwolf", type=click.Choice(["deconwolf", "dark_sectioning"]), help="Per-tile transform to apply (default: deconwolf)")
@click.option("gpu", "--gpu", is_flag=True, default=True, help="Whether to use GPU for deconvolution (default: True)")
@click.option("z_step", "--z_step", default=1.5, type=float, help="Z-step size for 3D images (default: 1.5)")
@click.option("z_stack_expand", "--z_stack_expand", is_flag=True, default=True, help="Whether to expand z-stack in 3D case (default: True)")
@click.option("z_stack_expand_slices", "--z_stack_expand_slices", default=3, type=int, help="Number of slices to expand z-stack (default: 3)")
@click.option("continue_stalled", "--continue_stalled", is_flag=True, default=False, help="Whether to continue stalled processes (default: False)")
@click.option("thr_tiles", "--thr_tiles/--no_thr_tiles", is_flag=True, default=True, help="Whether to threshold tiles after deconvolution (default: True)")
@click.option("plot_thr", "--plot_thr/--no_plot_thr", is_flag=True, default=False, help="Whether to plot the thresholding intermediary (default: False)")
@click.option("select_tiles_method", "--select-tiles", "--select_tiles", "--image-filter", "--image_filter", default=None, type=click.Choice(["otsu", "tissue_mask", "none"]), help="Tile-selection strategy; overrides --thr_tiles (default: None -> otsu if thr_tiles else none)")
@click.option("mask_path", "--mask_path", default=None, type=click.Path(), help="Path to tissue_mask.tif for --image-filter tissue_mask (default: <image_path>/tissue_mask.tif)")
@click.option("match_pre", "--match_pre", is_flag=True, default=False, help="Whether to match pre-existing tiles histogram distribution (default: False)")
@click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)")
@click.option("deconwolf_config_file", "--deconwolf_config_file", default=None, type=click.Path(), help="deconwolf .config.ini: binary, PSFs, [Optics] (default: DECONWOLF_CONFIG in the config)")
@click.pass_context
def transform_image(
    ctx,
    image_path: str | Path | None = None,
    data_org_path: str | Path = "{input}/dataorganization.csv",
    channels: str | list[str] = "DAPI",
    tile_size: int | tuple = 2960,
    overlap: int | tuple = 100,
    output_dir: str | Path = "tiles_output",
    stitched_dir: str | Path | None = None,
    transform_args: dict = None,
    transform: str = "deconwolf",
    gpu: bool = True,
    z_step: float = 1.5,
    z_stack_expand: bool = True,
    z_stack_expand_slices: int = 3,
    continue_stalled: bool = False,
    thr_tiles: bool = True,
    plot_thr: bool = False,
    select_tiles_method: str | None = None,
    mask_path: str | Path | None = None,
    match_pre: bool = False,
    image_store: str | Path | None = None,
    deconwolf_config_file: str | Path | None = None,
    **kwargs,
):
    """
    Deconvolve an image using Deconwolf.
    """
    if image_path is None:
        raise ValueError("image_path must be provided.")
        
    transform_args = (
        {k: parse_literal(v) for k, v in dict(transform_args).items()}
        if transform_args is not None
        else None
    )
    logger.info(f"Running transform_image")
    for key, value in ctx.params.items():
        logger.info(f"{key}: {value} ({type(value)})")


    _cfg = run_config(ctx, image_store=image_store, deconwolf_config_file=deconwolf_config_file)
    from spida.S.image_prep.orchestrator import transform_image as func
    func(
        image_path=image_path,
        data_org_path=data_org_path,
        channels=channels,
        tile_size=tile_size,
        overlap=overlap,
        output_dir=output_dir,
        transform_args=transform_args,
        transform=transform,
        gpu=gpu,
        z_step=z_step,
        z_stack_expand=z_stack_expand,
        z_stack_expand_slices=z_stack_expand_slices,
        continue_stalled=continue_stalled,
        thr_tiles=thr_tiles,
        plot_thr=plot_thr,
        select_tiles_method=select_tiles_method,
        mask_path=mask_path,
        stitched_dir=stitched_dir,
        match_pre=match_pre,
        image_store=_cfg.image_store,
        deconwolf_config=_cfg.deconwolf_config_file,
        **kwargs,
    )

# cli.add_command(test_command)
cli.add_command(ingest_region)
cli.add_command(ingest_all)
cli.add_command(load_segmentation_region)
cli.add_command(load_segmentation_all)
cli.add_command(load_deconvolution_region)
cli.add_command(align_geometries)
cli.add_command(transform_image)

if __name__ == "__main__":
    
    # TODO: Functionalize this configuration setup? 
    # - Technically once SPIDA is fully built this will move into the main entry point of the package.
    # - and this won't be done once in each module! 

    # SPIDA_LOG_LEVEL raises verbosity for one run without a code change. There is
    # otherwise no way to see DEBUG output from a SLURM job, where the process is
    # launched by a batch script rather than a person. Unset -> INFO, i.e. the
    # previous behaviour exactly.
    _level = getattr(logging, os.environ.get("SPIDA_LOG_LEVEL", "INFO").upper(),
                     logging.INFO)

    # Handlers must accept the chosen level for records to get through
    env = configure_logging_for_runtime(level=_level)
    # Root stays at WARNING to suppress other packages regardless
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)

    # The spida package itself follows the requested level
    spida_logger = logging.getLogger('spida')
    spida_logger.setLevel(_level)
    logger.setLevel(_level)
    logger.info(f"Logging configured for environment: {env}")
            
    cli()