### Integration / Annotations main file
# author: Amit Klein, a3klein@ucsd.edu

from pathlib import Path
import warnings
import logging
from rich.logging import RichHandler
import click
from rich_click import RichCommand  # type: ignore
from spida.utilities.script_utils import parse_click_kwargs, JSONParam, parse_list
from spida.settings import configure_logging_for_runtime
from spida.config import ConfigDefaultGroup
from spida.io.options import dataset_options, dataset_overrides, table_options, run_config, dataset_from_options

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

def setup_logging(**kwargs):
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s - %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
        handlers=[RichHandler(rich_tracebacks=True, show_time=False, markup=True)],
    )

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


def _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path, **roots):
    return dataset_from_options(run_config(ctx, **roots), exp, reg, zarr_path, raw_path, seg_path, figures_path,
                                anndata_path, annotations_path)


def _options(opts):
    def deco(f):
        for opt in reversed(opts):
            f = opt(f)
        return f
    return deco


def _extra_kwargs(ctx) -> dict:
    kwargs = parse_click_kwargs(ctx.args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")
    return kwargs


_TABLE_INPUT = [
    click.option('--suffix', type=click.STRING, default='_filt', help='Table suffix, e.g. _filt for the filtered cells. Defaults to "_filt".'),
    click.option('--adata_file', type=click.Path(exists=True, path_type=str), default=None, help='Annotate this AnnData file (h5ad or zarr) instead of the table of the region; labels are written back to it.'),
    click.option('--refresh_table', is_flag=True, default=False, help="Re-read the table from the zarr store even if its AnnData file already exists (drops annotations in that file)."),
]
_STORE_ROOTS = [
    click.option("zarr_store", "--zarr_store", default=None, type=click.Path(), help="Root of the zarr stores (default: ZARR_STORAGE_PATH in the config)"),
    click.option("anndata_store", "--anndata_store", default=None, type=click.Path(), help="Root of the AnnData (.h5ad) files (default: ANNDATA_STORE_PATH in the config)"),
    click.option("annotation_store", "--annotation_store", default=None, type=click.Path(), help="Root of the annotation intermediates (default: ANNOTATION_STORE_PATH in the config)"),
    click.option("image_store", "--image_store", default=None, type=click.Path(), help="Root of the figures (default: IMAGE_STORE_PATH in the config)"),
]


#### BACKUP ####
_SYNC = [
    click.option('--suffix', type=click.STRING, default='_filt', help='Table suffix. Defaults to "_filt".'),
    click.option('--reference', type=click.STRING, default=None, help='Only sync this reference (mmc_<reference>_* / allcools_<reference>_*). Default: every annotation column.'),
    click.option('--columns', type=click.STRING, multiple=True, help='Sync exactly these obs columns instead (repeatable).'),
    click.option('--overwrite', is_flag=True, default=False, help='Replace columns that already exist in the zarr table.'),
]


@cli.command(name="sync-annotations", cls=RichCommand)
@click.argument('prefix_name', type=click.STRING)
@_options(_SYNC)
@click.option('--adata_file', type=click.Path(exists=True, path_type=str), default=None, help="AnnData to sync from (default: the region's file in the AnnData store).")
@dataset_options
@_options(_STORE_ROOTS)
@dataset_overrides
@table_options
@click.pass_context
def sync_annotations_region(
    ctx,
    prefix_name: str,
    suffix: str = "_filt",
    reference: str | None = None,
    columns: tuple = (),
    overwrite: bool = False,
    adata_file: str | None = None,
    exp: str | None = None,
    reg: str | None = None,
    zarr_path: str | None = None,
    raw_path: str | None = None,
    seg_path: str | None = None,
    figures_path: str | None = None,
    anndata_path: str | None = None,
    annotations_path: str | None = None,
    **roots,
):
    """
    Copy the annotation columns of segmentation PREFIX_NAME's AnnData file into the region's zarr table (obs only).
    Runs in any env with anndata (the mmc env included); does not need spatialdata.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
    """
    from spida.I._io import sync_annotations

    d = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path, **roots)
    sync_annotations(d, prefix_name, suffix, adata_file=adata_file, reference=reference,
                     columns=columns or None, overwrite=overwrite)


@cli.command(name='sync-annotations-experiment', cls=RichCommand)
@click.argument('prefix_name', type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@_options(_SYNC)
@_options(_STORE_ROOTS)
@click.pass_context
def sync_annotations_experiment(
    ctx,
    prefix_name: str,
    exp_name: str,
    suffix: str = "_filt",
    reference: str | None = None,
    columns: tuple = (),
    overwrite: bool = False,
    **roots,
):
    """
    sync-annotations for every region of --exp; regions whose columns are already in the zarr table
    are skipped unless --overwrite.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        exp_name: Name of the experiment.\n
    """
    from spida.io.layout import TreeLayout
    from spida.I._io import sync_annotations

    layout = TreeLayout(run_config(ctx, **roots))
    for reg in layout.units(exp_name):
        try:
            sync_annotations(layout.dataset(exp_name, reg), prefix_name, suffix, reference=reference,
                             columns=columns or None, overwrite=overwrite)
        except FileExistsError as e:
            logger.warning("Skipping %s: %s", reg, e)


#### ALLCools INTEGRATION ####
_ALLCOOLS_OPTIONS = [
    click.option('--reference_name', type=click.STRING, default=None, help='Name of this reference in column names and paths. Defaults to the REF_PATH file stem.'),
    click.option("--gene_rename_dict", type=JSONParam(), default=None, help='JSON string (representing a dict), path to a JSON file, or a dict when calling programmatically. Defaults to None.'),
    click.option("--top_deg_genes", type=int, default=-1, help='Number of top differentially expressed genes to use from each cluster, set to -1 to turn off. Defaults to -1.'),
    click.option("--max_cells_per_cluster", type=int, default=3000, help='Maximum number of cells per cluster when downsampling the reference AnnData object. Defaults to 3000.'),
    click.option("--min_cells_per_cluster", type=int, default=0, help='Minimum number of cells per cluster when downsampling the reference AnnData object. Defaults to 0.'),
    click.option("--label_transfer_k", type=int, default=100, help='Number of neighbors to use for label transfer. Defaults to 100.'),
    click.option("--downsample_ref", is_flag=True, default=False, help='Whether to downsample the reference AnnData object. Defaults to False.'),
    click.option("--min_ref_cells", type=int, default=50000, help='Minimum number of cells to downsample to in the reference AnnData object. Defaults to 50000.'),
    click.option("--save_integrator/--no-save_integrator", default=True, help='Whether to save the integrator object. Defaults to True.'),
    click.option("--save_adata_comb/--no-save_adata_comb", default=True, help='Whether to save the combined AnnData object. Defaults to True.'),
    click.option("--rna_cell_type_column", type=click.STRING, default="supercluster_name", help='Column name in the AnnData object for RNA cell types. Defaults to "supercluster_name".'),
    click.option("--qry_cluster_column", type=click.STRING, default="leiden", help='Column name in the AnnData object for query clusters. Defaults to "leiden".'),
    click.option("--run_joint_embeddings", is_flag=True, default=False, help='Whether to run joint embeddings on the integrated data. Defaults to False.'),
    click.option("--joint_embedding_leiden_res", type=float, default=1.5, help='Resolution for clustering the joint embeddings when running joint embeddings. Defaults to 1.5.'),
    click.option("--run_clust_label_transfer/--no-run_clust_label_transfer", default=True, help='Whether to run cluster-to-cluster label transfer on the integrated data. Defaults to True.'),
    click.option("--confusion_matrix_cluster_min_value", type=float, default=0.25, help='Minimum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.25.'),
    click.option("--confusion_matrix_cluster_max_value", type=float, default=0.9, help='Maximum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.9.'),
    click.option("--confusion_matrix_cluster_resolution", type=float, default=1.5, help='Resolution for clustering the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 1.5.'),
    click.option("--qry_only_cluster_threshold", type=int, default=50, help='Min Reference cell number for defining query-only clusters during cluster-to-cluster label transfer. Defaults to 50.'),
    click.option("--ref_only_cluster_threshold", type=int, default=20, help='Min Query cell number for defining reference-only clusters during cluster-to-cluster label transfer. Defaults to 20.'),
    click.option("--plot", is_flag=True, default=False, help='Whether to plot the results of the ALLCools integration. Defaults to False.'),
]


def _allcools_one(
    d,
    prefix_name: str,
    ref_path: str,
    suffix: str,
    adata_file: str | None,
    refresh_table: bool,
    reference_name: str | None = None,
    plot: bool = False,
    **opts,
):
    """
    Integrate one region (or one AnnData file) and optionally plot it.

    Parameters:
        d: The region's locations (Dataset).
        reference_name: Name of the reference in column names and paths (default: REF_PATH file stem).
        **opts: ALLCools options, passed through to run_allcools_seurat.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.allcools import integrate_region, plot_integration

    logger.info("RUNNING ALLCOOLS INTEGRATION, %s, PREFIX %s" % (d.reg or adata_file, prefix_name))
    qry, _, query_path, output_dir, _ = integrate_region(
        d, prefix_name, ref_path, suffix=suffix, reference=reference_name,
        adata_file=adata_file, refresh_table=refresh_table, **opts)
    logger.info("DONE WITH ALLCOOLS INTEGRATION")
    if plot:
        reference = reference_name or Path(ref_path).stem
        plot_integration(
            query_path, output_dir / "adata_comb_rna_merfish.h5ad",
            d.figures_dir(prefix_name) / f"allcools_integration_{reference}.pdf", reference,
            plot_joint_embeddings=opts.get("run_joint_embeddings", False),
            plot_c2c_transfer=opts.get("run_clust_label_transfer", True),
            ref_cell_type_column=opts.get("rna_cell_type_column", "supercluster_name"),
        )


_CTX = dict(ignore_unknown_options=True, allow_extra_args=True)


@cli.command(name='allcools-integration-region', cls=RichCommand, context_settings=_CTX)
@click.argument('prefix_name', type=click.STRING)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@_options(_TABLE_INPUT)
@dataset_options
@_options(_STORE_ROOTS)
@_options(_ALLCOOLS_OPTIONS)
@dataset_overrides
@table_options
@click.pass_context
def allcools_integration_region(
    ctx,
    prefix_name: str,
    ref_path: str,
    suffix: str = "_filt",
    adata_file: str | None = None,
    refresh_table: bool = False,
    exp: str | None = None,
    reg: str | None = None,
    zarr_path: str | None = None,
    raw_path: str | None = None,
    seg_path: str | None = None,
    figures_path: str | None = None,
    anndata_path: str | None = None,
    annotations_path: str | None = None,
    zarr_store: str | None = None,
    anndata_store: str | None = None,
    annotation_store: str | None = None,
    image_store: str | None = None,
    **opts,
):
    """
    Run ALLCools integration of segmentation PREFIX_NAME for one region (--exp/--reg), or of --adata_file,
    against the reference RNA AnnData REF_PATH.

    Labels are added to the region's AnnData file (in ANNDATA_STORE_PATH) as allcools_<reference>_* columns.
    The integrator and combined AnnData go to ANNOTATION_STORE_PATH/<exp>/<reg>/<table>/allcools/<reference>.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        ref_path: Path to the reference RNA AnnData object.\n
        **opts: Additional keyword arguments for ALLCools integration (e.g. --deg_type=cef).\n
    """
    d = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path,
                 zarr_store=zarr_store, anndata_store=anndata_store, annotation_store=annotation_store, image_store=image_store)
    _allcools_one(d, prefix_name, ref_path, suffix, adata_file, refresh_table, **opts, **_extra_kwargs(ctx))


@cli.command(name='allcools-integration-experiment', cls=RichCommand, context_settings=_CTX)
@click.argument('prefix_name', type=click.STRING)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option('--suffix', type=click.STRING, default='_filt', help='Table suffix, e.g. _filt for the filtered cells. Defaults to "_filt".')
@click.option('--refresh_table', is_flag=True, default=False, help="Re-read each table from the zarr store even if its AnnData file already exists.")
@_options(_STORE_ROOTS)
@_options(_ALLCOOLS_OPTIONS)
@click.pass_context
def allcools_integration_experiment(
    ctx,
    prefix_name: str,
    ref_path: str,
    exp_name: str,
    suffix: str = "_filt",
    refresh_table: bool = False,
    zarr_store: str | None = None,
    anndata_store: str | None = None,
    annotation_store: str | None = None,
    image_store: str | None = None,
    **opts,
):
    """
    Run ALLCools integration of segmentation PREFIX_NAME for every region of an entire experiment --exp.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        ref_path: Path to the reference RNA AnnData object.\n
        exp_name: Name of the experiment.\n
        **opts: Additional keyword arguments for ALLCools integration.\n
    """
    from spida.io.layout import TreeLayout

    layout = TreeLayout(run_config(ctx, zarr_store=zarr_store, anndata_store=anndata_store,
                                   annotation_store=annotation_store, image_store=image_store))
    extra = _extra_kwargs(ctx)
    for reg in layout.units(exp_name):
        _allcools_one(layout.dataset(exp_name, reg), prefix_name, ref_path, suffix, None, refresh_table, **opts, **extra)


@cli.command(name='plot-allcools-integration', cls=RichCommand)
@click.argument('prefix_name', type=click.STRING)
@click.argument('reference_name', type=click.STRING)
@click.option('--suffix', type=click.STRING, default='_filt', help='Table suffix. Defaults to "_filt".')
@click.option('--adata_file', type=click.Path(exists=True, path_type=str), default=None, help="Annotated query (default: the region's file in the AnnData store).")
@click.option('--output_path', type=click.Path(path_type=str), default=None, help='Output pdf (default: IMAGE_STORE_PATH/<exp>/<reg>/<prefix>/allcools_integration_<reference>.pdf).')
@dataset_options
@_options(_STORE_ROOTS)
@click.option("--plot_joint_embeddings", is_flag=True, default=False, help='Whether to plot the joint embeddings. Defaults to False.')
@click.option("--plot_c2c_transfer", is_flag=True, default=False, help='Whether to plot the cluster-to-cluster label transfer results. Defaults to False.')
@click.option("--ref_cell_type_column", type=click.STRING, default="supercluster_name", help='Column name in the AnnData object for RNA cell types. Defaults to "supercluster_name".')
@dataset_overrides
@table_options
@click.pass_context
def plot_allcools_integration(
    ctx,
    prefix_name: str,
    reference_name: str,
    suffix: str = "_filt",
    adata_file: str | None = None,
    output_path: str | None = None,
    exp: str | None = None,
    reg: str | None = None,
    zarr_path: str | None = None,
    raw_path: str | None = None,
    seg_path: str | None = None,
    figures_path: str | None = None,
    anndata_path: str | None = None,
    annotations_path: str | None = None,
    zarr_store: str | None = None,
    anndata_store: str | None = None,
    annotation_store: str | None = None,
    image_store: str | None = None,
    plot_joint_embeddings: bool = False,
    plot_c2c_transfer: bool = False,
    ref_cell_type_column: str = "supercluster_name",
):
    """
    Plot the results of ALLCools integration REFERENCE_NAME for segmentation PREFIX_NAME of one region.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        reference_name: The reference name used by the integration (default there: the reference file stem).\n
    """
    from spida.io.layout import table_key
    from spida.I.allcools import plot_integration

    d = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path,
                 zarr_store=zarr_store, anndata_store=anndata_store, annotation_store=annotation_store, image_store=image_store)
    key = table_key(prefix_name, suffix)
    plot_integration(
        adata_file or d.table_h5ad(key),
        d.annotation_dir(key, "allcools", reference_name) / "adata_comb_rna_merfish.h5ad",
        output_path or d.figures_dir(prefix_name) / f"allcools_integration_{reference_name}.pdf",
        reference_name,
        plot_joint_embeddings=plot_joint_embeddings,
        plot_c2c_transfer=plot_c2c_transfer,
        ref_cell_type_column=ref_cell_type_column,
    )


#### MAPMYCELLS ####
_MMC_RUN = [
    click.option('--mmc_store_path', type=click.Path(file_okay=False, path_type=str, dir_okay=True), default=None, help="Store of MapMyCells references (default: MMC_DIR in the config)"),
    click.option('--n_cpu', type=int, default=1, help='Number of CPUs to use. Defaults to 1.'),
    click.option('--bootstrap_factor', type=float, default=0.8, help='Bootstrap factor for annotation. Defaults to 0.8.'),
    click.option('--bootstrap_iterations', type=int, default=100, help='Number of bootstrap iterations. Defaults to 100.'),
    click.option('--rng_seed', type=int, default=13, help='Random number generator seed. Defaults to 13.'),
    click.option("--plot", is_flag=True, default=False, help='Whether to plot the annotation results. Defaults to False.'),
    click.option("--plot_name", type=click.STRING, default="mmc_annotation_plots", help='Plot file name, suffixed with the identifier. Defaults to "mmc_annotation_plots".'),
    click.option("--palette_path", type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False), default=None, help='Path to a custom palette file for plotting. Defaults to None.'),
]


def _mmc_config(
    cfg,
    mmc_store_path: str | None,
    n_cpu: int,
    bootstrap_factor: float,
    bootstrap_iterations: int,
    rng_seed: int,
):
    """MMCConfig with paths from the SPIDA config (``mmc_store_path`` overrides MMC_DIR)."""
    from spida.I.mmc import MMCConfig
    return MMCConfig.from_spida(cfg, mmc_store_path=mmc_store_path, n_cpu=n_cpu, bootstrap_factor=bootstrap_factor,
                                bootstrap_iterations=bootstrap_iterations, rng_seed=rng_seed)


@cli.command(name='mmc-setup', cls=RichCommand, context_settings=_CTX)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@click.argument('identifier', type=click.STRING)
@click.option('--hierarchy_list', type=parse_list, default="Class,Subclass,Group", help='Comma-separated list of hierarchy levels (e.g., "level1,level2,level3")')
@click.option('--gene_names_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False), default=None, help='Path to gene_names.txt file (default: searched in GENE_PANEL_PATH by identifier).')
@click.option('--gene_name_mapping_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False), default=None, help='Path to gene_name_mapping.json file (default: searched in GENE_PANEL_PATH by identifier).')
@click.option('--mmc_store_path', type=click.Path(file_okay=False, path_type=str, dir_okay=True), default=None, help="Store of MapMyCells references (default: MMC_DIR in the config)")
@click.option('--gene_panel_path', type=click.Path(file_okay=False, path_type=str, dir_okay=True), default=None, help="Directory of gene panel files (default: GENE_PANEL_PATH in the config)")
@click.option('--ref_norm', type=click.STRING, default='log2CPM', help='Normalization method for the reference AnnData.X. Defaults to "log2CPM".')
@click.option('--n_cpu', type=int, default=8, help='Number of CPUs to use. Defaults to 8.')
@click.option('--n_valid', type=int, default=10, help='Number of valid markers. Defaults to 10.')
@click.option('--n_per_utility', type=int, default=10, help='Number of markers per utility. Defaults to 10.')
@click.pass_context
def mmc_setup(
    ctx,
    ref_path: str,
    identifier: str,
    hierarchy_list: list = ["Class", "Subclass", "Group"],
    gene_names_path: str | None = None,
    gene_name_mapping_path: str | None = None,
    mmc_store_path: str | None = None,
    gene_panel_path: str | None = None,
    ref_norm: str = "log2CPM",
    n_cpu: int = 8,
    n_valid: int = 10,
    n_per_utility: int = 10,
):
    """
    Setup function for MapMyCells: precompute the reference stats and markers for REF_PATH under IDENTIFIER.

    Parameters:\n
        ref_path: Path to the reference single-cell RNA-seq AnnData file.\n
        identifier: Unique identifier for this reference (e.g. 'mouse_motor_cortex').\n
        hierarchy_list: Comma-separated hierarchy levels (e.g. "level1,level2,level3").\n
        gene_names_path: gene_names.txt (auto-detected in GENE_PANEL_PATH if not provided).\n
        gene_name_mapping_path: gene_name_mapping.json (auto-detected in GENE_PANEL_PATH if not provided).\n
        n_valid: Number of valid markers.\n
        n_per_utility: Number of markers per utility.\n
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import MMCConfig, setup_annotation_pipeline

    cfg = run_config(ctx, mmc_store_path=mmc_store_path, gene_panel_path=gene_panel_path)
    config = MMCConfig.from_spida(cfg, ref_norm=ref_norm, n_cpu=n_cpu)
    logger.info(f"Hierarchy list: {hierarchy_list}")
    result = setup_annotation_pipeline(
        ref_path, identifier, hierarchy_list, config,
        gene_names_path=gene_names_path, gene_name_mapping_path=gene_name_mapping_path,
        n_valid=n_valid, n_per_utility=n_per_utility,
    )
    click.echo(f"Setup complete. Results: {result}")


@cli.command(name="mmc-annotation-region", cls=RichCommand, context_settings=_CTX)
@click.argument('prefix_name', type=click.STRING)
@click.argument('identifier', type=click.STRING)
@_options(_TABLE_INPUT)
@dataset_options
@_options(_STORE_ROOTS)
@_options(_MMC_RUN)
@click.option("--filter_annotations/--no-filter_annotations", default=False, help='Whether to filter annotations based on scores and counts. Defaults to False.')
@dataset_overrides
@table_options
@click.pass_context
def mmc_annotation_region(
    ctx,
    prefix_name: str,
    identifier: str,
    suffix: str = "_filt",
    adata_file: str | None = None,
    refresh_table: bool = False,
    exp: str | None = None,
    reg: str | None = None,
    zarr_path: str | None = None,
    raw_path: str | None = None,
    seg_path: str | None = None,
    figures_path: str | None = None,
    anndata_path: str | None = None,
    annotations_path: str | None = None,
    zarr_store: str | None = None,
    anndata_store: str | None = None,
    annotation_store: str | None = None,
    image_store: str | None = None,
    mmc_store_path: str | None = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    plot: bool = False,
    plot_name: str = "mmc_annotation_plots",
    palette_path: str | None = None,
    filter_annotations: bool = False,
):
    """
    Run MapMyCells annotation of segmentation PREFIX_NAME for one region (--exp/--reg), or of --adata_file,
    against reference IDENTIFIER.

    Labels are added to the region's AnnData file (in ANNDATA_STORE_PATH) as mmc_<identifier>_<level> columns.
    MapMyCells results go to ANNOTATION_STORE_PATH/<exp>/<reg>/<table>/mmc/<identifier>.
    Extra --lN_thr / --lN_count_thr options go to the annotation filter.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        identifier: Reference identifier (e.g. 'mouse_motor_cortex'), as given to mmc-setup.\n
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import annotate_region

    cfg = run_config(ctx, zarr_store=zarr_store, anndata_store=anndata_store, annotation_store=annotation_store,
                     image_store=image_store)
    d = dataset_from_options(cfg, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path)
    annotate_region(
        d, prefix_name, suffix, identifier,
        _mmc_config(cfg, mmc_store_path, n_cpu, bootstrap_factor, bootstrap_iterations, rng_seed),
        adata_file=adata_file, refresh_table=refresh_table, filter_annot=filter_annotations,
        plot=plot, palette_path=palette_path, plot_name=plot_name, **_extra_kwargs(ctx),
    )
    click.echo(f"Annotation complete for {d.reg or adata_file}")


@cli.command(name="mmc-annotation-experiment", cls=RichCommand, context_settings=_CTX)
@click.argument('prefix_name', type=click.STRING)
@click.argument('identifier', type=click.STRING)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option('--suffix', type=click.STRING, default='_filt', help='Table suffix, e.g. _filt for the filtered cells. Defaults to "_filt".')
@click.option('--refresh_table', is_flag=True, default=False, help="Re-read each table from the zarr store even if its AnnData file already exists.")
@_options(_STORE_ROOTS)
@_options(_MMC_RUN)
@click.option("--filter_annotations/--no-filter_annotations", default=True, help='Whether to filter annotations based on scores and counts. Defaults to True.')
@click.pass_context
def mmc_annotation_experiment(
    ctx,
    prefix_name: str,
    identifier: str,
    exp_name: str,
    suffix: str = "_filt",
    refresh_table: bool = False,
    zarr_store: str | None = None,
    anndata_store: str | None = None,
    annotation_store: str | None = None,
    image_store: str | None = None,
    mmc_store_path: str | None = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    plot: bool = False,
    plot_name: str = "mmc_annotation_plots",
    palette_path: str | None = None,
    filter_annotations: bool = True,
):
    """
    Run MapMyCells annotation of segmentation PREFIX_NAME for every region of an entire experiment --exp.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        identifier: Reference identifier (e.g. 'mouse_motor_cortex'), as given to mmc-setup.\n
        exp_name: Name of the experiment.\n
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import annotate_experiment

    cfg = run_config(ctx, zarr_store=zarr_store, anndata_store=anndata_store, annotation_store=annotation_store,
                     image_store=image_store)
    annotate_experiment(
        cfg, exp_name, prefix_name, suffix, identifier,
        _mmc_config(cfg, mmc_store_path, n_cpu, bootstrap_factor, bootstrap_iterations, rng_seed),
        refresh_table=refresh_table, filter_annot=filter_annotations,
        plot=plot, palette_path=palette_path, plot_name=plot_name, **_extra_kwargs(ctx),
    )
    click.echo(f"Annotation complete for experiment {exp_name}")


#### MOSCOT INTEGRATION ####
@cli.command(name="moscot-integration", cls=RichCommand, context_settings=_CTX)
def moscot_integration():
    raise NotImplementedError("MOSCOT integration is not implemented yet.")


if __name__ == "__main__":

    # Configure root logger with INFO level handlers (to allow INFO messages through)
    # but set the root logger level to WARNING (to suppress other modules)
    env = configure_logging_for_runtime(
        level=logging.INFO,  # Handlers need to accept INFO level
    )

    # Set root logger level to WARNING to suppress other modules
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)

    # Set the entire spida package to INFO level as well
    spida_logger = logging.getLogger('spida')
    spida_logger.setLevel(logging.INFO)

    # Configure the spida.I module logger (parent for all files in this module)
    module_logger = logging.getLogger('spida.I')
    module_logger.setLevel(logging.INFO)

    logger.setLevel(logging.INFO) # Set the level for the current logger

    # You can also set specific third-party modules to different levels if needed
    # For example, to suppress verbose output from specific libraries:
    # logging.getLogger('matplotlib').setLevel(logging.WARNING)
    # logging.getLogger('anndata').setLevel(logging.WARNING)
    # logging.getLogger('scanpy').setLevel(logging.WARNING)

    logger.info(f"Logging configured for environment: {env}")

    cli()
