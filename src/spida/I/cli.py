### Integration / Annotations main file
# author: Amit Klein, a3klein@ucsd.edu

import os
from dotenv import load_dotenv  # type: ignore
from pathlib import Path
import warnings
import logging
from rich.logging import RichHandler
# import fire  # type: ignore
import click
from rich_click import RichCommand  # type: ignore
from spida.utilities.script_utils import parse_click_kwargs, JSONParam, parse_list
from spida.settings import configure_logging_for_runtime
from spida.config import ConfigDefaultGroup

load_dotenv()
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
    default=".env",
    type=click.Path(exists=True),
    help="Path to the configuration file (default: .env)",
)
@click.pass_context
def cli(ctx, config):
    ctx.ensure_object(dict)
    pass

@cli.command(name="backup-adata", cls=RichCommand)
@click.argument('exp_name', type=click.STRING)
@click.argument('reg_name', type=click.STRING)
@click.argument('prefix_name', type=click.STRING)
@click.option('--adata_path',
              type=click.Path(file_okay=True, path_type=str, dir_okay=True),
              default=None,
              help='Path to the store of AnnData objects. Defaults to None.'
            )
def backup_adata_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    adata_path: str = None
):
    """
    Backup the AnnData objects for experiment EXP_NAME, region REG_NAME with prefix PREFIX_NAME.
    Can serve as input for integration methods incompatible with the spatialdata package (i.e. mapmycells).
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I._utils import backup_adata_region as func
    func(exp_name, reg_name, prefix_name, adata_path)


@cli.command(name='backup-adata-experiment', cls=RichCommand)
@click.argument('exp_name', type=str)
@click.argument('prefix_name', type=str)
@click.option(
    '--adata_path',
    type=click.Path(file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
def backup_adata_experiment(exp_name: str, prefix_name: str, adata_path: str = None):
    """
    Backup function for AnnData objects for an entire experiment EXP_NAME with prefix PREFIX_NAME.

    Arguments:\n
    exp_name (str): Name of the experiment.\n
    prefix_name (str): Prefix for the keys in the spatialdata object.\n
    """

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I._utils import backup_adata_experiment as func
    func(exp_name, prefix_name, adata_path)

#### ALLCools INTEGRATION ####
@cli.command(
    name='allcools-integration-region',
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('exp_name', type=click.STRING)
@click.argument('reg_name', type=click.STRING)
@click.argument('prefix_name', type=click.STRING)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@click.option('--suffix', type=click.STRING, default='_filt', help='Suffix for the keys in the spatialdata object. Defaults to "_filt".')
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.option("--gene_rename_dict",
              type=JSONParam(),
              default=None,
              help='JSON string (representing a dict), path to a JSON file, or a dict when calling programmatically. Defaults to None.'
)
@click.option(
    "--top_deg_genes",
    type=int,
    default=-1,
    help='Number of top differentially expressed genes to use from each cluster, set to -1 to turn off. Defaults to -1.')
@click.option("--max_cells_per_cluster", type=int, default=3000, help='Maximum number of cells per cluster when downsampling the reference AnnData object. Defaults to 3000.')
@click.option("--min_cells_per_cluster", type=int, default=0, help='Minimum number of cells per cluster when downsampling the reference AnnData object. Defaults to 0.')
@click.option("--label_transfer_k", type=int, default=100, help='Number of neighbors to use for label transfer. Defaults to 100.')
@click.option("--downsample_ref", type=bool, is_flag=True, default=False, help='Whether to downsample the reference AnnData object. Defaults to False.')
@click.option("--min_ref_cells", type=int, default=50000, help='Minimum number of cells to downsample to in the reference AnnData object. Defaults to 50000.')
@click.option("--save_integrator", type=bool, is_flag=True, default=True, help='Whether to save the integrator object. Defaults to True.')
@click.option("--save_adata_comb", type=bool, is_flag=True, default=True, help='Whether to save the combined AnnData object. Defaults to True.')
@click.option("--rna_cell_type_column", type=click.STRING, default="supercluster_name", help='Column name in the AnnData object for RNA cell types. Defaults to "supercluster_name".')
@click.option("--qry_cluster_column", type=click.STRING, default="leiden", help='Column name in the AnnData object for query clusters. Defaults to "leiden".')
@click.option("--run_joint_embeddings", type=bool, is_flag=True, default=False, help='Whether to run joint embeddings on the integrated data. Defaults to False.')
@click.option("--joint_embedding_leiden_res", type=float, default=1.5, help='Resolution for clustering the joint embeddings when running joint embeddings. Defaults to 1.5.')
@click.option("--run_clust_label_transfer", type=bool, is_flag=True, default=True, help='Whether to run cluster-to-cluster label transfer on the integrated data. Defaults to True.')
@click.option("--confusion_matrix_cluster_min_value", type=float, default=0.25, help='Minimum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.25.')
@click.option("--confusion_matrix_cluster_max_value", type=float, default=0.9, help='Maximum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.9.')
@click.option("--confusion_matrix_cluster_resolution", type=float, default=1.5, help='Resolution for clustering the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 1.5.')
@click.option("--qry_only_cluster_threshold", type=int, default=50, help='Min Reference cell number for defining query-only clusters during cluster-to-cluster label transfer. Defaults to 50.')
@click.option("--ref_only_cluster_threshold", type=int, default=20, help='Min Query cell number for defining reference-only clusters during cluster-to-cluster label transfer. Defaults to 20.')
@click.option("--plot", type=bool, is_flag=True, default=False, help='Whether to plot the results of the ALLCools integration. Defaults to False.')
@click.option("--image_path", type=click.Path(exists=False, file_okay=False, path_type=str, dir_okay=True), default=None, help='Path to the image file for plotting. Defaults to None.')
@click.option("--backup_to_spatialdata", type=bool, is_flag=True, default=True, help='Whether to backup the results to the spatialdata object. Defaults to True.')
@click.pass_context
def allcools_integration_region(
    ctx,
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    ref_path: str,
    suffix:str = "_filt",
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    gene_rename_dict: dict | Path | str = None,
    top_deg_genes: int = -1, 
    max_cells_per_cluster = 3000,
    min_cells_per_cluster = 0,
    label_transfer_k: int = 100,
    downsample_ref: bool = False,
    min_ref_cells: int = 50000,
    save_integrator:bool = True, 
    save_adata_comb: bool = True,
    rna_cell_type_column: str = "supercluster_name",
    qry_cluster_column:str = "leiden",
    run_joint_embeddings:bool = False,
    joint_embedding_leiden_res: float = 1.5,
    run_clust_label_transfer: bool = True,
    confusion_matrix_cluster_min_value: float = 0.25,
    confusion_matrix_cluster_max_value: float = 0.9,
    confusion_matrix_cluster_resolution: float = 1.5,
    qry_only_cluster_threshold: int = 50,
    ref_only_cluster_threshold: int = 20,
    plot:bool = False,
    image_path:str | Path | None = None,
    backup_to_spatialdata:bool = True, 
    **kwargs,
):
    """
    Run ALLCools integration on a given EXP_NAME, REG_NAME, and PREFIX_NAME

    Arguments:\n
    exp_name (str): Name of the experiment.\n
    reg_name (str): Name of the region.\n
    prefix_name (str): Prefix for the keys in the spatialdata object.\n
    ref_path (str): Path to the reference RNA AnnData object .\n
    **kwargs: Additional keyword arguments for ALLCools integration.\n
    """
    import json
    import anndata as ad # type: ignore
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.allcools import run_allcools_seurat
    
    logger.info(
        "RUNNING ALLCOOLS INTEGRATION, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )
    extra_args = ctx.args
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")

    # # Getting the adata object
    zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
    ad_path = f"{zarr_path}/tables/{prefix_name}_table{suffix}"
    logger.info(f"ad_path: {ad_path}")
    
    adata = ad.read_zarr(ad_path)
    logger.info(f"adata shape: {adata.shape}")

    # Handling the gene renaming: 
    if gene_rename_dict is not None:
        logger.info(f"Renaming genes using gene_rename_dict: {gene_rename_dict}")
        adata.var_names = [gene_rename_dict.get(name, name) for name in adata.var_names]

    ref_adata = ad.read_h5ad(ref_path, backed='r')
    logger.info(f"ref_adata shape: {ref_adata.shape}")
    if downsample_ref:
        logger.info("Downsampling the reference AnnData object")
        # Never less than the number in the reference, or the min_ref_cells params
        num_cells = min(ref_adata.shape[0], max(adata.shape[0], min_ref_cells))
        ref_adata = ref_adata[ref_adata.obs.sample(num_cells).index]
        logger.info(f"ref_adata shape after downsampling: {ref_adata.shape}")
    
    # Subsetting to shared genes before loading into memory:
    ref_adata = ref_adata[:, ref_adata.var_names[ref_adata.var_names.isin(adata.var_names)]]
    ref_adata = ref_adata.to_memory()

    logger.info(f"Missing Gene Names: {adata.var_names.difference(ref_adata.var_names)}")

    # Loading global parameters for storing anndata and annotations related files 
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

    logger.info("read data and calling the function")
    adata, _ = run_allcools_seurat(
        ref_adata=ref_adata, 
        qry_adata=adata, 
        anndata_store_path=anndata_store_path,
        annotations_store_path=annotations_store_path,
        rna_cell_type_column=rna_cell_type_column,
        qry_cluster_column=qry_cluster_column,
        top_deg_genes=top_deg_genes,
        max_cells_per_cluster=max_cells_per_cluster,
        min_cells_per_cluster=min_cells_per_cluster,
        label_transfer_k=label_transfer_k,
        run_joint_embeddings=run_joint_embeddings,
        joint_embedding_leiden_res=joint_embedding_leiden_res,
        run_clust_label_transfer=run_clust_label_transfer,
        confusion_matrix_cluster_min_value=confusion_matrix_cluster_min_value,
        confusion_matrix_cluster_max_value=confusion_matrix_cluster_max_value,
        confusion_matrix_cluster_resolution=confusion_matrix_cluster_resolution,
        qry_only_cluster_threshold=qry_only_cluster_threshold,
        ref_only_cluster_threshold=ref_only_cluster_threshold,
        save_integrator=save_integrator,
        save_adata_comb=save_adata_comb,
        **kwargs,
    )

    logger.info("DONE WITH ALLCOOLS INTEGRATION")

    if backup_to_spatialdata: 
        logger.info("Backing up the AnnData object to the spatialdata object")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.utilities.sd_utils import _gen_keys,_backup_adata
            from spida._constants import TABLE_KEY
            KEYS = _gen_keys(prefix_name, exp_name, reg_name)
            _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}")
    
    if plot:
        logger.info("Plotting the results of the ALLCools integration")
        exp = adata.uns.get("experiment")
        seg_name = adata.uns.get("segmentation")
        donor = adata.uns.get("donor")
        ctx.invoke(
            plot_allcools_integration,
            exp_name=exp,
            seg_name=seg_name,
            donor=donor,
            anndata_store_path=anndata_store_path,
            annotations_store_path=annotations_store_path,
            output_path=image_path,
            plot_joint_embeddings=run_joint_embeddings,
            plot_c2c_transfer=run_clust_label_transfer,
            ref_cell_type_column=rna_cell_type_column,
        )

@cli.command(
    name='plot-allcools-integration', 
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('exp_name', type=click.STRING)
@click.argument('seg_name', type=click.STRING)
@click.argument('donor', type=click.STRING)
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.option(
    '--output_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the output directory for results. Defaults to None.'
)
@click.option("--plot_joint_embeddings", type=bool, is_flag=True, default=False, help='Whether to plot the joint embeddings. Defaults to False.')
@click.option("--plot_c2c_transfer", type=bool, is_flag=True, default=False, help='Whether to plot the cluster-to-cluster label transfer results. Defaults to False.')
@click.option("--ref_cell_type_column", type=click.STRING, default="supercluster_name", help='Column name in the AnnData object for RNA cell types. Defaults to "supercluster_name".')
@click.pass_context
def plot_allcools_integration(
    ctx,
    exp_name: str,
    seg_name: str,
    donor: str,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    output_path: str = None,
    plot_joint_embeddings:bool = False,
    plot_c2c_transfer:bool = False,
    ref_cell_type_column:str = "supercluster_name",
):
    """
    Plot the results of ALLCools integration for a given EXP_NAME, SEG_NAME, and DONOR.

    Arguments:\n
    exp_name (str): Name of the experiment.\n
    seg_name (str): Name of the segmentation.\n
    donor (str): Name of the donor.\n
    """
    import anndata as ad # type: ignore
    from spida.pl import plot_allcools_joint_embeddings, plot_allcools_spatial_annot
    from matplotlib.backends.backend_pdf import PdfPages

    # TODO: Click configuration files (handled within the cli group function to get all of these filepaths)
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

    if output_path is None:
        image_store = os.getenv("IMAGE_STORE_PATH")
        output_path = Path(
            f"{image_store}/{exp_name}/{seg_name}/region_{donor}/allcools_integration.pdf"
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)    
    

    qry_path = Path(f"{anndata_store_path}/{exp_name}/{seg_name}/adata_{donor}.h5ad")
    qry_adata = ad.read_h5ad(qry_path)

    adata_comb_path = f"{annotations_store_path}/{exp_name}/{seg_name}/allcools/{donor}"
    adata_comb = ad.read_h5ad(f"{adata_comb_path}/adata_comb_rna_merfish.h5ad")

    pdf_file = PdfPages(output_path)
    if plot_joint_embeddings:
        plot_allcools_joint_embeddings(adata_comb, exp_name=exp_name, seg_name=seg_name, donor=donor, pdf_file=pdf_file)
    plot_allcools_spatial_annot(qry_adata, exp_name=exp_name, seg_name=seg_name, donor=donor, pdf_file=pdf_file)
    if plot_c2c_transfer:
        from spida.pl.I_plots import plot_allcools_c2c
        ref_group = adata_comb.uns['c2c_allcools_integration_results'].get("ref_group", None)
        ref_group = ref_group['ref_group'] if ref_group is not None else None
        query_group = adata_comb.uns['c2c_allcools_integration_results'].get("qry_group", None)
        query_group = query_group['qry_group'] if query_group is not None else None
        conf_mat = adata_comb.uns['c2c_allcools_integration_results'].get("confusion_matrix", None)
        plot_allcools_c2c(adata_comb, qry_adata, ref_group, query_group, conf_mat, ref_cell_type_column=ref_cell_type_column, pdf_file=pdf_file)
    pdf_file.close()

@cli.command(
    name='allcools-integration-experiment',
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('exp_name', type=click.STRING)
@click.argument('prefix_name', type=click.STRING)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@click.option('--suffix', type=click.STRING, default='_filt', help='Suffix for the keys in the spatialdata object. Defaults to "_filt".')
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.pass_context
def allcools_integration_experiment(
    ctx,
    exp_name: str,
    prefix_name: str,
    ref_path: str,
    suffix:str = "_filt",
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Run ALLCools integration for an entire EXP_NAME with PREFIX_NAME.

    Arguments:\n
    exp_name (str): Name of the experiment.\n
    prefix_name (str): Prefix for the keys in the spatialdata object.\n
    ref_path (str): Path to the reference RNA AnnData object .\n
    to None.
    **kwargs: Additional keyword arguments for ALLCools integration.\n
    """
    extra_args = ctx.args
    if extra_args:
        click.echo(f"Received extra arguments: {extra_args}")
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        click.echo(f"Parsed {key} = {value};  {type(value)}")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.allcools import allcools_integration_experiment as func
    func(
        exp_name,
        prefix_name,
        ref_path,
        suffix,
        anndata_store_path,
        annotations_store_path,
        **kwargs,
    )

@cli.command(
    name='mmc-setup', 
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('ref_path', type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False))
@click.argument('identifier', type=click.STRING)
@click.option('--hierarchy_list', type=parse_list, help='Comma-separated list of hierarchy levels (e.g., "level1,level2,level3")')
@click.option(
    '--gene_names_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False),
    default=None,
    help='Path to gene_names.txt file. Defaults to None.'
)
@click.option(
    '--gene_name_mapping_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False),
    default=None,
    help='Path to gene_name_mapping.json file. Defaults to None.'
)
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--ref_norm',
    type=click.STRING,
    default='log2CPM',
    help='Normalization method for the reference AnnData.X. Defaults to "log2CPM".'
)
@click.option(
    '--n_cpu',
    type=int,
    default=8,
    help='Number of CPUs to use. Defaults to 8.'
)
@click.option(
    '--n_valid',
    type=int,
    default=10,
    help='Number of valid markers. Defaults to 10.'
)
@click.option(
    '--n_per_utility',
    type=int,
    default=10,
    help='Number of markers per utility. Defaults to 10.'
)
@click.pass_context
def mmc_setup(
    ctx,
    ref_path: str,
    identifier: str,
    hierarchy_list: list = ["Class", "Subclass", "Group"],
    gene_names_path: str = None,
    gene_name_mapping_path: str = None,
    mmc_store_path: str = None,
    ref_norm: str = "log2CPM",
    n_cpu: int = 8,
    n_valid: int = 10,
    n_per_utility: int = 10,
    **kwargs,
):
    """
    Setup function for MapMyCells integration using refactored architecture.

    Arguments:
    ref_path: Path to the reference single-cell RNA-seq AnnData file
    identifier: Unique identifier for this reference (e.g., 'mouse_motor_cortex')

    Options:
    --hierarchy_list: Comma-separated list of hierarchy levels (e.g., "level1,level2,level3")
    --gene_names_path: Path to gene_names.txt file (auto-detected if not provided)
    --gene_name_mapping_path: Path to gene_name_mapping.json file (auto-detected if not provided)
    --mmc_store_path: Path to MapMyCells store (uses MMC_DIR env var if not provided)
    --ref_norm: Normalization method (default: log2CPM)
    --n_cpu: Number of CPUs (default: 8)
    --n_valid: Number of valid markers (default: 10)
    --n_per_utility: Number of markers per utility (default: 10)
    """
    extra_args = ctx.args
    if extra_args:
        logger.warning(f"Received extra arguments: {extra_args}")
    logger.info(f"Hierarchy list: {hierarchy_list}, type: {type(hierarchy_list)}")
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc_cli import mmc_setup_cli
    
    result = mmc_setup_cli(
        ref_path=ref_path,
        identifier=identifier,
        hierarchy_list=hierarchy_list,
        gene_names_path=gene_names_path,
        gene_name_mapping_path=gene_name_mapping_path,
        mmc_store_path=mmc_store_path,
        ref_norm=ref_norm,
        n_cpu=n_cpu,
        n_valid=n_valid,
        n_per_utility=n_per_utility,
    )
    
    logger.info(f"Setup complete. Results: {result}")
    click.echo(f"Setup complete. Results: {result}")


@cli.command(
    name="mmc-annotation-region",
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('exp_name', type=click.STRING)
@click.argument('reg_name', type=click.STRING)
@click.argument('prefix_name', type=click.STRING)
@click.argument('identifier', type=click.STRING)
@click.option('--suffix', type=click.STRING, default='_filt', help='Suffix for the keys in the spatialdata object. Defaults to "_filt".')
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--anndata_store',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotation_store',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.option(
    '--zarr_store',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the zarr store of spatialdata objects. Defaults to None.'
)
@click.option(
    '--n_cpu',
    type=int,
    default=1,
    help='Number of CPUs to use. Defaults to 1.'
)
@click.option(
    '--bootstrap_factor',
    type=float,
    default=0.8,
    help='Bootstrap factor for annotation. Defaults to 0.8.'
)
@click.option(
    '--bootstrap_iterations',
    type=int,
    default=100,
    help='Number of bootstrap iterations. Defaults to 100.'
)
@click.option(
    '--rng_seed',
    type=int,
    default=13,
    help='Random number generator seed. Defaults to 13.'
)
@click.option("--filter_annotations", type=bool, is_flag=True, default=False, help='Whether to filter annotations based on the number of markers. Defaults to False.')
@click.option("--plot", type=bool, is_flag=True, default=False, help='Whether to plot the annotation results. Defaults to False.')
@click.option("--plot_name", type=click.STRING, default="mmc_annotation_plots", help='Name for the annotation plots. Defaults to "mmc_annotation_plots".')
@click.option("--palette_path", type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False), default=None, help='Path to a custom palette file for plotting. Defaults to None.')
@click.option("--image_store", type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True), default=None, help='Path to the image store for saving plots. Defaults to None..')
@click.pass_context
def mmc_annotation_region(
    ctx, 
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    identifier: str,
    suffix: str = "_filt",
    mmc_store_path: str = None,
    anndata_store: str = None,
    annotation_store: str = None,
    zarr_store: str = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    filter_annotations: bool = False,
    plot: bool = False,
    plot_name: str = "mmc_annotation_plots",
    palette_path: str = None,
    image_store: str = None,
    **kwargs,
):
    """
    Run MapMyCells annotation on a given experiment and region using refactored architecture.

    Arguments:
    exp_name: Name of the experiment
    reg_name: Name of the region
    prefix_name: Prefix for keys in the spatialdata object
    identifier: Reference identifier (e.g., 'mouse_motor_cortex')

    Options:
    --mmc_store_path: Path to MapMyCells store (uses MMC_DIR env var if not provided)
    --anndata_store_path: Path to AnnData store (uses ANNDATA_STORE_PATH env var if not provided)
    --annotations_store_path: Path to annotations store (uses ANNOTATIONS_STORE_PATH env var if not provided)
    --n_cpu: Number of CPUs (default: 1)
    --bootstrap_factor: Bootstrap factor (default: 0.8)
    --bootstrap_iterations: Bootstrap iterations (default: 100)
    --rng_seed: Random seed (default: 13)
    --filter_annotations: Whether to filter annotations based on the number of markers (default: False)
    --plot: Whether to plot the annotation results (default: False)
    --palette_path: Path to a custom palette file for plotting (default: None)
    --image_store: Path to the image store for saving plots (default: None)
    """
    extra_args = ctx.args
    # if extra_args:
    #     logger.warning(f"Received extra arguments: {extra_args}")
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        logger.info(f"Parsed {key} = {value};  {type(value)}")
    
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc_cli import mmc_annotation_region_cli
    
    result = mmc_annotation_region_cli(
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        suffix=suffix,
        identifier=identifier,
        mmc_store_path=mmc_store_path,
        anndata_store_path=anndata_store,
        annotations_store_path=annotation_store,
        zarr_store_path=zarr_store,
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
        filter_annot=filter_annotations,
        plot=plot,
        plot_name=plot_name,
        palette_path=palette_path,
        image_store_path=image_store,
        **kwargs,
    )
    
    logger.info(f"Annotation complete for {exp_name}/{reg_name}")
    click.echo(f"Annotation complete for {exp_name}/{reg_name}")

@cli.command(
    name="mmc-annotation-experiment", 
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument('exp_name', type=click.STRING)
@click.argument('prefix_name', type=click.STRING)
@click.argument('identifier', type=click.STRING)
@click.option(
    '--suffix',
    type=click.STRING,
    default='_filt',
    help='Suffix for keys in spatialdata. Defaults to "_filt".'
)
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.option(
    '--n_cpu',
    type=int,
    default=1,
    help='Number of CPUs to use. Defaults to 1.'
)
@click.option(
    '--bootstrap_factor',
    type=float,
    default=0.8,
    help='Bootstrap factor for annotation. Defaults to 0.8.'
)
@click.option(
    '--bootstrap_iterations',
    type=int,
    default=100,
    help='Number of bootstrap iterations. Defaults to 100.'
)
@click.option(
    '--rng_seed',
    type=int,
    default=13,
    help='Random number generator seed. Defaults to 13.'
)
@click.option("--filter_annotations", type=bool, is_flag=True, default=True, help='Whether to filter annotations based on the number of markers. Defaults to True.')
@click.option("--plot", type=bool, is_flag=True, default=False, help='Whether to plot the annotation results. Defaults to False.')
@click.option("--plot_name", type=click.STRING, default="mmc_annotation_plots", help='Name for the annotation plots. Defaults to "mmc_annotation_plots".')
@click.option("--palette_path", type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False), default=None, help='Path to a custom palette file for plotting. Defaults to None.')
@click.option("--image_store", type=click.Path(exists=True, file_okay=False, path_type=str, dir_okay=True), default=None, help='Path to the image store for saving plots. Defaults to None..')
@click.pass_context
def mmc_annotation_experiment(
    ctx,
    exp_name: str,
    prefix_name: str,
    identifier: str,
    suffix: str = "_filt",
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    n_cpu: int = 1,
    bootstrap_factor: float = 0.8,
    bootstrap_iterations: int = 100,
    rng_seed: int = 13,
    filter_annotations: bool = True,
    plot: bool = False,
    plot_name: str = "mmc_annotation_plots",
    palette_path: str = None,
    image_store: str = None,
    **kwargs,
):
    """
    Run MapMyCells annotation for an entire experiment using refactored architecture.

    Arguments:
    exp_name: Name of the experiment
    prefix_name: Prefix for keys in the spatialdata object
    identifier: Reference identifier (e.g., 'mouse_motor_cortex')

    Options:
    --suffix: Suffix for keys in spatialdata (default: "_filt")
    --mmc_store_path: Path to MapMyCells store (uses MMC_DIR env var if not provided)
    --anndata_store_path: Path to AnnData store (uses ANNDATA_STORE_PATH env var if not provided)
    --annotations_store_path: Path to annotations store (uses ANNOTATIONS_STORE_PATH env var if not provided)
    --n_cpu: Number of CPUs (default: 1)
    --bootstrap_factor: Bootstrap factor (default: 0.8)
    --bootstrap_iterations: Bootstrap iterations (default: 100)
    --rng_seed: Random seed (default: 13)
    """
    extra_args = ctx.args
    if extra_args:
        logger.warning(f"Received extra arguments: {extra_args}")
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc_cli import mmc_annotation_experiment_cli
    
    result = mmc_annotation_experiment_cli(
        exp_name=exp_name,
        prefix_name=prefix_name,
        suffix=suffix,
        identifier=identifier,
        mmc_store_path=mmc_store_path,
        anndata_store_path=anndata_store_path,
        annotations_store_path=annotations_store_path,
        n_cpu=n_cpu,
        bootstrap_factor=bootstrap_factor,
        bootstrap_iterations=bootstrap_iterations,
        rng_seed=rng_seed,
        filter_annot=filter_annotations,
        plot=plot,
        plot_name=plot_name,
        palette_path=palette_path,
        image_store_path=image_store,
    )
    
    logger.info(f"Annotation complete for experiment {exp_name}")
    click.echo(f"Annotation complete for experiment {exp_name}")

#### MOSCOT INTEGRATION ####
@cli.command(
    name="moscot-integration",
    cls=RichCommand,
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
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
