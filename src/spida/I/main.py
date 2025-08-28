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
from rich_click import RichCommand, RichGroup # type: ignore
from spida.utilities.script_utils import parse_click_kwargs, JSONParam
from spida.settings import configure_logging_for_runtime

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

@click.group(cls=RichGroup)
def cli():
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
@click.option("--run_clust_label_transfer", type=bool, is_flag=True, default=True, help='Whether to run cluster-to-cluster label transfer on the integrated data. Defaults to True.')
@click.option("--confusion_matrix_cluster_min_value", type=float, default=0.25, help='Minimum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.25.')
@click.option("--confusion_matrix_cluster_max_value", type=float, default=0.9, help='Maximum value for the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 0.9.')
@click.option("--confusion_matrix_cluster_resolution", type=float, default=1.5, help='Resolution for clustering the confusion matrix when performing cluster-to-cluster label transfer. Defaults to 1.5.')
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
    run_clust_label_transfer: bool = True,
    confusion_matrix_cluster_min_value: float = 0.25,
    confusion_matrix_cluster_max_value: float = 0.9,
    confusion_matrix_cluster_resolution: float = 1.5,
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
        run_joint_embeddings=run_joint_embeddings,
        run_clust_label_transfer=run_clust_label_transfer,
        confusion_matrix_cluster_min_value=confusion_matrix_cluster_min_value,
        confusion_matrix_cluster_max_value=confusion_matrix_cluster_max_value,
        confusion_matrix_cluster_resolution=confusion_matrix_cluster_resolution,
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
@click.argument('heirarchy_list', type=click.STRING, nargs=-1)
@click.argument('BRAIN_REGION', type=click.STRING)
@click.argument('CODEBOOK', type=click.STRING)
@click.option(
    '--codebook_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=False),
    default=None,
    help='Path to the codebook file. Defaults to None.'
)
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--ref_norm',
    type=click.STRING,
    default='log2CPM',
    help='Normalization method for the reference AnnData.X. Defaults to "log2CPM".'
)
@click.pass_context
def mmc_setup(
    ctx,
    ref_path: str,
    heirarchy_list: list,
    BRAIN_REGION: str,
    CODEBOOK: str,
    codebook_path: str = None,
    mmc_store_path: str = None,
    ref_norm: str = "log2CPM",
    **kwargs,
):
    """
    Setup function for MapMyCells integration.

    Parameters:
    ref_path (str): Path to the reference data.
    heirarchy_list (list): List of hierarchy levels for the annotation.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    codebook_path (str, optional): Path to the codebook file. Defaults to None.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    ref_norm (str, optional): Normalization method for the reference AnnData.X. Defaults to "log2CPM".
    **kwargs: Additional keyword arguments.
    """
    extra_args = ctx.args
    if extra_args:
        click.echo(f"Received extra arguments: {extra_args}")
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        click.echo(f"Parsed {key} = {value};  {type(value)}")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import mmc_setup as func
    func(
        ref_path,
        BRAIN_REGION,
        CODEBOOK,
        codebook_path,
        heirarchy_list,
        mmc_store_path,
        ref_norm,
        **kwargs,
    )


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
@click.argument('BRAIN_REGION', type=click.STRING)
@click.argument('CODEBOOK', type=click.STRING)
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.pass_context
def mmc_annotation_region(
    ctx, 
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    BRAIN_REGION: str,
    CODEBOOK: str,
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Run MapMyCells annotation on a given experiment and region.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    extra_args = ctx.args
    if extra_args:
        click.echo(f"Received extra arguments: {extra_args}")
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        click.echo(f"Parsed {key} = {value};  {type(value)}")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import mmc_annotation_region as func
    func(
        exp_name,
        reg_name,
        prefix_name,
        BRAIN_REGION,
        CODEBOOK,
        mmc_store_path,
        anndata_store_path,
        annotations_store_path,
        **kwargs,
    )

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
@click.argument('BRAIN_REGION', type=click.STRING)
@click.argument('CODEBOOK', type=click.STRING)
@click.option(
    '--mmc_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of MapMyCells markers. Defaults to None.'
)
@click.option(
    '--anndata_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of AnnData objects. Defaults to None.'
)
@click.option(
    '--annotations_store_path',
    type=click.Path(exists=True, file_okay=True, path_type=str, dir_okay=True),
    default=None,
    help='Path to the store of annotation specific files. Defaults to None.'
)
@click.pass_context
def mmc_annotation_experiment(
    ctx,
    exp_name: str,
    prefix_name: str,
    BRAIN_REGION: str,
    CODEBOOK: str,
    mmc_store_path: str = None,
    anndata_store_path: str = None,
    annotations_store_path: str = None,
    **kwargs,
):
    """
    Run MapMyCells annotation for an entire experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    BRAIN_REGION (str): Brain region for the annotation.
    CODEBOOK (str): Codebook for the annotation.
    mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
    anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
    annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
    **kwargs: Additional keyword arguments.
    """
    extra_args = ctx.args
    if extra_args:
        click.echo(f"Received extra arguments: {extra_args}")
    kwargs = parse_click_kwargs(extra_args)
    for key, value in kwargs.items():
        click.echo(f"Parsed {key} = {value};  {type(value)}")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spida.I.mmc import mmc_annotation_experiment as func
    func(
        exp_name,
        prefix_name,
        BRAIN_REGION,
        CODEBOOK,
        mmc_store_path,
        anndata_store_path,
        annotations_store_path,
        **kwargs,
    )

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
