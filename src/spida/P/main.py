import os
import glob
from pathlib import Path
from dotenv import load_dotenv  # type: ignore
import warnings
import logging

from spida.utilities.sd_utils import (
    _gen_keys,
    _region_to_donor,
    _write_adata,
    _backup_adata,
    _get_adata,
    _assign_new_table,
)
from spida._constants import TABLE_KEY, IMAGE_KEY, SHAPES_KEY  # type: ignore

load_dotenv()
logger = logging.getLogger(__package__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


def filter_cells_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
):
    """
    Filter cells for a specific region in an experiment.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file. If None, uses a default path.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    Raises:
    NotImplementedError: If plotting is requested but not implemented.
    """
    from spida.P.filtering import run_filtering

    logger.info(
        "FILTERING CELLS, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    # default cutoffs path
    if cutoffs_path is None:
        cutoffs_path = os.getenv(
            "DEF_CUTOFFS_PATH",
            "/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json",
        )

    # determining donor from region name
    donor_name = _region_to_donor(reg_name)

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    adata = _get_adata(exp_name, reg_name, prefix_name)
    # Run the filtering
    adata = run_filtering(
        adata, exp_name, reg_name, prefix_name, donor_name, cutoffs_path
    )
    # backup the AnnData object
    # _assign_new_table(exp_name, reg_name, adata, KEYS[TABLE_KEY], suffix="_filt") # double the storage but allows for iteration on filts.
    _backup_adata(exp_name, reg_name, adata, KEYS[TABLE_KEY])

    logger.info(
        f"Passed QC Cells: {adata.obs['pass_qc'].sum()} out of {adata.n_obs} total cells"
    )
    logger.info("DONE")
    if plot:
        plot_filtering_region(exp_name, reg_name, prefix_name, image_path=image_path)


def plot_filtering_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
):
    """
    Plot the filtering results for a specific region in an experiment.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    """
    from spida.pl import plot_filtering
    from matplotlib.backends.backend_pdf import PdfPages

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix)

    if image_path is None:
        image_store = os.getenv(
            "IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images"
        )
        image_path = Path(
            f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_filt.pdf"
        )
        image_path.parent.mkdir(parents=True, exist_ok=True)

    pdf_file = PdfPages(image_path)
    fig, ax = plot_filtering(adata, exp_name, reg_name, prefix_name)
    pdf_file.savefig(fig)
    pdf_file.close()


def filter_cells_all(
    exp_name: str,
    prefix_name: str,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
):
    """
    Filter cells for all regions in an experiment (calls filter_cells_all).
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file. If None, uses a default path.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions:
        reg_name = reg.split("/")[-1]
        filter_cells_region(
            exp_name,
            reg_name,
            prefix_name,
            cutoffs_path=cutoffs_path,
            plot=plot,
            image_path=image_path,
        )


def write_adata(
    exp_name,
    reg_name: str = None,
    prefix_names: list[str] = None,
    output_path: Path = "/ceph/cephatlas/aklein/bican/data/anndatas/",
):
    """ """
    print(exp_name, reg_name, prefix_names, output_path)
    # Iterating over all prefixes
    for p in prefix_names:
        logger.info("PREFIX: %s" % p)
        # if a specific region is provided, write only that region
        if reg_name:
            fout = f"{output_path}/{exp_name}/{p}"
            _write_adata(exp_name, reg_name, p, fout)
        # if no region is provided, write all regions
        else:
            zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
            region_list = glob.glob(f"{zarr_store}/{exp_name}/region_*")
            for reg in region_list:
                rname = reg.split("/")[-1]
                fout = f"{output_path}/{exp_name}/{p}"
                _write_adata(exp_name, rname, p, fout)


def setup_adata_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str = "",
    plot=False,
    image_path: Path = None,
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
    from spida.P.setup_adata import run_setup

    logger.info(
        "SETTING UP ADATA, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    # determining donor from region name
    donor_name = _region_to_donor(reg_name)
    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    adata = _get_adata(
        exp_name, reg_name, prefix_name
    )  # Use the filtered data if available
    # Run the setup
    adata = run_setup(adata, exp_name, reg_name, prefix_name, donor_name)
    # backup the adata object
    try:
        _assign_new_table(
            exp_name, reg_name, adata, KEYS[TABLE_KEY], suffix=suffix
        )  # double the storage but allows for iteration on filts.
    except Exception as e:
        logger.warning("Failed to assign new table, using backup instead")
        logger.warning(f"Error: {e}", exc_info=True)
        _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}")

    logger.info("DONE SETUP")
    if plot:
        plot_setup_region(
            exp_name, reg_name, prefix_name, image_path=image_path, suffix=suffix
        )


def setup_adata_all(
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
):
    """
    Setup the AnnData objects for all regions in an experiment.
    This function iterates over all regions and calls setup_adata for each region.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    for reg in regions:
        reg_name = reg.split("/")[-1]
        setup_adata_region(
            exp_name,
            reg_name,
            prefix_name,
            suffix=suffix,
            plot=plot,
            image_path=image_path,
        )


def plot_setup_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
):
    """
    Plot the setup results for a specific region in an experiment.
    This function generates a PDF with the setup plots.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    """
    from spida.pl import plot_setup
    from matplotlib.backends.backend_pdf import PdfPages

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix)

    if image_path is None:
        image_store = os.getenv(
            "IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images"
        )
        image_path = Path(
            f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_setup.pdf"
        )
        image_path.parent.mkdir(parents=True, exist_ok=True)

    pdf_file = PdfPages(image_path)
    plot_setup(adata, exp_name, reg_name, prefix_name, pdf_file=pdf_file)
    pdf_file.close()


def remove_doublets_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    **model_kwargs,
):
    """
    Remove doublets from the AnnData object for a specific region in an experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    threshold (float, optional): Threshold for doublet detection. Defaults to 0.5.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """
    from spida.P.scvi_toolkit import identify_doublets  # , remove_doublets

    logger.info(
        "REMOVING DOUBLETS, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    adata = _get_adata(
        exp_name, reg_name, prefix_name, suffix=suffix
    )  # Use the filtered data if available

    # Identify doublets
    adata = identify_doublets(adata, threshold=threshold)
    logger.info(
        f"Doublets identified: {adata.obs['doublet_bool'].sum()} out of {adata.n_obs} total cells"
    )

    # backup the adata object
    _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}")

    logger.info("DONE REMOVING DOUBLETS")

    if plot:  # There has to be a better way to do this!
        logger.info("PLOTTING DOUBLETS")
        # Importing here for efficiency
        from matplotlib.backends.backend_pdf import PdfPages
        from spida.pl import plot_doublets
        import spatialdata as sd
        from spida.pl import plot_example_doublets

        KEYS = _gen_keys(prefix_name, exp_name, reg_name)

        zarr_store = os.getenv("ZARR_STORAGE_PATH")
        zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            sdata = sd.read_zarr(zarr_path)

        if image_path is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
            image_path = Path(
                f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/SOLO_doublets.png"
            )
            image_path.parent.mkdir(parents=True, exist_ok=True)

        plot_doublets(adata, image_path)
        image_path = image_path.parent / "SOLO_examples.pdf"
        with PdfPages(image_path) as pdf:
            plot_example_doublets(
                sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[TABLE_KEY], pdf_file=pdf
            )

    # Remove doublets
    # adata = remove_doublets(adata)


def remove_doublets_all(
    exp_name: str,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    **model_kwargs,
):
    """
    Remove doublets from all regions in an experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    threshold (float, optional): Threshold for doublet detection. Defaults to 0.5.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    """

    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    logger.info("REMOVING DOUBLETS FOR ALL REGIONS IN EXPERIMENT %s" % exp_name)
    logger.info("ZARR_PATH: %s" % zarr_store)

    for reg in regions:
        reg_name = reg.split("/")[-1]
        remove_doublets_region(
            exp_name,
            reg_name,
            prefix_name,
            threshold=threshold,
            plot=plot,
            image_path=image_path,
            suffix=suffix,
            **model_kwargs,
        )


def resolvi_cluster_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    model_save_path : Path = None,
    trained : bool = False, 
    **model_kwargs,
):
    """
    Perform RESOLVI clustering on the AnnData object for a specific region in an experiment.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    suffix (str, optional): Suffix for the AnnData object. Defaults to "".
    plot (bool, optional): Whether to plot the results. Defaults to False.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    """
    from spida.P.scvi_toolkit import resolvi_cluster
    from spida.P.setup_adata import _calc_embeddings

    logger.info(
        "RESOLVI CLUSTERING, EXPERIMENT %s, REGION %s, PREFIX %s"
        % (exp_name, reg_name, prefix_name)
    )

    model_kwargs = model_kwargs['model_kwargs']
    for key, val in model_kwargs.items():
        logger.info(f"Model kwargs: {key} = {val}")

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    adata = _get_adata(
        exp_name, reg_name, prefix_name, suffix=suffix
    )  # Use the filtered data if available

    logger.info("arguments passing into resolvi")
    logger.info(f"max_epochs: {model_kwargs.get('max_epochs', 200)}")
    logger.info(f"layer: {model_kwargs.get('layer', 'raw')}")
    logger.info(f"batch_key: {model_kwargs.get('batch_key', 'donor')}")
    logger.info(f"categorical_covariates: {model_kwargs.get('categorical_covariates', None)}")
    logger.info(f"model_save_path: {model_save_path}")
    logger.info(f"model_save_ext: {f'{exp_name}_{reg_name}_{prefix_name}'}")
    
    # Perform RESOLVI clustering
    adata = resolvi_cluster(
        adata,
        layer=model_kwargs.get("layer", "raw"),
        batch_key=model_kwargs.get("batch_key", "donor"),
        categorical_covariates=model_kwargs.get("categorical_covariates", None),
        model_save_path=model_save_path,
        model_save_ext=f"{exp_name}_{reg_name}_{prefix_name}",
        **model_kwargs,
    )

    logger.info("DONE WITH RESOLVI, NOW CLUSTERING")

    # _calc_embeddings(adata, key_added="base_")
    logger.info("embedding resolvi")
    _calc_embeddings(adata, use_rep="X_resolVI", key_added="resolvi_")
    logger.info("embedding generated expression")
    _calc_embeddings(adata, layer="generated_expression", key_added="corr_")

    # backup the adata object
    _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}")

    logger.info("DONE RESOLVI CLUSTERING")
    if plot:
        plot_resolvi_region(
            exp_name, reg_name, prefix_name, image_path=image_path, suffix="_filt"
        )


def resolvi_cluster_all(
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    model_save_path : Path = None, 
    **model_kwargs,
):
    """
    Perform RESOLVI clustering on all regions in an experiment.
    Parameters:
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    suffix (str, optional): Suffix for the AnnData object. Defaults to "".
    plot (bool, optional): Whether to plot the results. Defaults to False.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    """
    # Getting the regions for the experiment
    zarr_store = os.getenv("ZARR_STORAGE_PATH", "/data/aklein/bican_zarr")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")

    for reg in regions:
        reg_name = reg.split("/")[-1]
        resolvi_cluster_region(
            exp_name,
            reg_name,
            prefix_name,
            suffix=suffix,
            plot=plot,
            image_path=image_path,
            model_save_path=model_save_path,
            **model_kwargs,
        )


def plot_resolvi_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
):
    """
    Plot the setup results for a specific region in an experiment.
    This function generates a PDF with the setup plots.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    """
    from spida.pl import plot_resolvi
    from matplotlib.backends.backend_pdf import PdfPages

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix)

    if image_path is None:
        image_store = os.getenv(
            "IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images"
        )
        image_path = Path(
            f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/pixi_resolvi.pdf"
        )
        image_path.parent.mkdir(parents=True, exist_ok=True)

    pdf_file = PdfPages(image_path)
    plot_resolvi(adata, exp_name, reg_name, prefix_name, pdf_file=pdf_file)
    pdf_file.close()


def combine_datasets(
    experiment_list: list[str],
    prefix_name: str,
    suffix: str = "_filt",
    lab_name: str = "salk",
    project_name: str = "BG_salk",
    zarr_path: Path = None,
    output_path: Path = None,
    table_keys: list[str] = None,
):
    """
    Aggregate spatial data from multiple experiments into a single SpatialData object.

    Parameters:
        experiment_list (list[str]): List of experiment names to aggregate
        prefix_name (str): Prefix for the keys in the aggregated data
        suffix (str): Suffix to append to the table keys
        lab_name (str): Name of the lab, used for renaming experiments and regions
        project_name (str): Name of the project, used for naming the aggregated data
        zarr_path (Path): Path to the Zarr store containing the experiments
        output_path (Path): Path to save the aggregated AnnData object
        table_keys (list[str]): List of table keys to concatenate from the SpatialData object
    """
    from spida.P.project_level import aggregate_experiments, concatenate_tables

    if zarr_path is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
        zarr_path = Path(zarr_store)

    if output_path is None:
        output_path = Path(os.getenv("ANNDATA_STORE_PATH"))
        output_path = output_path / f"{project_name}_{lab_name}_{prefix_name}{suffix}.h5ad"

    # aggregate the spatialdata objects
    sdata = aggregate_experiments(
        zarr_path, experiment_list, prefix_name, suffix, lab_name, project_name
    )

    # create the combined anndata object
    adata = concatenate_tables(sdata, table_keys=table_keys, output_path=output_path)

    # Write the combined AnnData object to the spatialdata store
    sdata["combined_table"] = adata
    sdata.write_element("combined_table")


def setup_dataset(
    dataset_name : str,
    scale : bool = False,
    anndata_path : Path = None,
): 
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    from spida.P.setup_adata import combined_setup
    import anndata as ad

    logger.info(f"SETUP ON {dataset_name}")
    if anndata_path is None: 
        anndata_path = Path(os.getenv("ANNDATA_STORE_PATH"))
    anndata_path = anndata_path / dataset_name
    adata = ad.read_h5ad(anndata_path)    
    adata = combined_setup(adata, scale=scale)
    adata.write_h5ad(anndata_path)
    logger.info(f"DONE SETUP ON {dataset_name}")

def resolvi_dataset(
    dataset_name:str,
    anndata_path:Path = None,
    model_save_path:Path = None, 
    trained:bool = False,
    **model_kwargs
    ): 
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    from spida.P.scvi_toolkit import resolvi_cluster
    from spida.P.setup_adata import _calc_embeddings
    from spida.pl import plot_dataset
    import anndata as ad

    model_kwargs = model_kwargs['model_kwargs']
    
    logger.info(f"RESOLVI ON {dataset_name}")
    if anndata_path is None: 
        anndata_path = Path(os.getenv("ANNDATA_STORE_PATH"))
    anndata_path = anndata_path / dataset_name
    adata = ad.read_h5ad(anndata_path)    
    
    
    adata = resolvi_cluster(
        adata,
        layer=model_kwargs.get("layer", "counts"),
        batch_key=model_kwargs.get("batch_key", "dataset_id"),
        categorical_covariates=model_kwargs.get("categorical_covariates", None),
        model_save_path=model_save_path,
        model_save_ext=dataset_name,
        trained=trained,
        **model_kwargs
    )

    logger.info("DONE WITH RESOLVI, NOW CLUSTERING")

    # _calc_embeddings(adata, key_added="base_")
    logger.info("embedding resolvi")
    _calc_embeddings(adata, use_rep="X_resolvi", key_added="resolvi_")
    logger.info("embedding generated expression")
    _calc_embeddings(adata, layer="generated_expression", key_added="corr_")

    adata.write_h5ad(anndata_path)
    logger.info(f"DONE RESOLVI ON {dataset_name}")

    logger.info("PLOTTING DATASET")
    
    if image_path is None:
        image_store = os.getenv(
            "IMAGE_STORE_PATH", "/ceph/cephatlas/aklein/bican/images"
        )
        image_path = Path(
            f"{image_store}/{dataset_name}"
        )
        image_path.mkdir(parents=True, exist_ok=True)

    plot_dataset(adata, save_path=image_path, show=False)
    return 0
