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
from spida._constants import TABLE_KEY, IMAGE_KEY, SHAPES_KEY, POINTS_KEY  # type: ignore

load_dotenv()
logger = logging.getLogger(__package__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


def filter_cells_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    seg_fam:str = None,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None
):
    """
    Filter cells for a specific region in an experiment.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    seg_fam (str, optional): Which preset of column names to use when running filtering (default is None which uses a constant map)
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
        cutoffs_path = os.getenv("CUTOFFS_PATH",)

    # determining donor from region name
    donor_name = _region_to_donor(reg_name)

    # Get KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    # Get the AnnData object
    logger.info("ZARR STORE: %s" % zarr_store)
    adata = _get_adata(exp_name, reg_name, prefix_name, zarr_store=zarr_store)
    # Run the filtering
    adata = run_filtering(
        adata=adata,
        exp_name=exp_name,
        reg_name=reg_name,
        prefix_name=prefix_name,
        donor_name=donor_name,
        seg_fam=seg_fam, 
        cutoffs_path=cutoffs_path
    )
    # backup the AnnData object
    # _assign_new_table(exp_name, reg_name, adata, KEYS[TABLE_KEY], suffix="_filt") # double the storage but allows for iteration on filts.
    _backup_adata(exp_name, reg_name, adata, KEYS[TABLE_KEY], zarr_store=zarr_store)

    logger.info(
        f"Passed QC Cells: {adata.obs['pass_qc'].sum()} out of {adata.n_obs} total cells"
    )
    logger.info("DONE")
    if plot:
        plot_filtering_region(exp_name, reg_name, prefix_name, image_path=image_path, image_store=image_store, zarr_store=zarr_store)


def plot_filtering_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    image_store: Path = None,
    suffix: str = "",
    zarr_store: str | Path | None = None
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

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix, zarr_store=zarr_store)

    if image_path is None:
        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
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
    seg_fam:str = None,
    plot: bool = False,
    image_path: Path = None,
    zarr_store: str | Path | None = None,
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
            seg_fam=seg_fam,
            cutoffs_path=cutoffs_path,
            plot=plot,
            image_path=image_path,
            zarr_store=zarr_store
        )


def write_adata(
    exp_name,
    reg_name: str = None,
    prefix_names: list[str] = None,
    output_path: Path = "/ceph/cephatlas/aklein/bican/data/anndatas/",
    zarr_store: str | Path | None = None,
):
    """ """
    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")

    logger.info(f"{exp_name}, {reg_name}, {prefix_names}, {output_path}")
    # Iterating over all prefixes
    for p in prefix_names:
        logger.info("PREFIX: %s" % p)
        # if a specific region is provided, write only that region
        if reg_name:
            fout = f"{output_path}/{exp_name}/{p}"
            _write_adata(exp_name, reg_name, p, fout)
        # if no region is provided, write all regions
        else:
            region_list = glob.glob(f"{zarr_store}/{exp_name}/region_*")
            for reg in region_list:
                rname = reg.split("/")[-1]
                fout = f"{output_path}/{exp_name}/{p}"
                _write_adata(exp_name, rname, p, fout, zarr_store)


def setup_adata_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str = "",
    plot=False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None,
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
        exp_name, reg_name, prefix_name, zarr_store=zarr_store
    )  # Use the filtered data if available
    # Run the setup
    adata = run_setup(adata, exp_name, reg_name, prefix_name, donor_name)
    # backup the adata object
    try:
        _assign_new_table(
            exp_name, reg_name, adata, KEYS[TABLE_KEY], suffix=suffix, zarr_store=zarr_store
        )  # double the storage but allows for iteration on filts.
    except Exception as e:
        logger.warning("Failed to assign new table, using backup instead")
        logger.warning(f"Error: {e}", exc_info=True)
        _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}", zarr_store=zarr_store)

    logger.info("DONE SETUP")
    if plot:
        plot_setup_region(
            exp_name, reg_name, prefix_name, image_path=image_path, image_store=image_store, suffix=suffix, zarr_store=zarr_store
        )


def setup_adata_all(
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None,
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
    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
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
            image_store=image_store,
            zarr_store=zarr_store,
        )


def plot_setup_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    suffix: str = "",
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None
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

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix, zarr_store=zarr_store)

    if image_path is None:
        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
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
    zarr_store: str | Path | None = None,
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
        exp_name, reg_name, prefix_name, suffix=suffix, zarr_store=zarr_store
    )  # Use the filtered data if available

    # Identify doublets
    adata = identify_doublets(adata, threshold=threshold)
    logger.info(
        f"Doublets identified: {adata.obs['doublet_bool'].sum()} out of {adata.n_obs} total cells"
    )

    # backup the adata object
    _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}", zarr_store=zarr_store)

    logger.info("DONE REMOVING DOUBLETS")

    if plot:  # There has to be a better way to do this!
        logger.info("PLOTTING DOUBLETS")
        # Importing here for efficiency
        from matplotlib.backends.backend_pdf import PdfPages
        from spida.pl import plot_doublets
        import spatialdata as sd
        from spida.pl import plot_example_doublets

        KEYS = _gen_keys(prefix_name, exp_name, reg_name)

        if zarr_store is None:
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
    zarr_store: str | Path | None = None,
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
    if zarr_store is None:
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
            zarr_store=zarr_store,
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
    max_epochs: int = 100,
    layer: str = "raw",
    batch_key: str = None,
    categorical_covariates: list = None,
    zarr_store: str | Path | None = None,
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
        exp_name, reg_name, prefix_name, suffix=suffix, zarr_store=zarr_store
    )  # Use the filtered data if available

    logger.info("arguments passing into resolvi")
    logger.info(f"max_epochs: {max_epochs}")
    logger.info(f"layer: {layer}")
    logger.info(f"batch_key: {batch_key}")
    logger.info(f"categorical_covariates: {categorical_covariates}")
    logger.info(f"model_save_path: {model_save_path}")
    logger.info(f"model_save_ext: {f'{exp_name}_{reg_name}_{prefix_name}'}")
    
    # Perform RESOLVI clustering
    adata = resolvi_cluster(
        adata,
        layer=layer,
        batch_key=batch_key,
        max_epochs=max_epochs,
        categorical_covariates=categorical_covariates,
        model_save_path=model_save_path,
        model_save_ext=f"{exp_name}_{reg_name}_{prefix_name}",
        trained=trained,
        **model_kwargs,
    )

    logger.info("DONE WITH RESOLVI, NOW CLUSTERING")

    # _calc_embeddings(adata, key_added="base_")
    logger.info("embedding resolvi")
    _calc_embeddings(adata, use_rep="X_resolvi", key_added="resolvi_")
    logger.info("embedding generated expression")
    _calc_embeddings(adata, layer="generated_expression", key_added="corr_")

    # backup the adata object
    _backup_adata(exp_name, reg_name, adata, f"{KEYS[TABLE_KEY]}{suffix}", zarr_store=zarr_store)

    logger.info("DONE RESOLVI CLUSTERING")
    if plot:
        plot_resolvi_region(
            exp_name, reg_name, prefix_name, image_path=image_path, suffix="_filt", zarr_store=zarr_store
        )


def resolvi_cluster_all(
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    model_save_path : Path = None, 
    zarr_store: str | Path | None = None,
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
    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
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
            zarr_store=zarr_store,
            **model_kwargs,
        )


def plot_resolvi_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
    zarr_store: str | Path | None = None,
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

    adata = _get_adata(exp_name, reg_name, prefix_name, suffix=suffix, zarr_store=zarr_store)

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
    zarr_store: Path = None,
    anndata_store: Path = None,
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
        zarr_store (Path): Path to the Zarr store containing the experiments
        anndata_store (Path): Path to save the aggregated AnnData object
        table_keys (list[str]): List of table keys to concatenate from the SpatialData object
    """
    from spida.P.project_level import aggregate_experiments, concatenate_tables

    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = Path(zarr_store)

    if anndata_store is None:
        anndata_store = Path(os.getenv("ANNDATA_STORE_PATH"))
    output_path = Path(anndata_store) / f"{project_name}_{lab_name}_{prefix_name}{suffix}.h5ad"

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
    anndata_store : Path = None,
): 
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    from spida.P.setup_adata import combined_setup
    import anndata as ad

    logger.info(f"SETUP ON {dataset_name}")
    if anndata_store is None: 
        anndata_store = Path(os.getenv("ANNDATA_STORE_PATH"))
    anndata_path = Path(anndata_store) / dataset_name
    try: 
        adata = ad.read_h5ad(anndata_path.with_suffix(".h5ad"))
    except FileNotFoundError:
        logger.warning(f"Anndata file not found at {anndata_path}. Please check the path.")
        logger.info("trying to read from file with no suffix")
        adata = ad.read_h5ad(anndata_path)
    logger.info(f"Loaded AnnData object with {adata.n_obs} cells and {adata.n_vars} genes.")

    adata = combined_setup(adata, scale=scale)
    adata.write_h5ad(anndata_path.with_suffix(".h5ad"))
    logger.info(f"DONE SETUP ON {dataset_name}")

def resolvi_dataset(
    dataset_name:str,
    anndata_store:Path = None,
    model_save_path:Path = None, 
    trained:bool = False,
    image_path:Path = None,
    max_epochs:int = 100,
    layer:str = "counts", 
    batch_key:str = "dataset_id",
    categorical_covariates:list = ['donor', 'brain_region'],
    **model_kwargs
    ): 
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    from spida.P.scvi_toolkit import resolvi_cluster
    from spida.P.setup_adata import _calc_embeddings
    import anndata as ad

    model_kwargs = model_kwargs['model_kwargs']
    
    logger.info(f"RESOLVI ON {dataset_name}")
    if anndata_store is None: 
        anndata_store = Path(os.getenv("ANNDATA_STORE_PATH"))
    anndata_path = anndata_store / dataset_name
    try: 
        adata = ad.read_h5ad(anndata_path.with_suffix(".h5ad"))
    except FileNotFoundError:
        logger.warning(f"Anndata file not found at {anndata_path}. Please check the path.")
        logger.info("trying to read from file with no suffix")
        adata = ad.read_h5ad(anndata_path)
    logger.info(f"Loaded AnnData object with {adata.n_obs} cells and {adata.n_vars} genes.")
    
    
    adata = resolvi_cluster(
        adata,
        layer=layer,
        batch_key=batch_key,
        categorical_covariates=categorical_covariates,
        max_epochs=max_epochs,
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

    adata.write_h5ad(anndata_path.with_suffix(".h5ad"))
    logger.info(f"DONE RESOLVI ON {dataset_name}")

    logger.info("PLOTTING DATASET")
    plot_dataset_setup(dataset_name, anndata_store=anndata_store, image_path=image_path, show=False)
    
    return 0

def plot_dataset_setup(
    dataset_name : str,
    anndata_store: str | Path = None,
    image_path: str | Path | None = None,
    show: bool = False
):
    """
    Plot the dataset results for a specific region in an experiment.
    """
    from spida.pl import plot_dataset

    adata = _read_adata(dataset_name, anndata_store=anndata_store)

    if image_path is None:
        image_store = os.getenv(
            "IMAGE_STORE_PATH",
        )
        image_path = Path(
            f"{image_store}/{dataset_name}"
        )
        image_path.mkdir(parents=True, exist_ok=True)

    plot_dataset(adata, save_path=image_path, show=show)


def _read_adata(
    dataset_name: str | Path,
    anndata_store: str | Path = None
): 
    """ handle adata reading with different suffixes and paths """
    import anndata as ad

    if anndata_store is None: 
        anndata_store = Path(os.getenv("ANNDATA_STORE_PATH"))
    if isinstance(anndata_store, str):
        anndata_store = Path(anndata_store)
    path = anndata_store / dataset_name
    try: 
        adata = ad.read_h5ad(path.with_suffix(".h5ad"))
    except FileNotFoundError:
        logger.warning(f"Anndata file not found at {path}. Please check the path.")
        logger.info("trying to read from file with no suffix")
        adata = ad.read_h5ad(path)
    logger.info(f"Loaded AnnData object with {adata.n_obs} cells and {adata.n_vars} genes.")
    return adata

def call_region_tz(
    exp_name : str,
    reg_name : str,
    prefix_name : str = "default",
    geoms_name : str = "wm_region",
    use_genes : str | list[str] = "BCAS1",
    save_geoms_path : str | None = None,
    dsc_comp_min_size : int = 5,
    hex_size : int = 50,
    hex_overlap : int = 0,
    gmm_ncomp : int | str = "auto",
    gmm_cov_type : str = "full",
    gene_agreement_thr : float = 0.75,
    top_n_comp : int = 1,
    gen_plots : bool = False,
    image_path: Path = None,
    image_store: Path = None,
    zarr_store: str | Path | None = None
): 
    """
    Call the call_regions function from the SPIDA.P module for a specific region in an experiment.

    Specify which genes to define the target regions using use_genes. 
    """
    from spida.P.call_regions import call_regions 
    
    # Get the zarr path
    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = Path(f"{zarr_store}/{exp_name}/{reg_name}")

    # Get the KEYS
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    # Get the image path
    if image_path is None:
        if image_store is None:
            image_store = os.getenv("IMAGE_STORE_PATH")
        image_path = Path(
            f"{image_store}/{exp_name}/{prefix_name}/{reg_name}"
        )
        image_path.mkdir(parents=True, exist_ok=True)

    call_regions(
        zarr_path = zarr_path,
        geoms_name = geoms_name,
        points_key = KEYS[POINTS_KEY],
        transfer_genes = use_genes,
        save_geoms_path = save_geoms_path,
        dsc_comp_min_size = dsc_comp_min_size,
        hex_size = hex_size,
        hex_overlap = hex_overlap,
        gmm_ncomp = gmm_ncomp,
        gmm_cov_type = gmm_cov_type,
        gene_agreement_thr = gene_agreement_thr,
        top_n_comp = top_n_comp,
        gen_plots = gen_plots,
        plot_save_path = image_path,
        image_key = KEYS[IMAGE_KEY],
    )

def transcript_qc_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str = "default",
    qc_shapes_key: str = "transcript_qc_shapes",
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int | None = 100,
    min_density: float | None = None,
    plot: bool = False,
    image_path: Path | None = None,
    image_store: Path | None = None,
    zarr_store: str | Path | None = None,
):
    import spatialdata as sd
    from spida.P.transcript_qc import run_transcript_qc
    from spida.pl import plot_hex_qc

    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = Path(f"{zarr_store}/{exp_name}/{reg_name}")

    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    sdata = sd.read_zarr(zarr_path)

    adata_qc, grid, _ = run_transcript_qc(
        sdata,
        points_key=KEYS[POINTS_KEY],
        qc_shapes_key=qc_shapes_key,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        min_transcripts=min_transcripts,
        min_density=min_density,
    )

    if plot:
        logger.info("PLOTTING TRANSCRIPT QC")
        if image_path is None:
            if image_store is None:
                image_store = os.getenv("IMAGE_STORE_PATH")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/transcript_qc.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)
        fig = plot_hex_qc(grid)
        fig.savefig(image_path, bbox_inches="tight")

def cluster_hexes_region(
    exp_name: str,
    reg_name: str,
    prefix_name: str = "default",
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int | None = 100,
    min_density: float | None = None,
    leiden_resolution: float = 1.0,
    min_cells: int = 10,
    min_genes: int = 100,
    n_top_genes: int = 200,
    pca_comps: int = 50,
    plot: bool = False,
    image_path: Path | None = None,
    image_store: Path | None = None,
    zarr_store: str | Path | None = None,
):
    import spatialdata as sd
    import geopandas as gpd
    from spida.P.transcript_qc import run_cluster_hexes, _obs_to_grid_geodf
    from spida.pl import plot_hex_clusters, save_cluster_panels_pdf

    if zarr_store is None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = Path(f"{zarr_store}/{exp_name}/{reg_name}")

    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    sdata = sd.read_zarr(zarr_path)

    adata_clustered, _ = run_cluster_hexes(
        sdata,
        points_key=KEYS[POINTS_KEY],
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
    )

    if plot:
        logger.info("PLOTTING CLUSTER HEXES")
        if image_path is None:
            if image_store is None:
                image_store = os.getenv("IMAGE_STORE_PATH")
            image_path = Path(f"{image_store}/{exp_name}/{prefix_name}/{reg_name}/cluster_hexes.pdf")
            image_path.parent.mkdir(parents=True, exist_ok=True)
        grid = _obs_to_grid_geodf(adata_clustered.obs.copy())
        figs = plot_hex_clusters(grid, adata_clustered)
        save_cluster_panels_pdf(image_path, figs)