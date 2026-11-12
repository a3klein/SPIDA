from pathlib import Path
import warnings
import logging

from spida.config import SpidaConfig
from spida.io.layout import Dataset, TreeLayout, dataset_from, table_key
from spida.io.store import is_writeable, obs_columns, open_store, read_table, require, write_table
from spida.utilities.sd_utils import _write_adata
from spida._constants import TABLE_KEY, IMAGE_KEY, SHAPES_KEY, POINTS_KEY  # type: ignore

logger = logging.getLogger(__package__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


def _require_dir(path, key: str) -> Path:
    if path is None:
        raise ValueError(f"no directory given; pass it or set {key} in the --config JSON or the environment")
    return Path(path)


def _figure(dataset: Dataset, label: str, name: str, image_path=None) -> Path:
    """Explicit image_path, else ``{figures}/{label}/{name}``; parent dir created."""
    path = Path(image_path) if image_path is not None else dataset.figures_dir(label) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


_FILTER_COLUMNS = ["nCount_RNA", "nFeature_RNA", "nBlank", "nCount_RNA_per_Volume", "pass_qc_pre", "pass_qc"]
_DOUBLET_COLUMNS = ["doublet_score", "doublet_bool"]
_RESOLVI_COLUMNS = ["true_proportion", "diffusion_proportion", "background_proportion", "resolvi_leiden", "corr_leiden"]


def _skip_if_written(reg_name, what, **check) -> bool:
    # -all loops: True (with a warning) when this region's outputs already exist
    try:
        is_writeable(**check)
        return False
    except FileExistsError:
        logger.warning("Skipping %s: %s already done (pass --overwrite to redo)", reg_name, what)
        return True


def filter_cells_region(
    dataset: Dataset,
    prefix_name: str,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
):
    """
    Filter the cells of one segmentation table.
    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    overwrite (bool, optional): Re-filter a table that was already filtered. Defaults to False.
    """
    from spida.P.filtering import run_filtering

    logger.info("FILTERING CELLS, %s, PREFIX %s" % (dataset.zarr, prefix_name))

    KEYS = dataset.keys(prefix_name)
    zarr_path = dataset.require("zarr")
    require(zarr_path, [KEYS[TABLE_KEY]])
    is_writeable(zarr_path, columns={KEYS[TABLE_KEY]: _FILTER_COLUMNS}, overwrite=overwrite)
    adata = read_table(zarr_path, KEYS[TABLE_KEY])
    # Run the filtering
    adata = run_filtering(
        adata=adata,
        meta=dataset.meta(prefix_name),
        cutoffs_path=cutoffs_path
    )
    # backup the AnnData object
    # _assign_new_table(zarr_path, adata, f"{KEYS[TABLE_KEY]}_filt") # double the storage but allows for iteration on filts.
    write_table(zarr_path, KEYS[TABLE_KEY], adata, overwrite=True)

    logger.info(
        f"Passed QC Cells: {adata.obs['pass_qc'].sum()} out of {adata.n_obs} total cells"
    )
    logger.info("DONE")
    if plot:
        plot_filtering_region(dataset, prefix_name, image_path=image_path)


def plot_filtering_region(
    dataset: Dataset,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
):
    """
    Plot the filtering results of one segmentation table.
    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot (default ``{figures}/{prefix}/pixi_filt.pdf``).
    """
    from spida.pl import plot_filtering
    from matplotlib.backends.backend_pdf import PdfPages

    adata = read_table(dataset.require("zarr"), f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}")
    image_path = _figure(dataset, prefix_name, "pixi_filt.pdf", image_path)

    pdf_file = PdfPages(image_path)
    fig, ax = plot_filtering(adata, dataset.meta(prefix_name))
    pdf_file.savefig(fig)
    pdf_file.close()


def filter_cells_all(
    config: SpidaConfig,
    exp_name: str,
    prefix_name: str,
    cutoffs_path: Path = None,
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
):
    """
    Filter cells for every region of an experiment in the zarr store.
    Regions already filtered are skipped unless ``overwrite``.
    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    cutoffs_path (Path, optional): Path to the cutoffs JSON file.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    overwrite (bool, optional): Re-filter regions already filtered. Defaults to False.
    """
    for reg_name in TreeLayout(config).units(exp_name):
        dataset = dataset_from(config, exp_name, reg_name)
        if _skip_if_written(reg_name, f"{prefix_name} filtering", zarr_path=dataset.zarr, overwrite=overwrite,
                            columns={dataset.keys(prefix_name)[TABLE_KEY]: _FILTER_COLUMNS}):
            continue
        filter_cells_region(
            dataset,
            prefix_name,
            cutoffs_path=cutoffs_path,
            plot=plot,
            image_path=image_path,
            overwrite=overwrite,
        )


def write_adata(
    config: SpidaConfig,
    exp_name,
    reg_name: str = None,
    prefix_names: list[str] = None,
    suffix: str = "",
    overwrite: bool = False,
):
    """
    Export tables ``{prefix}_table{suffix}`` to their h5ad copies, ``{anndata}/{table_key}.h5ad``.
    One region if ``reg_name`` is given, else every region of the experiment.
    An existing h5ad is replaced only with ``overwrite`` (one region: raises; all regions: skipped).
    """
    regions = [reg_name] if reg_name else TreeLayout(config).units(exp_name)
    logger.info(f"{exp_name}, {regions}, {prefix_names}")
    for p in prefix_names:
        logger.info("PREFIX: %s" % p)
        key = table_key(p, suffix)
        for reg in regions:
            dataset = dataset_from(config, exp_name, reg)
            out = dataset.table_h5ad(key)
            require(dataset.require("zarr"), [key])
            if reg_name:
                is_writeable(files=[out], overwrite=overwrite)
            elif _skip_if_written(reg, out.name, files=[out], overwrite=overwrite):
                continue
            _write_adata(dataset.require("zarr"), key, out)


def setup_adata_region(
    dataset: Dataset,
    prefix_name: str,
    suffix: str = "",
    plot=False,
    do_tsne: bool = False,
    image_path: Path = None,
    consensus_cluster: bool = True,
    n_jobs: int | None = None,
    overwrite: bool = False,
):
    """
    Setup the AnnData object for downstream analysis.
    This involves normalizing data, calculating PCA, umap, tsne, and leiden clusters

    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    do_tsne (bool, optional): Also compute the tSNE embedding (needed by the setup plots). Defaults to False.
    consensus_cluster (bool, optional): Use consensus clustering rather than a single
        Leiden run. Defaults to True. See `run_setup` for the memory/runtime tradeoff.
    n_jobs (int | None, optional): Consensus clustering workers. Defaults to None: the job's cores, max 24.
    overwrite (bool, optional): Rebuild ``{prefix}_table{suffix}`` if it exists (drops doublet / resolvi
        results on it). Defaults to False.
    """
    from spida.P.setup_adata import run_setup

    logger.info("SETTING UP ADATA, %s, PREFIX %s" % (dataset.zarr, prefix_name))

    KEYS = dataset.keys(prefix_name)
    zarr_path = dataset.require("zarr")
    out = f"{KEYS[TABLE_KEY]}{suffix}"
    require(zarr_path, columns={KEYS[TABLE_KEY]: ["pass_qc"]})
    is_writeable(zarr_path, [out], overwrite=overwrite)
    lost = [c for c in obs_columns(zarr_path, out) if c in _DOUBLET_COLUMNS + _RESOLVI_COLUMNS]
    if lost:
        logger.warning("Rebuilding %s drops these columns from later steps: %s", out, ", ".join(lost))
    adata = read_table(zarr_path, KEYS[TABLE_KEY])
    # Run the setup
    adata = run_setup(adata, dataset.meta(prefix_name), do_tsne=do_tsne,
                      consensus_cluster=consensus_cluster, n_jobs=n_jobs)
    # backup the adata object
    write_table(zarr_path, out, adata, overwrite=True)  # new table or swap if overwrite

    logger.info("DONE SETUP")
    if plot:
        plot_setup_region(dataset, prefix_name, image_path=image_path, suffix=suffix)


def setup_adata_all(
    config: SpidaConfig,
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
):
    """
    Setup the AnnData objects for every region of an experiment in the zarr store.
    Regions whose ``{prefix}_table{suffix}`` exists are skipped unless ``overwrite``.

    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    overwrite (bool, optional): Rebuild existing tables. Defaults to False.
    """
    for reg_name in TreeLayout(config).units(exp_name):
        dataset = dataset_from(config, exp_name, reg_name)
        out = f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}"
        if _skip_if_written(reg_name, out, zarr_path=dataset.zarr, elements=[out], overwrite=overwrite):
            continue
        setup_adata_region(
            dataset,
            prefix_name,
            suffix=suffix,
            plot=plot,
            image_path=image_path,
            overwrite=overwrite,
        )


def plot_setup_region(
    dataset: Dataset,
    prefix_name: str,
    suffix: str = "",
    image_path: Path = None,
):
    """
    Plot the setup results of one segmentation table.

    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot (default ``{figures}/{prefix}/pixi_setup.pdf``).
    """
    from spida.pl import plot_setup
    from matplotlib.backends.backend_pdf import PdfPages

    adata = read_table(dataset.require("zarr"), f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}")
    image_path = _figure(dataset, prefix_name, "pixi_setup.pdf", image_path)

    pdf_file = PdfPages(image_path)
    plot_setup(adata, dataset.meta(prefix_name), pdf_file=pdf_file)
    pdf_file.close()


def remove_doublets_region(
    dataset: Dataset,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
    **model_kwargs,
):
    """
    Flag doublets (SOLO) in one segmentation table.

    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    threshold (float, optional): Threshold for doublet detection. Defaults to 0.5.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    overwrite (bool, optional): Re-run on a table that already has doublet calls. Defaults to False.
    """
    from spida.P.scvi_toolkit import identify_doublets  # , remove_doublets

    logger.info("REMOVING DOUBLETS, %s, PREFIX %s" % (dataset.zarr, prefix_name))

    KEYS = dataset.keys(prefix_name)
    zarr_path = dataset.require("zarr")
    require(zarr_path, [f"{KEYS[TABLE_KEY]}{suffix}"])
    is_writeable(zarr_path, columns={f"{KEYS[TABLE_KEY]}{suffix}": _DOUBLET_COLUMNS}, overwrite=overwrite)
    adata = read_table(zarr_path, f"{KEYS[TABLE_KEY]}{suffix}")

    # Identify doublets
    adata = identify_doublets(adata, threshold=threshold)
    logger.info(
        f"Doublets identified: {adata.obs['doublet_bool'].sum()} out of {adata.n_obs} total cells"
    )

    # backup the adata object
    write_table(zarr_path, f"{KEYS[TABLE_KEY]}{suffix}", adata, overwrite=True)

    logger.info("DONE REMOVING DOUBLETS")

    if plot:  # There has to be a better way to do this!
        logger.info("PLOTTING DOUBLETS")
        # Importing here for efficiency
        from matplotlib.backends.backend_pdf import PdfPages
        from spida.pl import plot_doublets
        from spida.pl import plot_example_doublets

        sdata = open_store(zarr_path)

        image_path = _figure(dataset, prefix_name, "SOLO_doublets.png", image_path)

        plot_doublets(adata, image_path)
        image_path = image_path.parent / "SOLO_examples.pdf"
        with PdfPages(image_path) as pdf:
            plot_example_doublets(
                sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[TABLE_KEY], pdf_file=pdf
            )

    # Remove doublets
    # adata = remove_doublets(adata)


def remove_doublets_all(
    config: SpidaConfig,
    exp_name: str,
    prefix_name: str,
    threshold: float = 0.5,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    overwrite: bool = False,
    **model_kwargs,
):
    """
    Flag doublets for every region of an experiment in the zarr store.
    Regions that already have doublet calls are skipped unless ``overwrite``.

    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    threshold (float, optional): Threshold for doublet detection. Defaults to 0.5.
    plot (bool, optional): Whether to plot the results. Defaults to False.
    overwrite (bool, optional): Re-run regions that already have doublet calls. Defaults to False.
    """
    logger.info("REMOVING DOUBLETS FOR ALL REGIONS IN EXPERIMENT %s" % exp_name)
    for reg_name in TreeLayout(config).units(exp_name):
        dataset = dataset_from(config, exp_name, reg_name)
        key = f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}"
        if _skip_if_written(reg_name, f"doublets on {key}", zarr_path=dataset.zarr,
                            columns={key: _DOUBLET_COLUMNS}, overwrite=overwrite):
            continue
        remove_doublets_region(
            dataset,
            prefix_name,
            threshold=threshold,
            plot=plot,
            image_path=image_path,
            suffix=suffix,
            overwrite=overwrite,
            **model_kwargs,
        )


def resolvi_cluster_region(
    dataset: Dataset,
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
    overwrite: bool = False,
    **model_kwargs,
):
    """
    Perform RESOLVI clustering on one segmentation table.

    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    suffix (str, optional): Suffix for the AnnData object. Defaults to "".
    plot (bool, optional): Whether to plot the results. Defaults to False.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    overwrite (bool, optional): Re-run on a table that already has RESOLVI results. Defaults to False.
    """
    from spida.P.scvi_toolkit import resolvi_cluster
    from spida.utilities.ad_utils import _calc_embeddings

    logger.info("RESOLVI CLUSTERING, %s, PREFIX %s" % (dataset.zarr, prefix_name))

    model_kwargs = model_kwargs['model_kwargs']
    for key, val in model_kwargs.items():
        logger.info(f"Model kwargs: {key} = {val}")

    KEYS = dataset.keys(prefix_name)
    zarr_path = dataset.require("zarr")
    require(zarr_path, [f"{KEYS[TABLE_KEY]}{suffix}"])
    is_writeable(zarr_path, columns={f"{KEYS[TABLE_KEY]}{suffix}": _RESOLVI_COLUMNS}, overwrite=overwrite)
    adata = read_table(zarr_path, f"{KEYS[TABLE_KEY]}{suffix}")
    model_save_ext = "_".join(p for p in (dataset.namespace, prefix_name) if p)

    logger.info("arguments passing into resolvi")
    logger.info(f"max_epochs: {max_epochs}")
    logger.info(f"layer: {layer}")
    logger.info(f"batch_key: {batch_key}")
    logger.info(f"categorical_covariates: {categorical_covariates}")
    logger.info(f"model_save_path: {model_save_path}")
    logger.info(f"model_save_ext: {model_save_ext}")

    # Perform RESOLVI clustering
    adata = resolvi_cluster(
        adata,
        layer=layer,
        batch_key=batch_key,
        max_epochs=max_epochs,
        categorical_covariates=categorical_covariates,
        model_save_path=model_save_path,
        model_save_ext=model_save_ext,
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
    write_table(zarr_path, f"{KEYS[TABLE_KEY]}{suffix}", adata, overwrite=True)

    logger.info("DONE RESOLVI CLUSTERING")
    if plot:
        plot_resolvi_region(dataset, prefix_name, image_path=image_path, suffix="_filt")


def resolvi_cluster_all(
    config: SpidaConfig,
    exp_name: str,
    prefix_name: str,
    suffix: str = "",
    plot: bool = False,
    image_path: Path = None,
    model_save_path : Path = None,
    overwrite: bool = False,
    **model_kwargs,
):
    """
    Perform RESOLVI clustering for every region of an experiment in the zarr store.
    Regions that already have RESOLVI results are skipped unless ``overwrite``.
    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    suffix (str, optional): Suffix for the AnnData object. Defaults to "".
    plot (bool, optional): Whether to plot the results. Defaults to False.
    image_path (Path, optional): Path to save the plot. If None, uses a default path.
    overwrite (bool, optional): Re-run regions that already have RESOLVI results. Defaults to False.
    """
    for reg_name in TreeLayout(config).units(exp_name):
        dataset = dataset_from(config, exp_name, reg_name)
        key = f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}"
        if _skip_if_written(reg_name, f"RESOLVI on {key}", zarr_path=dataset.zarr,
                            columns={key: _RESOLVI_COLUMNS}, overwrite=overwrite):
            continue
        resolvi_cluster_region(
            dataset,
            prefix_name,
            suffix=suffix,
            plot=plot,
            image_path=image_path,
            model_save_path=model_save_path,
            overwrite=overwrite,
            **model_kwargs,
        )


def plot_resolvi_region(
    dataset: Dataset,
    prefix_name: str,
    image_path: Path = None,
    suffix: str = "",
):
    """
    Plot the RESOLVI results of one segmentation table.

    Parameters:
    dataset (Dataset): The region's locations (zarr store, figures dir) and names.
    prefix_name (str): Prefix for the keys in the spatialdata object.
    image_path (Path, optional): Path to save the plot (default ``{figures}/{prefix}/pixi_resolvi.pdf``).
    """
    from spida.pl import plot_resolvi
    from matplotlib.backends.backend_pdf import PdfPages

    adata = read_table(dataset.require("zarr"), f"{dataset.keys(prefix_name)[TABLE_KEY]}{suffix}")
    image_path = _figure(dataset, prefix_name, "pixi_resolvi.pdf", image_path)

    pdf_file = PdfPages(image_path)
    plot_resolvi(adata, dataset.meta(prefix_name), pdf_file=pdf_file)
    pdf_file.close()


def combine_datasets(
    experiment_list: list[str],
    prefix_name: str,
    suffix: str = "_filt",
    lab_name: str = "salk",
    project_name: str = "BG_salk",
    config: SpidaConfig | None = None,
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
        config (SpidaConfig): Store roots (zarr_store is read, anndata_store receives the h5ad)
        table_keys (list[str]): List of table keys to concatenate from the SpatialData object
    """
    from spida.P.project_level import aggregate_experiments, concatenate_tables

    config = config or SpidaConfig()
    zarr_path = config.require("zarr_store")
    output_path = config.require("anndata_store") / f"{project_name}_{lab_name}_{prefix_name}{suffix}.h5ad"

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
    anndata_path = _require_dir(anndata_store, "ANNDATA_STORE_PATH") / dataset_name
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
    image_store: Path | None = None,
    **model_kwargs
    ):
    """
    This function runs resolvi on a dataset level object. Takes more time and memory than resolvi_region.
    """
    from spida.P.scvi_toolkit import resolvi_cluster
    from spida.utilities.ad_utils import _calc_embeddings
    import anndata as ad

    model_kwargs = model_kwargs['model_kwargs']
    
    logger.info(f"RESOLVI ON {dataset_name}")
    anndata_store = _require_dir(anndata_store, "ANNDATA_STORE_PATH")
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
    plot_dataset_setup(dataset_name, anndata_store=anndata_store, image_path=image_path,
                       image_store=image_store, show=False)
    
    return 0

def plot_dataset_setup(
    dataset_name : str,
    anndata_store: str | Path = None,
    image_path: str | Path | None = None,
    image_store: str | Path | None = None,
    show: bool = False
):
    """
    Plot a dataset-level h5ad; figures go to ``image_path`` or ``{image_store}/{dataset_name}``.
    """
    from spida.pl import plot_dataset

    adata = _read_adata(dataset_name, anndata_store=anndata_store)

    if image_path is None:
        image_path = _require_dir(image_store, "IMAGE_STORE_PATH") / dataset_name
        image_path.mkdir(parents=True, exist_ok=True)

    plot_dataset(adata, save_path=image_path, show=show)


def _read_adata(
    dataset_name: str | Path,
    anndata_store: str | Path = None
): 
    """ handle adata reading with different suffixes and paths """
    import anndata as ad

    path = _require_dir(anndata_store, "ANNDATA_STORE_PATH") / dataset_name
    try: 
        adata = ad.read_h5ad(path.with_suffix(".h5ad"))
    except FileNotFoundError:
        logger.warning(f"Anndata file not found at {path}. Please check the path.")
        logger.info("trying to read from file with no suffix")
        adata = ad.read_h5ad(path)
    logger.info(f"Loaded AnnData object with {adata.n_obs} cells and {adata.n_vars} genes.")
    return adata

def call_region_tz(
    dataset: Dataset,
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
    overwrite: bool = False,
):
    """
    Call spatial regions (e.g. white matter) from transcript density of ``use_genes``.

    ``image_path`` is the plot directory (default ``{figures}/{prefix}``).
    Existing ``geoms_name`` shapes (and ``save_geoms_path``) are replaced only with ``overwrite``.
    """
    from spida.P.call_regions import call_regions

    KEYS = dataset.keys(prefix_name)
    require(dataset.require("zarr"), [KEYS[POINTS_KEY]] + ([KEYS[IMAGE_KEY]] if gen_plots else []))
    is_writeable(dataset.require("zarr"), [geoms_name],
                 files=[save_geoms_path] if save_geoms_path else [], overwrite=overwrite)
    image_path = Path(image_path) if image_path is not None else dataset.figures_dir(prefix_name)
    image_path.mkdir(parents=True, exist_ok=True)

    call_regions(
        zarr_path = dataset.require("zarr"),
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
        overwrite = overwrite,
    )

def transcript_qc_region(
    dataset: Dataset,
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
    overwrite: bool = False,
):
    from spida.P.transcript_qc import run_transcript_qc
    from spida.pl import plot_hex_qc

    KEYS = dataset.keys(prefix_name)
    require(dataset.require("zarr"), [KEYS[POINTS_KEY]])
    is_writeable(dataset.require("zarr"), [qc_shapes_key], overwrite=overwrite)  # hex table is shared, reused
    sdata = open_store(dataset.require("zarr"))

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
        overwrite=overwrite,
    )

    if plot:
        logger.info("PLOTTING TRANSCRIPT QC")
        image_path = _figure(dataset, prefix_name, "transcript_qc.pdf", image_path)
        fig = plot_hex_qc(grid)
        fig.savefig(image_path, bbox_inches="tight")

def tissue_mask_region(
    dataset: Dataset,
    prefix_name: str = "default",
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    rho_min: float = 5.0,
    ref_um: float = 25.0,
    erode_sigma: float = 1.0,
    a_min_um2: float = 50000.0,
    gb: float = 5.0,
    write_raster: bool = True,
    shapes_key: str = "qc_mask",
    labels_key: str = "tissue_mask",
    plot: bool = False,
    image_path: Path | None = None,
    mask_dir: Path | None = None,
    overwrite: bool = False,
):
    """Compute the KDE tissue mask; persist polygon + labels, and the tif to ``mask_dir`` (default: raw images dir).
    Existing polygon / labels / tif are replaced only with ``overwrite``."""
    from spida.P.tissue_mask import run_tissue_mask, QC_MASK_TIF_NAME
    from spida.io.formats.merscope import TRANSFORM_FILE
    from spida.io.store import is_writeable, require

    images_dir = dataset.raw_images()
    KEYS = dataset.keys(prefix_name)
    tif = Path(mask_dir or images_dir) / QC_MASK_TIF_NAME
    require(dataset.require("zarr"), [KEYS[POINTS_KEY]],
            files=[images_dir / TRANSFORM_FILE, images_dir] if write_raster else [])
    is_writeable(dataset.require("zarr"), [shapes_key, labels_key],
                 files=[tif] if write_raster else [], overwrite=overwrite)
    sdata = open_store(dataset.require("zarr"))

    stages, grid, x0, y0 = run_tissue_mask(
        sdata,
        points_key=KEYS[POINTS_KEY],
        images_dir=images_dir,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        gb=gb,
        rho_min=rho_min,
        ref_um=ref_um,
        a_min_um2=a_min_um2,
        erode_sigma=erode_sigma,
        shapes_key=shapes_key,
        labels_key=labels_key,
        write_raster=write_raster,
        mask_dir=mask_dir,
        overwrite=overwrite,
    )

    if plot:
        logger.info("PLOTTING TISSUE MASK")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        image_path = _figure(dataset, prefix_name, "tissue_mask.png", image_path)
        fig, ax = plt.subplots(1, 3, figsize=(18, 6))
        ax[0].imshow(np.log1p(stages["dens100"]), origin="lower", cmap="magma", aspect="equal")
        ax[0].set_title(f"log transcript density\n{dataset.exp or ''}\n{dataset.reg or dataset.zarr.name}", fontsize=9)
        ax[1].imshow(stages["final"], origin="lower", cmap="gray_r", aspect="equal")
        ax[1].set_title(f"tissue mask (rho>={rho_min:g} erode={erode_sigma:g}sig)\n"
                        f"{int(stages['final'].sum())} px ({100*stages['final'].mean():.1f}%)", fontsize=9)
        ax[2].imshow(np.log1p(stages["dens100"]), origin="lower", cmap="magma", aspect="equal")
        for geom in grid.geometry:
            for poly in (geom.geoms if geom.geom_type == "MultiPolygon" else [geom]):
                xs, ys = poly.exterior.xy
                ax[2].plot((np.asarray(xs) - x0) / gb, (np.asarray(ys) - y0) / gb, color="lime", lw=0.6)
        ax[2].set_title("qc_mask polygon over density", fontsize=9)
        for a in ax:
            a.set_xticks([]); a.set_yticks([])
        fig.tight_layout()
        fig.savefig(image_path, dpi=120, bbox_inches="tight")
        logger.info("Saved tissue mask figure to %s", image_path)

def cluster_hexes_region(
    dataset: Dataset,
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
    overwrite: bool = False,
):
    from spida.P.transcript_qc import run_cluster_hexes, _obs_to_grid_geodf, _hex_table_key
    from spida.pl import plot_hex_clusters, save_cluster_panels_pdf

    KEYS = dataset.keys(prefix_name)
    require(dataset.require("zarr"), [KEYS[POINTS_KEY]])
    is_writeable(dataset.require("zarr"), [f"{_hex_table_key(hex_size, hex_overlap)}_clustered"],
                 overwrite=overwrite)  # hex table is shared, reused
    sdata = open_store(dataset.require("zarr"))

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
        overwrite=overwrite,
    )

    if plot:
        logger.info("PLOTTING CLUSTER HEXES")
        image_path = _figure(dataset, prefix_name, "cluster_hexes.pdf", image_path)
        grid = _obs_to_grid_geodf(adata_clustered.obs.copy())
        figs = plot_hex_clusters(grid, adata_clustered)
        save_cluster_panels_pdf(image_path, figs)