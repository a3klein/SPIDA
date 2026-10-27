from spida._constants import *
from spida.config import SpidaConfig
from spida.io.layout import Dataset, TreeLayout, dataset_from, element_keys
from spida.io.store import is_writeable, open_store, require
from spida.io.load import (
    read_merscope,
    read_segmentation,
    write_segmentation,
    load_decon_images,
)
from .seg_qc_filter import apply_qc_filter_to_segmentation
from spida.pl import plot_images, plot_shapes, plot_points, plot_overlap, plot_seg_load

from pathlib import Path
import warnings
import logging

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

from matplotlib.backends.backend_pdf import PdfPages

logger = logging.getLogger(__package__)

# load-segmentation options passed to the segmentation QC filter, not the reader
_QC_FILTER_KWARGS = ("qc_regions", "qc_table_key", "qc_shapes_key", "qc_filter_col", "qc_pass_value",
                     "table_cell_id_col", "points_cell_id_col")


def ingest_region(
    dataset: Dataset,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    overwrite: bool = False,
    **kwargs,
):
    """
    Ingest one MERSCOPE region into a new zarr store.

    Parameters:
    dataset (Dataset): Where the region's raw data, zarr store and figures live, and its names.
    type (str): Type of the data to ingest (default is "merscope").
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    source (str): Source of the data (default is "machine").
    plot (bool): Whether to generate plots (default is False).
    overwrite (bool): Replace an existing zarr store for this region (default is False).
    """
    namespace = f"{dataset.exp}_{dataset.reg}" if dataset.named else ""
    logger.info("INGESTING REGION %s -> %s (element names suffix %r)", dataset.raw, dataset.zarr, namespace)
    KEYS = element_keys(prefix_name, namespace)
    require(files=[dataset.require("raw")])
    is_writeable(files=[dataset.require("zarr")], overwrite=overwrite)

    if type == "merscope" and source == "machine":
        sdata = read_merscope(
            dataset.require("raw"),
            dataset.require("zarr"),
            keys=KEYS,
            namespace=namespace,
            exp_name=dataset.exp,
            reg_name=dataset.reg,
            prefix_name=prefix_name,
        )

        image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

        logger.info(sdata.tables.keys())
        if plot:
            image_path = dataset.figures_dir("default") / "pixi-ing.pdf"
            image_path.parent.mkdir(parents=True, exist_ok=True)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                with PdfPages(image_path) as pdf:
                    plot_images(
                        sdata,
                        KEYS[IMAGE_KEY],
                        image_scale_keys,
                        image_channels,
                        cs="global",
                        pdf_file=pdf,
                    )
                    plot_shapes(
                        sdata,
                        KEYS[SHAPES_KEY],
                        table_name=KEYS[TABLE_KEY],
                        cs="pixel",
                        pdf_file=pdf,
                    )
                    plot_points(
                        sdata,
                        KEYS[POINTS_KEY],
                        KEYS[TABLE_KEY],
                        cs="pixel",
                        cmap="tab10",
                        pdf_file=pdf,
                    )
                    plot_overlap(
                        sdata,
                        KEYS[IMAGE_KEY],
                        KEYS[SHAPES_KEY],
                        KEYS[POINTS_KEY],
                        KEYS[TABLE_KEY],
                        image_scale_keys,
                        cs="pixel",
                        pdf_file=pdf,
                    )


def ingest_all(
    config: SpidaConfig,
    exp_name: str,
    type: str = "merscope",
    prefix_name: str = "default",
    source: str = "machine",
    plot: bool = False,
    overwrite: bool = False,
    **kwargs,
):
    """
    Ingest every region of an experiment found under the raw root.
    Regions that already have a zarr store are skipped unless ``overwrite``.

    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    overwrite (bool): Replace existing zarr stores (default is False).
    """
    logger.info("INGESTING ALL REGIONS; EXPERIMENT %s" % (exp_name))
    if type == "merscope" and source == "machine":
        for reg in TreeLayout(config).units(exp_name, under="raw"):
            dataset = dataset_from(config, exp_name, reg)
            try:
                is_writeable(files=[dataset.require("zarr")], overwrite=overwrite)
            except FileExistsError:
                logger.warning("Skipping %s: zarr store %s already exists (pass --overwrite to replace)", reg, dataset.zarr)
                continue
            ingest_region(
                dataset,
                type=type,
                prefix_name=prefix_name,
                source=source,
                plot=plot,
                overwrite=overwrite,
            )


def _load_inputs(
    spec,
    seg_region: Path,
    load_kwargs: dict,
) -> list[Path]:
    """
    Files load-segmentation reads from a segmentation's folder (names overridable via ``load_kwargs``).
    Parameters:
    spec (SegmentationClass): The method.
    seg_region (Path): The segmentation's folder.
    load_kwargs (dict): The loader's keyword arguments.
    """
    from spida.io.formats.segmenters import PROSEG_POLYGONS
    from spida.io.schema.segmentation import (
        SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_CELL_METADATA, SCHEMA_TRANSCRIPTS, resolve,
    )

    meta = load_kwargs.get("cell_metadata_fname", SCHEMA_CELL_METADATA)
    if spec.provides_assignment:
        if spec.spatialdata_file:                                    # proseg v3: the zarr (metadata optional)
            return [seg_region / spec.spatialdata_file]
        return [seg_region / spec.transcripts_file, seg_region / spec.counts_file,
                seg_region / PROSEG_POLYGONS,
                resolve(seg_region, meta, "cell_metadata_with_signals.csv")]
    boundaries = load_kwargs.get("boundaries_fname") or SCHEMA_BOUNDARIES
    return [resolve(seg_region, boundaries, spec.legacy_boundaries_file),
            seg_region / load_kwargs.get("cell_by_gene_fname", SCHEMA_CELL_BY_GENE),
            seg_region / meta,
            seg_region / load_kwargs.get("detected_transcripts_fname", SCHEMA_TRANSCRIPTS)]


def load_segmentation_region(
    dataset: Dataset,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    overwrite: bool = False,
    **load_kwargs,
):
    """
    Load one segmentation's schema files into the dataset's store.

    The segmentation directory is ``dataset.seg_dir(label)`` (label = ``prefix_name``, or the
    method for "default"); the legacy ``{exp}/{label}/{reg}`` folder structure is read as a fallback.

    Parameters:
    dataset (Dataset): The region's locations (zarr store, segmentation dir, figures dir) and names.
    type (str): Type of the segmentation data to load (default is "vpt").
    prefix_name (str): Segmentation name; used both as the storage label and the spatialdata
        key prefix (default is "default").
    overwrite (bool): Replace this segmentation's points, shapes and table if already loaded (default is False).
    """

    load_kwargs = load_kwargs.get("load_kwargs", load_kwargs)
    qc_kwargs = {k: load_kwargs.pop(k) for k in _QC_FILTER_KWARGS if k in load_kwargs}
    zarr_path = dataset.require("zarr")
    logger.info("LOADING SEGMENTATION %s (%s) into %s", prefix_name, type, zarr_path)

    # Map the legacy `type` to a segmentation method; the spec then drives the
    # read path (boundary-only segmentation schema vs proseg-native) and legacy filename fallback.
    from spida.S.segmentation.backends import get_spec

    _TYPE_TO_METHOD = {
        "vpt": "cellpose",        # legacy: VPT-format segmentation schema of a boundary-only segmenter
        "cellpose": "cellpose",
        "mesmer": "mesmer",
        "proseg": "proseg",
    }
    method = _TYPE_TO_METHOD.get(type, type)
    version = load_kwargs.pop("version", None)    # proseg v2/v3 (defaults to v3)
    spec = get_spec(method, version)

    # Label = prefix_name (the method for the "default" sentinel).
    label = prefix_name if prefix_name and prefix_name != "default" else method
    seg_region = dataset.seg_dir(label, must_exist=True)
    KEYS = dataset.keys(prefix_name)
    qc_elements = [qc_kwargs[k] for k in ("qc_shapes_key", "qc_table_key") if qc_kwargs.get(k) is not None]
    require(zarr_path, [dataset.keys("default")[SHAPES_KEY], *qc_elements],
            files=_load_inputs(spec, seg_region, load_kwargs))
    is_writeable(zarr_path, [KEYS[POINTS_KEY], KEYS[SHAPES_KEY], KEYS[TABLE_KEY]], overwrite=overwrite)
    sdata = open_store(zarr_path)

    sdata = read_segmentation(
        sdata,
        spec,
        keys=KEYS,
        default_keys=dataset.keys("default"),
        seg_dir=str(seg_region),
        exp_name=dataset.exp,
        reg_name=dataset.reg,
        dataset_id=dataset.namespace,
        prefix_name=prefix_name,
        **load_kwargs,
    )
    if any(qc_kwargs.get(k) is not None for k in ("qc_regions", "qc_table_key", "qc_shapes_key")):
        logger.info("Applying transcript QC filter to loaded segmentation elements.")
        sdata = apply_qc_filter_to_segmentation(
            sdata,
            table_key=KEYS[TABLE_KEY],
            shapes_key=KEYS[SHAPES_KEY],
            points_key=KEYS[POINTS_KEY],
            **qc_kwargs,
        )
    sdata = write_segmentation(sdata, KEYS, overwrite=overwrite)

    # TODO: add a function to filter the spatialdata object to something very permissive (like at least 5 transcripts per cell?)

    if plot:
        # image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
        image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

        image_path = dataset.figures_dir(prefix_name) / "pixi-load.pdf"
        image_path.parent.mkdir(parents=True, exist_ok=True)

        # plot
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with PdfPages(image_path) as pdf:
                # plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                # plot_shapes(
                #     sdata,
                #     KEYS[SHAPES_KEY],
                #     table_name=KEYS[TABLE_KEY],
                #     cs="pixel",
                #     pdf_file=pdf,
                # )
                plot_seg_load(
                    sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], cs="global", pdf_file=pdf
                )
                # plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                # plot_overlap(
                #     sdata,
                #     KEYS[IMAGE_KEY],
                #     KEYS[SHAPES_KEY],
                #     KEYS[POINTS_KEY],
                #     KEYS[TABLE_KEY],
                #     image_scale_keys,
                #     cs="pixel",
                #     pdf_file=pdf,
                # )


def load_segmentation_all(
    config: SpidaConfig,
    exp_name: str,
    type: str = "vpt",
    prefix_name: str = "default",
    plot: bool = False,
    overwrite: bool = False,
    **load_kwargs,
):
    """
    Load one segmentation for every region of an experiment that has it under the segmentation store.
    Regions where it is already loaded are skipped unless ``overwrite``.

    Parameters:
    config (SpidaConfig): Store roots.
    exp_name (str): Name of the experiment.
    type (str): Type of the segmentation data to load (default is "vpt").
    prefix_name (str): Segmentation name (storage label + key prefix) (default is "default").
    overwrite (bool): Replace segmentations already loaded (default is False).
    """
    _METHOD = {"vpt": "cellpose", "cellpose": "cellpose", "mesmer": "mesmer", "proseg": "proseg"}
    label = prefix_name if prefix_name and prefix_name != "default" else _METHOD.get(type, type)
    for reg_name in TreeLayout(config).units(exp_name, under="seg", label=label):
        dataset = dataset_from(config, exp_name, reg_name)
        try:
            K = dataset.keys(prefix_name)
            is_writeable(dataset.zarr, [K[POINTS_KEY], K[SHAPES_KEY], K[TABLE_KEY]], overwrite=overwrite)
        except FileExistsError:
            logger.warning("Skipping %s: %s already loaded (pass --overwrite to replace)", reg_name, prefix_name)
            continue
        load_segmentation_region(
            dataset,
            type=type,
            prefix_name=prefix_name,
            plot=plot,
            overwrite=overwrite,
            **load_kwargs,
        )


def load_deconvolution_region(
    dataset: Dataset,
    image_dir: str | Path | None = None,
    image_name: str = "decon_image",
    suffix: str = ".decon.tif",
    z_layer: int | str = 3,
    plot: bool = False,
    overwrite: bool = False,
    **load_kwargs,
):
    """
    Load deconvolution images into the dataset's store.

    Parameters:
    dataset (Dataset): The dataset (store, raw images, figures dir).
    image_dir (str|Path|None): Directory holding the deconvolved images (default: the raw images dir).
    image_name (str): Name of the image key in the spatialdata object (default is "decon_image").
    suffix (str): Suffix of the deconvolution image files (default is ".decon.tif").
    plot (bool): Whether to generate plots (default is False).
    overwrite (bool): Replace ``image_name`` if it already exists (default is False).
    load_kwargs: Additional keyword arguments for the load_decon_images function.
    """
    image_dir = Path(image_dir) if image_dir is not None else dataset.raw_images()
    zarr_path = dataset.require("zarr")
    require(files=[image_dir])
    is_writeable(zarr_path, [image_name], overwrite=overwrite)
    logger.info("LOADING DECONVOLUTION IMAGES from %s into %s", image_dir, zarr_path)
    sdata = open_store(zarr_path)

    sdata = load_decon_images(
        sdata,
        image_dir=image_dir,
        image_name=image_name,
        suffix=suffix,
        z_layer=z_layer,
        overwrite=overwrite,
        **load_kwargs
    )

    if plot:
        image_channels = sd.models.get_channel_names(sdata[image_name])
        image_scale_keys = list(sdata[image_name].keys())

        image_path = dataset.figures_dir("default") / "pixi-decon.pdf"
        image_path.parent.mkdir(parents=True, exist_ok=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with PdfPages(image_path) as pdf:
                plot_images(
                    sdata,
                    image_name,
                    image_scale_keys,
                    image_channels,
                    cs="global",
                    pdf_file=pdf,
                )
