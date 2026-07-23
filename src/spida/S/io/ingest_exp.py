import subprocess
import warnings
import logging
from pathlib import Path

import pandas as pd
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd
    from spatialdata_io import merscope

    # Transformations
    from spatialdata.transformations import (
        Identity,
        get_transformation,
        set_transformation,
    )

from spida.utilities.sd_utils import _gen_keys
from spida._constants import SHAPES_KEY, POINTS_KEY, TABLE_KEY, IMAGE_KEY
from .filt_to_qc import apply_qc_filter_to_segmentation

logger = logging.getLogger(__package__)


### TODO:
#  Implement new logic when it comes to removing disconnected polygon geometries (OLD)
# Ideas:
# - Only keep the biggest polygon (now)
# - Only remove polygons with shapes smaller than a certain cutoff (artifacts) + create new cell artifacts from remaining disconnected geometries
# - Same as above but bridge together somehow between the two disconnected geometries
# - Same as above, but remove cells that have more than two disconnected geometries post filtering of small artifacts

def read_merscope(
    path,
    zarr_path,
    exp_name: str = None,
    reg_name: str = None,
    prefix_name: str = None,
    **kwargs,
):
    """
    Read Merscope data from a given path as a spatialdata object.
    Saves the data to a zarr file at the specified zarr_path.

    The function then preprocesses some of the spatialdata attributes to make it compatible with downstream analysis.

    Parameters:
    path (str): Path to the Merscope data directory.
    zarr_path (str): Path to save the zarr file.
    **kwargs: Additional keyword arguments to pass to the function.
    """
    # Generating spatialdata keys
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    # Reading in the experiment and then writing and reloading it as zarr format
    # surpressing the read_only warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata = merscope(path, slide_name=f"{prefix_name}_{exp_name}")
        if sd.__version__ <= "0.6.0": 
            sdata.write(zarr_path, overwrite=True)
            # Renaming the default table from "table" to "default_table"
            subprocess.run(
                [
                    "mv",
                    f"{zarr_path}/tables/table",
                    f"{zarr_path}/tables/{KEYS[TABLE_KEY]}",
                ],
                check=True,
            )
            sdata = sd.read_zarr(zarr_path, on_bad_files="warn")
        else: 
            sdata['default_table'] = sd.deepcopy(sdata['table'])
            del sdata['table']
            sdata.write(zarr_path, overwrite=True)
            sdata = sd.read_zarr(zarr_path, on_bad_files="warn")


    # All MultiPolygons to Polygons
    sdata = _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=["EntityID"]
    )

    # sdata[KEYS[TABLE_KEY]] = sd.deepcopy(sdata['table'])
    # sdata.delete_element_from_disk("table")
    # Making sure that there is no duplicates (columns / index) in the table
    adata = sdata[KEYS[TABLE_KEY]].copy()
    adata.obs.index.name = "index"
    sdata.delete_element_from_disk(KEYS[TABLE_KEY])  # Remove the old table from disk
    sdata[KEYS[TABLE_KEY]] = adata
    sdata.write_element(KEYS[TABLE_KEY])
    # sdata[KEYS[TABLE_KEY]].obs.index.name = "index"

    # Creating the transformations
    identity = Identity()
    affine = get_transformation(sdata[KEYS[POINTS_KEY]])
    affine_inv = affine.inverse()

    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[IMAGE_KEY]], affine_inv, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")
    if sd.__version__ <= "0.6.0": 
        sd.save_transformations(sdata)
    else: 
        sdata.write_transformations()

    return sdata


def _cast_multipolygons_to_polygons(
    sdata, shapes_key, subset_field: list = ["EntityID"]
):
    """
    Casts multipolygons to polygons in the spatialdata object.
    Remoces disconnected polygons by exploding multipolygons and dropping the smaller duplicates.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object containing the shapes.
    shapes_key (str): The key for the shapes in the spatialdata object.

    Returns:
    spatialdata.SpatialData: The modified spatialdata object with polygons.
    """
    gdf = sdata[shapes_key]  # getting the GeodataFrame from the spatialdata object
    if (
        "MultiPolygon" in gdf.geom_type.unique()
    ):  # if there are multipolygons in the geometry
        gdf = gdf.explode()
        gdf["area"] = gdf.area  # sort by area to remove smaller duplicates
        gdf.sort_values(by="area", ascending=False, inplace=True)
        index_dup = gdf.index[gdf.index.duplicated()]
        logger.info(
            f"{len(index_dup)} duplicated indices found in gdf"
        )  # print how many duplicates were found
        gdf = gdf.drop_duplicates(subset=subset_field)
        sdata[shapes_key] = gdf
    return sdata


def _resolve_file(
    seg_region: Path,
    canonical: str,
    legacy: str | None = None,
) -> Path:
    """
    Resolve a segmentation schema filename, preferring the canonical name and falling back to a
    method's legacy name (for pre-redesign data).

    Parameters:
    seg_region (Path): The region's segmentation-output directory.
    canonical (str): The segmentation schema filename to look for first.
    legacy (str | None): A legacy filename to fall back to if the canonical one is absent (default is None).

    Returns:
    Path: The resolved path (the canonical path if neither exists, so downstream errors point at the expected name).
    """
    p = seg_region / canonical
    if p.exists():
        return p
    if legacy is not None and (seg_region / legacy).exists():
        logger.info("using legacy file %r (canonical %r absent)", legacy, canonical)
        return seg_region / legacy
    return p


def _finalize_and_write(
    sdata: sd.SpatialData,
    KEYS: dict,
    identity,
    qc_kwargs: dict,
):
    """
    Shared loader tail: dedup multipolygons, optionally apply a transcript-QC
    filter, set the pixel coordinate system, then idempotently write the three
    elements and their transformations.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object being populated.
    KEYS (dict): The element-key mapping for this prefix/experiment/region.
    identity (Identity): The identity transformation used for the pixel coordinate system.
    qc_kwargs (dict): Transcript-QC filter arguments (applied only if any qc region/key is set).

    Returns:
    spatialdata.SpatialData: The finalized, on-disk-written spatialdata object.
    """
    from ._seg_readers import _CELL_KEY

    sdata = _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=[_CELL_KEY]
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

    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    for key, kind in (
        (KEYS[POINTS_KEY], "points"),
        (KEYS[SHAPES_KEY], "shapes"),
        (KEYS[TABLE_KEY], "tables"),
    ):
        if f"{kind}/{key}" not in sdata.elements_paths_on_disk():
            sdata.write_element(key)
        else:
            logger.warning("%s %r already exists on disk; not overwriting.", kind, key)

    if sd.__version__ <= "0.6.0":
        sd.save_transformations(sdata)
    else:
        sdata.write_transformations()
    return sdata


def _load_proseg_native(
    sdata: sd.SpatialData,
    spec,
    seg_region: Path,
    KEYS: dict,
    affine,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    cell_metadata_fname: str,
):
    """
    Load proseg output: counts + transcripts come from proseg's own output
    (zarr for v3, CSVs for v2), and shapes are the native 2D union polygons
    (kept for plotting). Everything is keyed on EntityID (== proseg ``cell``).

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    spec (SegmentationClass): The proseg method spec (selects v2 CSVs vs v3 zarr).
    seg_region (Path): The region's proseg-output directory.
    KEYS (dict): The element-key mapping for this prefix/experiment/region.
    affine: The global (micron) transformation taken from the default shapes element.
    reg_name (str): The name of the region.
    exp_name (str): The name of the experiment.
    dataset_id (str): The dataset identifier ("{exp}_{region}") stored on the table.
    cell_metadata_fname (str): The canonical obs (cell metadata) filename.

    Returns:
    spatialdata.SpatialData: The spatialdata object with proseg points/shapes/table added.
    """
    import dask.dataframe as dd
    from ._seg_readers import (
        get_table, get_table_from_adata, parse_points, _CELL_KEY,
    )
    from shapely.geometry import MultiPolygon

    # obs: segmentation schema cell_metadata.csv (EntityID) else legacy signals file
    obs_path = _resolve_file(seg_region, cell_metadata_fname, "cell_metadata_with_signals.csv")
    cell_meta = pd.read_csv(obs_path, index_col=0) if obs_path.exists() else None
    if cell_meta is not None:
        cell_meta.index = cell_meta.index.astype(int)
        cell_meta.index.name = _CELL_KEY

    zarr_path = seg_region / (spec.spatialdata_file or "__none__")
    if spec.spatialdata_file and zarr_path.exists():
        # ---- v3: zarr bundles table (counts) + transcripts (points) + union shapes ----
        proseg = sd.read_zarr(str(zarr_path))
        pts = parse_points(proseg["transcripts"].copy(), None)
        set_transformation(pts, affine, to_coordinate_system="global")
        sdata[KEYS[POINTS_KEY]] = pts

        shp = proseg["cell_boundaries"].set_crs(None, allow_override=True).copy()
        shp[_CELL_KEY] = shp["cell"].astype(int)
        shp.index = shp[_CELL_KEY].astype(str)
        set_transformation(shp, affine, to_coordinate_system="global")
        sdata[KEYS[SHAPES_KEY]] = shp

        ztab = proseg["table"].copy()
        if cell_meta is not None:                       # align obs to zarr's cell order
            cell_meta = cell_meta.reindex(ztab.obs["cell"].astype(int).values)
        table = get_table_from_adata(
            ztab, cell_meta, reg_name, exp_name, dataset_id, KEYS[SHAPES_KEY],
            cell_id_col=_CELL_KEY,
        )
    else:
        # ---- v2: native CSVs (expected-counts / transcript-metadata) + union geojson ----
        import gzip
        from spatialdata.models import ShapesModel
        transformations = {"global": affine}
        tx = dd.read_csv(str(seg_region / spec.transcripts_file), compression="gzip")
        sdata[KEYS[POINTS_KEY]] = parse_points(tx, transformations)

        union = seg_region / "cell-polygons.geojson.gz"    # 2D union (native proseg)
        with gzip.open(str(union), "rt") as f:
            shp = gpd.read_file(f)
        shp[_CELL_KEY] = shp["cell"].astype(int)
        shp.index = shp[_CELL_KEY].astype(str)
        shp.geometry = shp.geometry.map(
            lambda g: MultiPolygon(g.geoms) if g.geom_type == "MultiPolygon" else MultiPolygon([g])
        )
        sdata[KEYS[SHAPES_KEY]] = ShapesModel.parse(shp, transformations=transformations)

        table = get_table(
            seg_region / spec.counts_file, obs_path,
            reg_name, exp_name, dataset_id, KEYS[SHAPES_KEY], cell_id_col=_CELL_KEY,
        )
    table.obs.index.name = "index"
    sdata[KEYS[TABLE_KEY]] = table
    return sdata


def load_segmentation(
    sdata: sd.SpatialData,
    spec,
    exp_name: str,
    reg_name: str,
    seg_dir: str,
    prefix_name: str = "default",
    boundaries_fname: str | None = None,
    cell_by_gene_fname: str = "cell_by_gene.csv",
    detected_transcripts_fname: str = "detected_transcripts.csv",
    cell_metadata_fname: str = "cell_metadata.csv",
    qc_regions=None,
    qc_table_key: str | None = None,
    qc_shapes_key: str | None = None,
    qc_filter_col: str = "filtered",
    qc_pass_value: bool | int | str = False,
    table_cell_id_col: str | None = None,
    points_cell_id_col: str | None = None,
    **kwargs,
):
    """
    Load a segmenter's post-processed output into a spatialdata object (spec-driven).

    Replaces the old per-method ``load_vpt_segmentation`` /
    ``load_proseg_segmentation_v2`` / ``load_proseg_segmentation_v3``.
    Boundary-only methods (cellpose/mesmer) read the segmentation schema files;
    assignment-native methods (proseg) read counts/transcripts from their native
    output and use the 2D union polygons for shapes. All elements are keyed on
    EntityID; canonical filenames fall back to each method's legacy name.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    spec (SegmentationClass): The method spec (from ``get_spec``) driving file names and the read path.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    seg_dir (str): The path to the segmentation-output directory (contains ``{reg_name}/``).
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    boundaries_fname (str | None): Override for the boundaries filename; defaults to the segmentation schema name (default is None).
    cell_by_gene_fname (str): The cell-by-gene filename for boundary-only methods (default is "cell_by_gene.csv").
    detected_transcripts_fname (str): The detected-transcripts filename for boundary-only methods (default is "detected_transcripts.csv").
    cell_metadata_fname (str): The obs (cell metadata) filename (default is "cell_metadata.csv").
    qc_regions: Transcript-QC regions to filter by (default is None).
    qc_table_key (str | None): Table key holding transcript-QC pass/fail info (default is None).
    qc_shapes_key (str | None): Shapes key holding transcript-QC regions (default is None).
    qc_filter_col (str): The QC column used to filter (default is "filtered").
    qc_pass_value (bool | int | str): The value in ``qc_filter_col`` marking a passing cell (default is False).
    table_cell_id_col (str | None): Override for the table's cell-id column during QC (default is None).
    points_cell_id_col (str | None): Override for the points' cell-id column during QC (default is None).

    Returns:
    spatialdata.SpatialData: The spatialdata object with the loaded segmentation elements written to disk.
    """
    from spida.S.segmentation.backends.base import SCHEMA_BOUNDARIES
    from ._seg_readers import get_polygons, get_table, _CELL_KEY

    logger.info("load_segmentation: %s v%s -> region %s", spec.name, spec.version, reg_name)
    for k, v in kwargs.items():
        logger.info("  extra kwarg %s=%s", k, v)

    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)
    identity = Identity()
    affine = get_transformation(sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global")
    transformations = {"global": affine}
    seg_region = Path(seg_dir) / reg_name
    dataset_id = f"{exp_name}_{reg_name}"

    qc_kwargs = dict(
        qc_regions=qc_regions, qc_table_key=qc_table_key, qc_shapes_key=qc_shapes_key,
        qc_filter_col=qc_filter_col, qc_pass_value=qc_pass_value,
        table_cell_id_col=table_cell_id_col, points_cell_id_col=points_cell_id_col,
    )

    if spec.provides_assignment:
        sdata = _load_proseg_native(
            sdata, spec, seg_region, KEYS, affine, reg_name, exp_name,
            dataset_id, cell_metadata_fname,
        )
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spatialdata_io.readers.merscope import _get_points

        boundaries = _resolve_file(
            seg_region, boundaries_fname or SCHEMA_BOUNDARIES, spec.legacy_boundaries_file
        )
        sdata[KEYS[POINTS_KEY]] = _get_points(
            str(seg_region / detected_transcripts_fname), transformations
        )
        sdata[KEYS[SHAPES_KEY]] = get_polygons(
            boundaries, transformations, cell_id_col=_CELL_KEY
        )
        table = get_table(
            seg_region / cell_by_gene_fname, seg_region / cell_metadata_fname,
            reg_name, exp_name, dataset_id, KEYS[SHAPES_KEY], cell_id_col=_CELL_KEY,
        )
        table.obs.index.name = "index"
        sdata[KEYS[TABLE_KEY]] = table

    return _finalize_and_write(sdata, KEYS, identity, qc_kwargs)


def load_decon_images(
    sdata: sd.SpatialData,
    image_dir: str | Path,
    z_layer: int | str = 3,
    image_name: str = "decon_images",
    suffix: str = ".decon.tif",
    **image_models_kwargs,
):
    """
    Load deconvolution images into the spatialdata object.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the images into.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    image_path (str|Path): The path to the deconvolution images.
    """

    image_models_kwargs = {}
    image_models_kwargs["chunks"] = (1, 4906, 4906)
    image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    sdata[image_name] = _read_decon_images(
        image_dir, suffix=suffix,  z_layer=z_layer, **image_models_kwargs
    )
    if f"images/{image_name}" not in sdata.elements_paths_on_disk():
        sdata.write_element(image_name)
        sdata.write_metadata(image_name)
    else:
        logger.warning(
            f"Images {image_name} already exists in the spatialdata object. Not overwriting it."
        )

    return sdata


def _read_decon_images(
    image_dir,
    suffix: str = ".decon.tif",
    z_layer: int | str = 3,
    **image_models_kwargs):
    """
    Read deconvolution images from the specified directory.

    Parameters:
    suffix (str): The suffix of the image files (default is ".decon.tif").
    z_layer (int): The z-layer to read (default is 3).
    image_dir (str|Path): The directory containing the deconvolution images.
    **image_models_kwargs: Additional keyword arguments for the image model.

    Returns:
    Image2DModel: The parsed image model.
    """
    import os
    import numpy as np
    from spatialdata.models import Image2DModel
    from dask import array as da
    from dask_image.imread import imread

    if isinstance(image_dir, str):
        image_dir = Path(image_dir)
    stainings = ["DAPI", "PolyT"]
    z_stack_format = np.any([True for _fp in os.listdir(image_dir) if "z_stack" in _fp])
    if z_stack_format:
        im = da.stack(
            [
                imread(image_dir / f"mosaic_{stain}_z_stack_z{z_layer}{suffix}").squeeze()
                for stain in stainings
            ],
            axis=0,
        )
    else:
        im = da.stack(
            [
                imread(image_dir / f"mosaic_{stain}_z{z_layer}{suffix}").squeeze()
                for stain in stainings
            ],
            axis=0,
        )
    

    return Image2DModel.parse(
        im, dims=("c", "y", "x"), c_coords=stainings, rgb=None, **image_models_kwargs
    )
