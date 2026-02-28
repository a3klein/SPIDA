import subprocess
import warnings
import logging
from pathlib import Path

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


# Loading VPT function
def load_vpt_segmentation(
    sdata: sd.SpatialData,
    exp_name: str,
    reg_name: str,
    vpt_path: str,
    prefix_name: str = "vpt",
    cell_metadata_fname: str = "cell_metadata.csv",
    cell_by_gene_fname: str = "cell_by_gene.csv",
    detected_transcripts_fname: str = "detected_transcripts.csv",
    cellpose_micron_space_fname: str = "cellpose_micron_space.parquet",
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
    Load the vpt segmentation into a spatialdata object.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    vpt_path (str): The path to the vpt segmentation output directory.
    cell_metadata_fname (str): The filename for the cell metadata (default is "cell_metadata.csv").
    cell_by_gene_fname (str): The filename for the cell by gene data (default is "cell_by_gene.csv").
    detected_transcripts_fname (str): The filename for the detected transcripts data (default is "detected_transcripts.csv").
    cellpose_micron_space_fname (str): The filename for the cellpose micron space data (default is "cellpose_micron_space.parquet").
    """

    logger.info(f"cell_metadata_fname={cell_metadata_fname}")
    logger.info(f"cell_by_gene_fname={cell_by_gene_fname}")
    logger.info(f"detected_transcripts_fname={detected_transcripts_fname}")
    logger.info(f"cellpose_micron_space_fname={cellpose_micron_space_fname}")
    for key, value in kwargs.items():
        logger.info(f"{key}={value}")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spatialdata_io.readers.merscope import _get_points, _get_table
    from .read_vpt import _get_polygons

    # KEYS
    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    identity = Identity()
    affine = get_transformation(
        sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global"
    )
    transformations = {"global": affine}

    # getting the shapes (but also the EntityID to regular ID mapping)
    # boundaries_path = f"{vpt_path}/{reg_name}/{cellpose_micron_space_fname}"
    # pols = _get_polygons(boundaries_path, transformations)
    # pols['ID'] = pols['ID'] + 1
    # entity_to_id_dict = pols['ID'].to_dict()
    # assert pols.shape[0] == len(entity_to_id_dict) == np.unique(pols['ID']).shape[0], \
    #     "The number of polygons does not match the number of unique IDs. Please check the data."

    # Getting the points
    points = {}
    transcripts_path = f"{vpt_path}/{reg_name}/{detected_transcripts_fname}"
    tz = _get_points(transcripts_path, transformations)
    # tz['cell_id'] = tz['cell_id'].astype(str).map(entity_to_id_dict).fillna(0).astype(int) # .compute()
    points[KEYS[POINTS_KEY]] = tz

    # Getting the shapes
    # shapes = {}
    # pols.set_index("ID", inplace=True)
    # shapes[KEYS[SHAPES_KEY]] = pols

    # # Getting the shapes
    shapes = {}
    boundaries_path = f"{vpt_path}/{reg_name}/{cellpose_micron_space_fname}"
    shapes[KEYS[SHAPES_KEY]] = _get_polygons(boundaries_path, transformations)

    # Getting the table
    tables = {}
    count_path = f"{vpt_path}/{reg_name}/{cell_by_gene_fname}"
    obs_path = f"{vpt_path}/{reg_name}/{cell_metadata_fname}"
    table = _get_table(
        count_path,
        obs_path,
        reg_name,
        exp_name,
        f"{exp_name}_{reg_name}",
        KEYS[SHAPES_KEY],
    )
    # table.obs.index = table.obs.index.map(entity_to_id_dict).astype(str)
    table.obs.index.name = "index"
    tables[KEYS[TABLE_KEY]] = table

    # adding the data to the spatialdata object
    sdata[KEYS[TABLE_KEY]] = tables[KEYS[TABLE_KEY]]
    sdata[KEYS[POINTS_KEY]] = points[KEYS[POINTS_KEY]]
    sdata[KEYS[SHAPES_KEY]] = shapes[KEYS[SHAPES_KEY]]
    sdata = _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=["EntityID"]
    )

    if (qc_regions is not None) or (qc_table_key is not None) or (qc_shapes_key is not None):
        logger.info("Applying transcript QC filter to loaded VPT segmentation elements.")
        sdata = apply_qc_filter_to_segmentation(
            sdata,
            table_key=KEYS[TABLE_KEY],
            shapes_key=KEYS[SHAPES_KEY],
            points_key=KEYS[POINTS_KEY],
            qc_regions=qc_regions,
            qc_table_key=qc_table_key,
            qc_shapes_key=qc_shapes_key,
            qc_filter_col=qc_filter_col,
            qc_pass_value=qc_pass_value,
            table_cell_id_col=table_cell_id_col,
            points_cell_id_col=points_cell_id_col,
        )

    # Doing the transformations:
    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    # saving data to disk
    # saving points
    if f"points/{KEYS[POINTS_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[POINTS_KEY])
    else:
        logger.warning(
            f"Points {KEYS[POINTS_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    # saving shapes
    if f"shapes/{KEYS[SHAPES_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[SHAPES_KEY])
    else:
        logger.warning(
            f"Shapes {KEYS[SHAPES_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    # saving table
    if f"tables/{KEYS[TABLE_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[TABLE_KEY])
    else:
        logger.warning(
            f"Table {KEYS[TABLE_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    if sd.__version__ <= "0.6.0": 
        sd.save_transformations(sdata)
    else:
        sdata.write_transformations()
    # sd.save_transformations(sdata)

    return sdata


@DeprecationWarning
def load_proseg_segmentation_v2(
    sdata: sd.SpatialData,
    exp_name: str,
    reg_name: str,
    proseg_path: str,
    prefix_name: str = "proseg",
    cell_metadata_fname: str = "cell-metadata.csv.gz",
    cell_by_gene_fname: str = "expected-counts.csv.gz",
    detected_transcripts_fname: str = "transcript-metadata.csv.gz",
    cell_polygons_fname: str = "cell-polygons.geojson.gz",
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
    Load the ProSeg segmentation data.
    This function is a placeholder for loading ProSeg segmentation data.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    proseg_path (str): The path to the ProSeg segmentation output directory.
    cell_metadata_fname (str): The filename for the cell metadata (default is "cell_metadata.csv.gz").
    cell_by_gene_fname (str): The filename for the cell by gene data (default is "expected-counts.csv.gz").
    detected_transcripts_fname (str): The filename for the detected transcripts
    cell_polygons_fname (str): The filename for the union polygons data (default is "cell-polygons.geojson.gz").
    """

    from spida.S.io.read_proseg import _get_polygons, _get_points, _get_table

    # KEYS
    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    identity = Identity()
    affine = get_transformation(
        sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global"
    )
    transformations = {"global": affine}

    # Getting the points
    points = {}
    transcripts_path = f"{proseg_path}/{reg_name}/{detected_transcripts_fname}"
    points[KEYS[POINTS_KEY]] = _get_points(transcripts_path, transformations)

    # Getting the shapes
    shapes = {}
    boundaries_path = f"{proseg_path}/{reg_name}/{cell_polygons_fname}"
    shapes[KEYS[SHAPES_KEY]] = _get_polygons(boundaries_path, transformations)

    # Getting the table
    tables = {}
    count_path = f"{proseg_path}/{reg_name}/{cell_by_gene_fname}"
    obs_path = f"{proseg_path}/{reg_name}/{cell_metadata_fname}"
    tables[KEYS[TABLE_KEY]] = _get_table(
        count_path,
        obs_path,
        reg_name,
        exp_name,
        f"{exp_name}_{reg_name}",
        KEYS[SHAPES_KEY],
    )
    tables[KEYS[TABLE_KEY]].obs.index.name = "index"

    # adding the data to the spatialdata object
    sdata[KEYS[TABLE_KEY]] = tables[KEYS[TABLE_KEY]]
    sdata[KEYS[POINTS_KEY]] = points[KEYS[POINTS_KEY]]
    sdata[KEYS[SHAPES_KEY]] = shapes[KEYS[SHAPES_KEY]]
    sdata = _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=["cell"]
    )

    if (qc_regions is not None) or (qc_table_key is not None) or (qc_shapes_key is not None):
        logger.info("Applying transcript QC filter to loaded ProSeg (v2) segmentation elements.")
        sdata = apply_qc_filter_to_segmentation(
            sdata,
            table_key=KEYS[TABLE_KEY],
            shapes_key=KEYS[SHAPES_KEY],
            points_key=KEYS[POINTS_KEY],
            qc_regions=qc_regions,
            qc_table_key=qc_table_key,
            qc_shapes_key=qc_shapes_key,
            qc_filter_col=qc_filter_col,
            qc_pass_value=qc_pass_value,
            table_cell_id_col=table_cell_id_col,
            points_cell_id_col=points_cell_id_col,
        )

    # Doing the transformations:
    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    # saving data to disk
    # saving data to disk
    # saving points
    if f"points/{KEYS[POINTS_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[POINTS_KEY])
    else:
        logger.warning(
            f"Points {KEYS[POINTS_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    # saving shapes
    if f"shapes/{KEYS[SHAPES_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[SHAPES_KEY])
    else:
        logger.warning(
            f"Shapes {KEYS[SHAPES_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    # saving table
    if f"tables/{KEYS[TABLE_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[TABLE_KEY])
    else:
        logger.warning(
            f"Table {KEYS[TABLE_KEY]} already exists in the spatialdata object. Not overwriting it."
        )

    if sd.__version__ <= "0.6.0":
        sd.save_transformations(sdata)
    else:
        sdata.write_transformations()

    return sdata

def load_proseg_segmentation_v3(
    sdata: sd.SpatialData,
    exp_name: str,
    reg_name: str,
    proseg_path: str,
    prefix_name: str = "proseg",
    zarr_name: str = "proseg_outputs.zarr",
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
    Load the ProSeg segmentation data.
    This function is a placeholder for loading ProSeg segmentation data.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    proseg_path (str): The path to the ProSeg segmentation output directory.
    zarr_name (str): The path to the zarr file to save the data.
    """
    from spida.S.io.read_proseg import _get_table_v3, parse_points
    logger.info("Loading ProSeg segmentation v3")

    # KEYS
    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    identity = Identity()
    affine = get_transformation(
        sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global"
    )
    transformations = {"global": affine}

    sdata_proseg = sd.read_zarr(f"{proseg_path}/{reg_name}/{zarr_name}")

    points = {}
    points[KEYS[POINTS_KEY]] = parse_points(sdata_proseg["transcripts"].copy(), None)
    sd.transformations.set_transformation(points[KEYS[POINTS_KEY]], transformation=affine, to_coordinate_system="global")

    shapes = {}
    shapes[KEYS[SHAPES_KEY]] = sdata_proseg["cell_boundaries"].set_crs(None, allow_override=True).copy() 
    shapes[KEYS[SHAPES_KEY]]['cell'] = shapes[KEYS[SHAPES_KEY]]['cell'].astype(str)
    sd.transformations.set_transformation(shapes[KEYS[SHAPES_KEY]], transformation=affine, to_coordinate_system="global")

    tables = {}
    tables[KEYS[TABLE_KEY]] = _get_table_v3(
        sdata_proseg["table"].copy(),
        reg_name,
        exp_name,
        f"{exp_name}_{reg_name}",
        KEYS[SHAPES_KEY],
    )
    tables[KEYS[TABLE_KEY]].obs.index.name = "index"

    sdata[KEYS[TABLE_KEY]] = tables[KEYS[TABLE_KEY]]
    sdata[KEYS[POINTS_KEY]] = points[KEYS[POINTS_KEY]]
    sdata[KEYS[SHAPES_KEY]] = shapes[KEYS[SHAPES_KEY]]
    sdata = _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=["cell"]
    )

    if (qc_regions is not None) or (qc_table_key is not None) or (qc_shapes_key is not None):
        logger.info("Applying transcript QC filter to loaded ProSeg (v3) segmentation elements.")
        sdata = apply_qc_filter_to_segmentation(
            sdata,
            table_key=KEYS[TABLE_KEY],
            shapes_key=KEYS[SHAPES_KEY],
            points_key=KEYS[POINTS_KEY],
            qc_regions=qc_regions,
            qc_table_key=qc_table_key,
            qc_shapes_key=qc_shapes_key,
            qc_filter_col=qc_filter_col,
            qc_pass_value=qc_pass_value,
            table_cell_id_col=table_cell_id_col,
            points_cell_id_col=points_cell_id_col,
        )

    # Doing the transformations:
    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")
    
    # saving points
    if f"points/{KEYS[POINTS_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[POINTS_KEY])
    else:
        logger.warning(
            f"Points {KEYS[POINTS_KEY]} already exists in the spatialdata object. Not overwriting it."
        )
    # saving shapes
    if f"shapes/{KEYS[SHAPES_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[SHAPES_KEY])
    else:
        logger.warning(
            f"Shapes {KEYS[SHAPES_KEY]} already exists in the spatialdata object. Not overwriting it."
        )
    # saving table
    if f"tables/{KEYS[TABLE_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[TABLE_KEY])
    else:
        logger.warning(
            f"Table {KEYS[TABLE_KEY]} already exists in the spatialdata object. Not overwriting it."
        )
    
    if sd.__version__ <= "0.6.0":   
        sd.save_transformations(sdata)
    else:
        sdata.write_transformations()

    return sdata


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
