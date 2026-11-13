"""
Build spatialdata elements from files.
Ingest an experiment, load in segmentations, deconvolved images, etc...
"""

from __future__ import annotations

import warnings
import logging
from pathlib import Path

import pandas as pd
import anndata as ad
import geopandas as gpd
from shapely.geometry import MultiPolygon
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd
    from spatialdata_io import merscope
    from spatialdata.models import PointsModel, ShapesModel, TableModel

    # Transformations
    from spatialdata.transformations import (
        BaseTransformation,
        Identity,
        get_transformation,
        set_transformation,
    )

from spida._constants import (
    SHAPES_KEY, POINTS_KEY, TABLE_KEY, IMAGE_KEY,
    DEFAULT_PRESET, PROSEG_PRESET, REGION_KEY, GEOMETRY_COL,
)
from spida.io.schema.segmentation import read_boundaries, resolve
from spida.io.store import write_element

logger = logging.getLogger(__package__)

# cell-id key across all segmentation schema files (proseg cell id == EntityID)
_CELL_KEY = DEFAULT_PRESET["METADATA_CELL_KEY"]  # "EntityID"


def _cast_to_mp(x):
    try:
        return MultiPolygon(x.geoms)
    except Exception:
        return MultiPolygon([x.convex_hull])


def get_polygons(
    boundaries_path: str | Path,
    transformations: dict[str, BaseTransformation],
    *,
    cell_id_col: str = _CELL_KEY,
) -> gpd.GeoDataFrame:
    """Parse a boundaries file into a ShapesModel, indexed by ``cell_id_col`` (str).

    Multi-z segmentation schema boundaries are collapsed downstream by
    ``_cast_multipolygons_to_polygons`` (dedup-to-largest per cell id).
    """
    gdf = read_boundaries(boundaries_path)
    if gdf.geometry.name != GEOMETRY_COL:
        gdf = gdf.rename_geometry(GEOMETRY_COL)
    gdf.geometry = gdf.geometry.make_valid()
    gdf.geometry = gdf.geometry.map(_cast_to_mp)
    gdf.index = gdf[cell_id_col].astype(str)
    # spatialdata promotes 2D polygons to 3D if any column is literally named "z"
    if "z" in gdf.columns:
        gdf = gdf.rename(columns={"z": "z_layer"})
    return ShapesModel.parse(gdf, transformations=transformations)


def get_table(
    count_path: str | Path,
    obs_path: str | Path,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    region: str,
    *,
    cell_id_col: str = _CELL_KEY,
) -> ad.AnnData:
    """Build the annotated table from a cell-by-gene CSV + an obs (metadata) CSV,
    keyed on ``cell_id_col``. Blank barcodes are moved to ``obsm['blank']``.
    (pandas infers ``.gz`` compression from the filename.)"""
    data = pd.read_csv(count_path, index_col=0)
    obs = pd.read_csv(obs_path, index_col=0)
    obs.index = obs.index.astype(str)
    data.index = data.index.astype(str)

    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = ad.AnnData(data.loc[:, is_gene], dtype=data.values.dtype, obs=obs)
    adata.obsm["blank"] = data.loc[:, ~is_gene]
    if {"center_x", "center_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
    elif {"centroid_x", "centroid_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["centroid_x", "centroid_y"]].values

    return _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col)


def get_table_from_adata(
    adata: ad.AnnData,
    obs_df: pd.DataFrame | None,
    reg_name: str,
    exp_name: str,
    dataset_id: str,
    region: str,
    *,
    cell_id_col: str = _CELL_KEY,
) -> ad.AnnData:
    """Build the table from a counts AnnData (e.g. proseg's zarr ``table``),
    optionally overriding obs with ``obs_df`` (the segmentation schema ``cell_metadata.csv``).
    Blank genes go to ``obsm['blank']``; keyed on ``cell_id_col``."""
    adata = adata.copy()
    if obs_df is not None:
        obs_df = obs_df.copy()
        obs_df.index = obs_df.index.astype(str)
        adata.obs = obs_df

    gene = adata.var["gene"] if "gene" in adata.var.columns else pd.Series(adata.var_names, index=adata.var_names)
    is_blank = gene.str.startswith("Blank-").to_numpy()
    blank = adata[:, is_blank].copy()
    adata = adata[:, ~is_blank].copy()
    adata.obs["transcript_count"] = adata.X.toarray().sum(axis=1).flatten()
    adata.obsm["blank"] = blank.X.toarray()
    if {"center_x", "center_y"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
    elif {"centroid_x", "centroid_y"}.issubset(adata.obs.columns):  # proseg-native obs (no schema metadata)
        adata.obsm["spatial"] = adata.obs[["centroid_x", "centroid_y"]].values
    adata.uns.pop("spatialdata_attrs", None)  # re-parse from scratch

    return _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col)


def _tag_and_parse(adata, reg_name, exp_name, dataset_id, region, cell_id_col):
    """Attach the standard obs bookkeeping columns and parse into a TableModel."""
    adata.obs["region"] = pd.Series(reg_name, index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(exp_name, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(dataset_id, index=adata.obs_names, dtype="category")
    adata.obs[REGION_KEY] = pd.Series(region, index=adata.obs_names, dtype="category")
    adata.obs[cell_id_col] = adata.obs.index
    return TableModel.parse(
        adata, region_key=REGION_KEY, region=region, instance_key=cell_id_col
    )


def parse_points(transcripts_df, transformations):
    """Parse a proseg transcript table (zarr ``transcripts`` or CSV) into a
    PointsModel using the proseg column preset."""
    gene_col = PROSEG_PRESET["tz_col"]["gene_col"]
    transcripts_df[gene_col] = transcripts_df[gene_col].astype("string")
    return PointsModel.parse(
        transcripts_df,
        coordinates={
            "x": PROSEG_PRESET["tz_col"]["x_col"],
            "y": PROSEG_PRESET["tz_col"]["y_col"],
        },
        transformations=transformations,
        feature_key=gene_col,
    )


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
    keys: dict,
    namespace: str,
    exp_name: str | None = None,
    reg_name: str | None = None,
    prefix_name: str = "default",
    **kwargs,
):
    """
    Read a MERSCOPE region into a new zarr store at ``zarr_path``.

    Element names follow ``keys``. Without --exp/--reg (empty ``namespace``) the ``<exp>_<reg>`` part
    is dropped from the spatialdata-io element names. The names used are recorded in ``sdata.attrs["spida"]``.

    Parameters:
    path (str): MERSCOPE region directory.
    zarr_path (str): Zarr store to write.
    keys (dict): Element keys (``Dataset.keys(prefix_name)``).
    namespace (str): The ``<exp>_<reg>`` part of element names ("" when ingesting without names).
    """
    KEYS = keys
    slide = f"{prefix_name}_{exp_name}" if exp_name else prefix_name
    region = reg_name or "region"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata = merscope(path, slide_name=slide, region_name=region)
        if not namespace:
            _drop_namespace(sdata, dataset_id=f"{slide}_{region}", prefix_name=prefix_name)
        sdata[KEYS[TABLE_KEY]] = sd.deepcopy(sdata["table"])
        del sdata["table"]
        sdata.attrs["spida"] = {"namespace": namespace, "exp": exp_name, "reg": reg_name}
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
    sdata.write_transformations()

    return sdata


def _drop_namespace(sdata, dataset_id: str, prefix_name: str) -> None:
    """Rename spatialdata-io's ``{dataset_id}_*`` elements to ``{prefix_name}_*`` (no ``<exp>_<reg>``), in memory."""
    for kind in ("images", "points", "shapes"):
        for name in list(getattr(sdata, kind)):
            if name.startswith(f"{dataset_id}_"):
                sdata[f"{prefix_name}{name[len(dataset_id):]}"] = sdata[name]
                del sdata[name]
    if "table" in sdata.tables:
        table, region = sdata["table"], f"{prefix_name}_polygons"
        attrs = table.uns["spatialdata_attrs"]
        table.obs[attrs["region_key"]] = pd.Categorical([region] * table.n_obs)
        sdata.set_table_annotates_spatialelement(
            "table", region=region,
            region_key=attrs["region_key"], instance_key=attrs["instance_key"],
        )


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


def write_segmentation(
    sdata: sd.SpatialData,
    KEYS: dict,
    overwrite: bool = False,
):
    """
    Set the pixel coordinate system on a read segmentation's points and shapes, then write
    the three elements (replaced with ``overwrite``).

    Parameters:
    sdata (spatialdata.SpatialData): The store's spatialdata object holding the segmentation elements.
    KEYS (dict): The element-key mapping for this prefix/experiment/region.
    overwrite (bool): Replace elements already on disk (default is False).

    Returns:
    spatialdata.SpatialData: The store re-read after the write.
    """
    identity = Identity()
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    for key in (KEYS[POINTS_KEY], KEYS[SHAPES_KEY], KEYS[TABLE_KEY]):
        write_element(sdata.path, key, sdata[key], overwrite=overwrite)
    return sd.read_zarr(sdata.path)


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

    # obs: segmentation schema cell_metadata.csv (EntityID) else legacy signals file
    obs_path = resolve(seg_region, cell_metadata_fname, "cell_metadata_with_signals.csv")
    cell_meta = pd.read_csv(obs_path, index_col=0) if obs_path.exists() else None
    if cell_meta is not None:
        cell_meta.index = cell_meta.index.astype(int)
        cell_meta.index.name = _CELL_KEY

    native = spec.native_reader("assigned")(seg_region, spec)
    if native["version"] == 3:
        # ---- v3: zarr bundles table (counts) + transcripts (points) + union shapes ----
        pts = parse_points(native["transcripts"].copy(), None)
        set_transformation(pts, affine, to_coordinate_system="global")
        sdata[KEYS[POINTS_KEY]] = pts

        shp = native["cells"].set_crs(None, allow_override=True).copy()
        shp[_CELL_KEY] = shp["cell"].astype(int)
        shp.index = shp[_CELL_KEY].astype(str)
        set_transformation(shp, affine, to_coordinate_system="global")
        sdata[KEYS[SHAPES_KEY]] = shp

        ztab = native["counts"].copy()
        if cell_meta is not None:                       # align obs to zarr's cell order
            cell_meta = cell_meta.reindex(ztab.obs["cell"].astype(int).values)
        table = get_table_from_adata(
            ztab, cell_meta, reg_name, exp_name, dataset_id, KEYS[SHAPES_KEY],
            cell_id_col=_CELL_KEY,
        )
    else:
        # ---- v2: native CSVs (expected-counts / transcript-metadata) + union geojson ----
        from spatialdata.models import ShapesModel
        transformations = {"global": affine}
        sdata[KEYS[POINTS_KEY]] = parse_points(native["transcripts"], transformations)

        shp = native["cells"]
        shp[_CELL_KEY] = shp["cell"].astype(int)
        shp.index = shp[_CELL_KEY].astype(str)
        shp.geometry = shp.geometry.map(
            lambda g: MultiPolygon(g.geoms) if g.geom_type == "MultiPolygon" else MultiPolygon([g])
        )
        sdata[KEYS[SHAPES_KEY]] = ShapesModel.parse(shp, transformations=transformations)

        table = get_table(
            native["counts_file"], obs_path,
            reg_name, exp_name, dataset_id, KEYS[SHAPES_KEY], cell_id_col=_CELL_KEY,
        )
    table.obs.index.name = "index"
    sdata[KEYS[TABLE_KEY]] = table
    return sdata


def read_segmentation(
    sdata: sd.SpatialData,
    spec,
    keys: dict,
    default_keys: dict,
    seg_dir: str,
    exp_name: str | None = None,
    reg_name: str | None = None,
    dataset_id: str = "",
    prefix_name: str = "default",
    boundaries_fname: str | None = None,
    cell_by_gene_fname: str = "cell_by_gene.csv",
    detected_transcripts_fname: str = "detected_transcripts.csv",
    cell_metadata_fname: str = "cell_metadata.csv",
    **kwargs,
):
    """
    Read a segmenter's post-processed output into a spatialdata object, in memory (spec-driven).
    Multipolygons are reduced to their largest polygon. Nothing is written; see ``write_segmentation``.

    Replaces the old per-method ``load_vpt_segmentation`` /
    ``load_proseg_segmentation_v2`` / ``load_proseg_segmentation_v3``.
    Boundary-only methods (cellpose/mesmer) read the segmentation schema files;
    assignment-native methods (proseg) read counts/transcripts from their native
    output and use the 2D union polygons for shapes. All elements are keyed on
    EntityID; canonical filenames fall back to each method's legacy name.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the segmentation into.
    spec (SegmentationClass): The method spec (from ``get_spec``) driving file names and the read path.
    keys (dict): Element keys for this segmentation (``Dataset.keys(prefix_name)``).
    default_keys (dict): Element keys of the ingested ``default`` elements (source of the transform).
    seg_dir (str): This segmentation's output directory.
    exp_name, reg_name (str | None): Written to the table's ``slide`` / ``region`` columns ("" when unnamed).
    dataset_id (str): Written to the table's ``dataset_id`` column.
    prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
    boundaries_fname (str | None): Override for the boundaries filename; defaults to the segmentation schema name (default is None).
    cell_by_gene_fname (str): The cell-by-gene filename for boundary-only methods (default is "cell_by_gene.csv").
    detected_transcripts_fname (str): The detected-transcripts filename for boundary-only methods (default is "detected_transcripts.csv").
    cell_metadata_fname (str): The obs (cell metadata) filename (default is "cell_metadata.csv").

    Returns:
    spatialdata.SpatialData: The spatialdata object with the segmentation's points, shapes and table (not written).
    """
    from spida.io.schema.segmentation import SCHEMA_BOUNDARIES

    logger.info("read_segmentation: %s v%s -> region %s", spec.name, spec.version, reg_name)
    for k, v in kwargs.items():
        logger.info("  extra kwarg %s=%s", k, v)

    KEYS = keys
    exp_name, reg_name = exp_name or "", reg_name or ""
    affine = get_transformation(sdata[default_keys[SHAPES_KEY]], to_coordinate_system="global")
    transformations = {"global": affine}

    seg_region = Path(seg_dir)

    if spec.provides_assignment:
        sdata = _load_proseg_native(
            sdata, spec, seg_region, KEYS, affine, reg_name, exp_name,
            dataset_id, cell_metadata_fname,
        )
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spatialdata_io.readers.merscope import _get_points

        boundaries = resolve(
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

    return _cast_multipolygons_to_polygons(
        sdata, KEYS[SHAPES_KEY], subset_field=[_CELL_KEY]
    )


def load_decon_images(
    sdata: sd.SpatialData,
    image_dir: str | Path,
    z_layer: int | str = 3,
    image_name: str = "decon_images",
    suffix: str = ".decon.tif",
    overwrite: bool = False,
    **image_models_kwargs,
):
    """
    Load deconvolution images into the spatialdata object.

    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object to load the images into.
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    image_path (str|Path): The path to the deconvolution images.
    overwrite (bool): Replace ``image_name`` if it already exists (default is False).
    """

    image_models_kwargs.setdefault("chunks", (1, 4906, 4906))
    image_models_kwargs.setdefault("scale_factors", [2, 2, 2, 2])

    image = _read_decon_images(
        image_dir, suffix=suffix,  z_layer=z_layer, **image_models_kwargs
    )
    write_element(sdata.path, image_name, image, overwrite=overwrite)
    return sd.read_zarr(sdata.path)


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
    from spida.io.schema.decon import Z_STACK_TAG, plane_name

    if isinstance(image_dir, str):
        image_dir = Path(image_dir)
    stainings = ["DAPI", "PolyT"]
    z_stack_format = np.any([True for _fp in os.listdir(image_dir) if Z_STACK_TAG in _fp])
    if z_stack_format:
        im = da.stack(
            [
                imread(image_dir / plane_name(stain, z_layer, suffix, z_stack=True)).squeeze()
                for stain in stainings
            ],
            axis=0,
        )
    else:
        im = da.stack(
            [
                imread(image_dir / plane_name(stain, z_layer, suffix)).squeeze()
                for stain in stainings
            ],
            axis=0,
        )
    

    return Image2DModel.parse(
        im, dims=("c", "y", "x"), c_coords=stainings, rgb=None, **image_models_kwargs
    )
