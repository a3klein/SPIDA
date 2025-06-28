import sys
import subprocess
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd
    from spatialdata_io import merscope

    # Transformations 
    from spatialdata.transformations import (
        Identity,
        get_transformation, 
        set_transformation
    )

spida_path = "/ceph/cephatlas/aklein/spida/src"
sys.path.append(spida_path)
from spida._utilities import _gen_keys
from spida._constants import SHAPES_KEY, POINTS_KEY, TABLE_KEY, IMAGE_KEY


### TODO: 
#  Implement new logic when it comes to removing disconnected polygon geometries
# Ideas: 
# - Only keep the biggest polygon (now) 
# - Only remove polygons with shapes smaller than a certain cutoff (artifacts) + create new cell artifacts from remaining disconnected geometries
# - Same as above but bridge together somehow between the two disconnected geometries 
# - Same as above, but remove cells that have more than two disconnected geometries post filtering of small artifacts


def read_merscope(path, zarr_path, exp_name:str=None, reg_name:str=None, prefix_name:str=None, **kwargs): 
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
        sdata.write(zarr_path, overwrite=True)
        # Renaming the default table from "table" to "default_table"
        subprocess.run(["mv", f"{zarr_path}/tables/table", f"{zarr_path}/tables/{KEYS[TABLE_KEY]}"], check=True)
        sdata = sd.read_zarr(zarr_path)

    # All MultiPolygons to Polygons
    sdata = _cast_multipolygons_to_polygons(sdata, KEYS[SHAPES_KEY], subset_field=['EntityID'])

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
    sd.save_transformations(sdata)

    return sdata

    
def _cast_multipolygons_to_polygons(sdata, shapes_key, subset_field:list=['EntityID']):
    """
    Casts multipolygons to polygons in the spatialdata object. 
    Remoces disconnected polygons by exploding multipolygons and dropping the smaller duplicates.
    
    Parameters:
    sdata (spatialdata.SpatialData): The spatialdata object containing the shapes.
    shapes_key (str): The key for the shapes in the spatialdata object.
    
    Returns:
    spatialdata.SpatialData: The modified spatialdata object with polygons.
    """
    gdf = sdata[shapes_key] # getting the GeodataFrame from the spatialdata object
    if "MultiPolygon" in gdf.geom_type.unique(): # if there are multipolygons in the geometry
        gdf = gdf.explode() 
        gdf['area'] = gdf.area # sort by area to remove smaller duplicates
        gdf.sort_values(by="area", ascending=False, inplace=True) 
        index_dup = gdf.index[gdf.index.duplicated()]
        print(len(index_dup), "duplicated indices found in gdf") # print how many duplicates were found 
        ### TODO: Depending on the preset might need to change the subset field to a parameter
        gdf = gdf.drop_duplicates(subset=subset_field)
        sdata[shapes_key] = gdf
    return sdata
    

# Loading VPT function
def load_vpt_segmentation(sdata:sd.SpatialData, 
                          exp_name:str, 
                          reg_name:str,
                          vpt_path:str,
                          prefix_name:str="vpt",
                          cell_metadata_fname:str="cell_metadata.csv",
                          cell_by_gene_fname:str="cell_by_gene.csv",
                          detected_transcripts_fname:str="detected_transcripts.csv",
                          cellpose_micron_space_fname:str="cellpose_micron_space.parquet",
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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        from spatialdata_io.readers.merscope import _get_polygons, _get_points, _get_table

    # KEYS
    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    identity = Identity()
    affine = get_transformation(sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global")
    transformations = {'global' : affine}

    # Getting the points 
    points = {}
    transcripts_path = f"{vpt_path}/{detected_transcripts_fname}"
    points[KEYS[POINTS_KEY]] = _get_points(transcripts_path, transformations)
    
    # Getting the shapes
    shapes = {}
    boundaries_path = f"{vpt_path}/{cellpose_micron_space_fname}"
    shapes[KEYS[SHAPES_KEY]] = _get_polygons(boundaries_path, transformations)
    
    # Getting the table
    tables = {}
    count_path = f"{vpt_path}/{cell_by_gene_fname}"
    obs_path = f"{vpt_path}/{cell_metadata_fname}"
    tables[KEYS[TABLE_KEY]] = _get_table(count_path, obs_path, reg_name, exp_name, f"{exp_name}_{reg_name}", KEYS[SHAPES_KEY])
    tables[KEYS[TABLE_KEY]].obs.index.name = "index"

    # adding the data to the spatialdata object
    sdata[KEYS[TABLE_KEY]] = tables[KEYS[TABLE_KEY]]
    sdata[KEYS[POINTS_KEY]] = points[KEYS[POINTS_KEY]]
    sdata[KEYS[SHAPES_KEY]] = shapes[KEYS[SHAPES_KEY]]
    sdata = _cast_multipolygons_to_polygons(sdata, KEYS[SHAPES_KEY], subset_field=['EntityID'])

    # Doing the transformations: 
    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    # saving data to disk
    if f"points/{KEYS[POINTS_KEY]}" not in sdata.elements_paths_on_disk(): 
        sdata.write_element(KEYS[POINTS_KEY])
    if f"shapes/{KEYS[SHAPES_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[SHAPES_KEY])
    if f"tables/{KEYS[TABLE_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[TABLE_KEY])
    
    sd.save_transformations(sdata)
    
    return sdata
    

def load_proseg_segmentation(sdata:sd.SpatialData,
                             exp_name:str, 
                             reg_name:str,
                             proseg_path:str,
                             prefix_name:str="proseg",
                             cell_metadata_fname:str="cell-metadata.csv.gz",
                             cell_by_gene_fname:str="expected-counts.csv.gz",
                             detected_transcripts_fname:str="transcript-metadata.csv.gz",
                             cell_polygons_fname:str="cell-polygons.geojson.gz",
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

    from spida.io.read_proseg import _get_polygons, _get_points, _get_table

    # KEYS
    DEF_KEYS = _gen_keys("default", exp_name, reg_name)
    KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    identity = Identity()
    affine = get_transformation(sdata[DEF_KEYS[SHAPES_KEY]], to_coordinate_system="global")
    transformations = {'global' : affine}

    # Getting the points 
    points = {}
    transcripts_path = f"{proseg_path}/{detected_transcripts_fname}"
    points[KEYS[POINTS_KEY]] = _get_points(transcripts_path, transformations)
    
    # Getting the shapes
    shapes = {}
    boundaries_path = f"{proseg_path}/{cell_polygons_fname}"
    shapes[KEYS[SHAPES_KEY]] = _get_polygons(boundaries_path, transformations)
    
    # Getting the table
    tables = {}
    count_path = f"{proseg_path}/{cell_by_gene_fname}"
    obs_path = f"{proseg_path}/{cell_metadata_fname}"
    tables[KEYS[TABLE_KEY]] = _get_table(count_path, obs_path, reg_name, exp_name, f"{exp_name}_{reg_name}", KEYS[SHAPES_KEY])
    tables[KEYS[TABLE_KEY]].obs.index.name = "index"

    # adding the data to the spatialdata object
    sdata[KEYS[TABLE_KEY]] = tables[KEYS[TABLE_KEY]]
    sdata[KEYS[POINTS_KEY]] = points[KEYS[POINTS_KEY]]
    sdata[KEYS[SHAPES_KEY]] = shapes[KEYS[SHAPES_KEY]]
    sdata = _cast_multipolygons_to_polygons(sdata, KEYS[SHAPES_KEY], subset_field=['cell'])

    # Doing the transformations: 
    # Setting the pixel space coordinate system
    set_transformation(sdata[KEYS[POINTS_KEY]], identity, to_coordinate_system="pixel")
    set_transformation(sdata[KEYS[SHAPES_KEY]], identity, to_coordinate_system="pixel")

    # saving data to disk
    if f"points/{KEYS[POINTS_KEY]}" not in sdata.elements_paths_on_disk(): 
        sdata.write_element(KEYS[POINTS_KEY])
    if f"shapes/{KEYS[SHAPES_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[SHAPES_KEY])
    if f"tables/{KEYS[TABLE_KEY]}" not in sdata.elements_paths_on_disk():
        sdata.write_element(KEYS[TABLE_KEY])
    
    sd.save_transformations(sdata)
    
    return sdata