import logging
from pathlib import Path
import warnings
import shutil
import pickle

import anndata as ad
import spatialdata as sd

from spida.utilities.sd_utils import _gen_keys
from spida._constants import (
    IMAGE_KEY,
    SHAPES_KEY,
    POINTS_KEY,
    TABLE_KEY,
    rename_exp_salk,
    rename_reg_salk,
    rename_exp_ucsd,
    rename_reg_ucsd,
    ren_to_exp_map,
)

logger = logging.getLogger(__package__)


def aggregate_experiments(
    zarr_path: Path,
    experiment_list: list[str],
    prefix_name: str,
    suffix: str = "_filt",
    lab_name: str = "salk",
    project_name: str = "BG_salk",
):
    """
    Aggregate spatial data from multiple experiments into a single SpatialData object.
    Args:
        zarr_path (Path): Path to the directory containing the zarr stores for each experiment
        experiments (list[str]): List of experiment names to aggregate
        prefix_ext (str): Prefix for the keys in the aggregated data
        suffix (str): Suffix to append to the table keys
        lab_name (str): Name of the lab, used for renaming experiments and regions
    Returns:
        SpatialData: Aggregated spatial data object containing images, shapes, points, and tables
    """

    if lab_name == "salk":
        # For renaming the experiments
        experiment_rename_func = rename_exp_salk
        region_rename_func = rename_reg_salk
    elif lab_name == "ucsd":
        # For renaming the experiments
        experiment_rename_func = rename_exp_ucsd
        region_rename_func = rename_reg_ucsd
    else:
        logger.warning(
            f"Lab name {lab_name} not recognized. Using default (salk) renaming functions."
        )
        experiment_rename_func = rename_exp_salk
        region_rename_func = rename_reg_salk

    images = {}
    shapes = {}
    points = {}
    tables = {}
    for e in experiment_list:
        experiment = Path(zarr_path) / e  # getting the full path to the experiment
        # for each experiment in the zarr store path
        exp_name = experiment.name  # getting the full experiment name
        ename = experiment_rename_func(
            exp_name
        )  # getting the shorthand experiment name
        if not experiment.is_dir():  # processing the experiment
            continue
        logger.info(f"Processing {exp_name}...")
        # for each region in the experiment
        for region in experiment.iterdir():
            region_name = region.name  # get the full region name
            rname = region_rename_func(region_name)  # get the shorthand region name
            if not region.is_dir():  # processing the region
                continue
            logger.info(f"  Processing {region_name}...")
            # Generate keys for the current region
            KEYS = _gen_keys(prefix_name, exp_name, region_name)
            # Load the spatial data for the current region
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # Read the spatial data
                sdata = sd.read_zarr(region)

            image_name = f"{ename}_{rname}_raw_images"
            decon_image_name = f"{ename}_{rname}_decon_images"
            shapes_name = f"{ename}_{rname}_shapes"
            table_name = f"{ename}_{rname}_tables"
            points_name = f"{ename}_{rname}_points"

            images[image_name] = sdata[KEYS[IMAGE_KEY]]
            images[decon_image_name] = sdata["decon_images"]
            shapes[shapes_name] = sdata[KEYS[SHAPES_KEY]]
            points[points_name] = sdata[KEYS[POINTS_KEY]]
            tables[table_name] = sdata[f"{KEYS[TABLE_KEY]}{suffix}"]

            tables[table_name].uns["spatialdata_attrs"]["region"] = shapes_name
            region_key = tables[table_name].uns["spatialdata_attrs"]["region_key"]
            tables[table_name].obs[region_key] = shapes_name
            tables[table_name].obs[region_key] = (
                tables[table_name].obs[region_key].astype("category")
            )
            
            if lab_name == "salk":
                tables[table_name].obs["brain_region"] = ename
            elif lab_name == "ucsd": # converting experiment number to region name
                tables[table_name].obs["ename"] = ename
                tables[table_name].obs["brain_region"] = ren_to_exp_map.get(ename, ename)
            tables[table_name].obs["donor"] = rname
            tables[table_name].obs["replicate"] = lab_name
            tables[table_name].obs['dataset_id'] = f"{ename}_{rname}_{lab_name}"
            # unique coordinate system for each experiment!
            for elem in [
                images[image_name],
                images[decon_image_name],
                shapes[shapes_name],
                points[points_name],
            ]:
                sd.transformations.set_transformation(
                    elem,
                    sd.transformations.get_transformation(
                        elem, to_coordinate_system="global"
                    ),
                    to_coordinate_system=f"{ename}_{rname}_global",
                )
                sd.transformations.remove_transformation(
                    elem, to_coordinate_system="global"
                )
                # If the transformation does not exist, we can ignore it
                try:
                    sd.transformations.remove_transformation(
                        elem, to_coordinate_system="pixel"
                    )
                except KeyError:
                    pass

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        sdata_comb = sd.SpatialData(
            images=images,
            shapes=shapes,
            points=points,
            tables=tables,
        )

        project_path = Path(zarr_path) / f"{project_name}_{lab_name}_{prefix_name}{suffix}"
        if project_path.exists():
            logger.warning(
                f"Project path {project_path} already exists. Overwriting it."
            )
            shutil.rmtree(project_path)
        project_path.parent.mkdir(parents=True, exist_ok=True)
        sdata_comb.write(project_path)
        sdata_comb = sd.read_zarr(project_path)

    return sdata_comb


def concatenate_tables(
    sdata: sd.SpatialData,
    table_keys: list[str] = None,
    output_path: Path = None,
):
    if table_keys is None:
        table_keys = list(sdata.tables.keys())

    uns_dict = {}
    concatenated_tables = []
    for key in table_keys:
        if key in sdata.tables:
            concatenated_tables.append(sdata.tables[key])
            uns_dict[key] = sdata.tables[key].uns.copy()
        else:
            logger.warning(
                f"Table key {key} not found in SpatialData object. Skipping."
            )
                        
    with open(output_path.parent / f"{output_path.name.split('.')[0]}_uns_dict.pkl", 'wb') as f:
        pickle.dump(uns_dict, f)
    
    # Need to double check
    adata = ad.concat(concatenated_tables)
    if output_path is not None:
        adata.write(output_path)
    else:
        logger.warning(
            "No output path provided. Returning concatenated tables without writing to disk."
        )

    return adata
