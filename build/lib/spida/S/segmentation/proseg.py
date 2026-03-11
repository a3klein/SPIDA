import os
from dotenv import load_dotenv  # type: ignore

import logging
import warnings
import pathlib
import subprocess
import re
import gzip

import pandas as pd
import geopandas as gpd

load_dotenv()
logger = logging.getLogger(__package__)


def _determine_proseg_version():
    """
    Determine the version of the proseg binary installed.
    Returns:
        str: The version of the proseg binary.
    """
    try:
        ret = subprocess.run(["proseg", "-V"], capture_output=True, check=True)
        if ret.returncode != 0:
            raise ValueError("Proseg failed to run. Check the installation and the PATH variable.")
        proseg_version = ret.stdout.decode('utf-8').split(" ")[-1].strip()
        return proseg_version
    except FileNotFoundError:
        raise FileNotFoundError("Proseg binary not found. Please ensure it is installed and in your PATH.")

# click? 
def get_proseg_version():
    """
    Get the version of the proseg binary installed.
    Returns:
        str: The version of the proseg binary.
    """
    try:
        proseg_version = _determine_proseg_version()
        logger.info(f"Proseg version: {proseg_version}")
        return proseg_version
    except Exception as e:
        logger.error(f"Error determining Proseg version: {e}")
        raise

def _add_proseg_binary(rust_bin_path : str | None = None):
    """
    Add the path to the proseg binary to the PATH environment variable.
    This is necessary to ensure that the proseg command can be found when running the script.
    Returns:
        bool: True if the proseg version is 3.x, False otherwise.
    """
    # Need to add the rust path to the PATH variable
    if rust_bin_path is not None:
        rust_path = pathlib.Path(rust_bin_path)
    else:
        rust_path = pathlib.Path(os.getenv("RUST_BIN_PATH", ""))
    if not rust_path.exists():
        raise ValueError(f"RUST_BIN_PATH {rust_path} does not exist. Please set the correct path.")
    os.environ["PATH"] += os.pathsep + str(rust_path)

    # Test that proseg can be run:
    ret = subprocess.run(["proseg", "--help"], capture_output=True, check=True)
    if ret.returncode != 0:
        raise ValueError(
            "Proseg failed to run. Check the installation and the PATH variable."
        )
    proseg_v3 = _determine_proseg_version()[0] == '3'
    logger.info(f"Proseg version 3.x detected: {proseg_v3}")

    proseg_out = ret.stdout.decode()
    proseg_options = proseg_out.split("Options:")[-1]

    # Extract the options from the proseg output
    global proseg_args
    escaped_start = re.escape("--")
    escaped_end = re.escape("\n")
    pattern = f"{escaped_start}(.*?){escaped_end}"
    proseg_args = re.findall(pattern, proseg_options.strip())
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    proseg_args = [ansi_escape.sub("", arg) for arg in proseg_args]
    proseg_args = [arg.split()[0] for arg in proseg_args]

    logger.info("proseg binary added to PATH successfully.")

    return proseg_v3


@DeprecationWarning
def _execute_cli_proseg_v2(root_dir: str, output_dir: str, region: str, **proseg_params):
    """
    Run the proseg command with the prespecified parameters using subprocess.
    """

    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \n
        --output-path {output_dir}/{region} \n
        --output-expected-counts expected-counts.csv.gz \n
        --output-cell-metadata cell-metadata.csv.gz \n
        --output-transcript-metadata transcript-metadata.csv.gz \n
        --output-cell-polygons cell-polygons.geojson.gz \n
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \n
        --output-union-cell-polygons union-cell-polygons.geojson.gz \n
        """
    for _key, _val in proseg_params.items():
        if _key not in proseg_args:
            logger.info("skipping unknown proseg parameter: %s" % _key)
            continue
        if isinstance(_val, bool):
            if _val:
                logger.info(f"Adding boolean proseg parameter: --{_key}")
                proseg_run_cmd += f"--{_key} \n"
        else:
            logger.info(f"Adding proseg parameter: --{_key} {_val}")
            proseg_run_cmd += f"--{_key} {_val} \n"

    # proseg_run_cmd += f"--nthreads {os.cpu_count()} \\\n"
    logger.info(f"Running proseg with command:\n{proseg_run_cmd}")

    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret


def _execute_cli_proseg_v3(root_dir: str, output_dir: str, region: str, **proseg_params):
    """
    Run the proseg command for the 3.0.1 version of proseg
    """
    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \n
        --output-path {output_dir}/{region} \n
        --output-spatialdata proseg_outputs.zarr \n
        --output-cell-polygons cell-polygons.geojson.gz \n
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \n
        --output-cell-metadata cell-metadata.csv.gz \n
        """
    for _key, _val in proseg_params.items():
        if _key not in proseg_args:
            logger.info("skipping unknown proseg parameter: %s" % _key)
            continue
        if isinstance(_val, bool):
            if _val:
                logger.info(f"Adding boolean proseg parameter: --{_key}")
                proseg_run_cmd += f"--{_key} \n"
        else:
            logger.info(f"Adding proseg parameter: --{_key} {_val}")
            proseg_run_cmd += f"--{_key} {_val} \n"

    # proseg_run_cmd += f"--nthreads {os.cpu_count()} \\\n"
    logger.info(f"Running proseg with command:\n{proseg_run_cmd}")

    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret



def run_proseg(root_dir: str, output_dir: str, region: str, rust_bin_path=None, **proseg_params):
    """
    Run the proseg command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    """

    for key, value in proseg_params.items():
        logger.info(f"Proseg parameter: {key} = {value}")

    # Ensure the output directory exists
    pathlib.Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)

    # Add proseg path to environment
    v3_flag = _add_proseg_binary(rust_bin_path=rust_bin_path)

    # Prepare the command to run proseg
    if v3_flag: 
        ret = _execute_cli_proseg_v3(root_dir, output_dir, region, **proseg_params)
    else: # run the old version of proseg
        ret = _execute_cli_proseg_v2(root_dir, output_dir, region, **proseg_params)

    print(ret.returncode)
    print(ret.stdout.decode("utf-8"))


def _proseg_3d_to_vpt_parquet(shapes_3d: gpd.GeoDataFrame, micron_per_z: float, output_path: str):
    shapes_3d = shapes_3d.rename_geometry("Geometry")
    shapes_3d['z'] = shapes_3d['layer']
    shapes_3d['ID'] = shapes_3d['cell']
    shapes_3d['ZIndex'] = shapes_3d['layer']
    shapes_3d['Zlevel'] = shapes_3d['ZIndex'] * micron_per_z
    shapes_3d['EntityID'] = pd.factorize(shapes_3d['ID'])[0] + 1
    shapes_3d.to_parquet(output_path, index=False)

def add_signals_meta_to_proseg(
    root_dir: str | pathlib.Path,
    proseg_dir: str | pathlib.Path,
    reg_name: str,
    shapes_filename: str = "cell-polygons-layers.geojson.gz",
    cell_meta_filename: str = "cell-metadata.csv.gz",
    micron_per_z: float = 1.5,
    vpt_bin_path: str | pathlib.Path | None = None,
    **kwargs
):
    """
    Add the sum signals metadata to the proseg cell metadata output.
    This is necessary for downstream analysis and visualization in VPT.
    """
    # from spida.S.segmentation.vpt_utils import read_micron_to_mosaic_transform, transform_geoms, update_geometry
    # import geopandas as gpd
    # import pandas as pd
    from spida.S.segmentation.vpt import _add_vpt_binary, _cli_sum_signals    
    _add_vpt_binary(vpt_bin_path=vpt_bin_path)

    if isinstance(proseg_dir, str):
        proseg_dir = pathlib.Path(proseg_dir)

    # Get paths and make sure they exist
    shapes_3d_path = proseg_dir / reg_name / shapes_filename
    cell_meta_path = proseg_dir / reg_name / cell_meta_filename
    assert shapes_3d_path.exists(), f"Shapes file not found at {shapes_3d_path}"
    assert cell_meta_path.exists(), f"Cell metadata file not found at {cell_meta_path}"
    # Read in the shapes and cell metadata
    with gzip.open(shapes_3d_path, "rb") as f:
        shapes_3d = gpd.read_file(f)
    logger.info(f"Read shapes from {shapes_3d_path} with shape {shapes_3d.shape}")
    with gzip.open(cell_meta_path, "rt") as f:
        cell_meta = pd.read_csv(f)
    logger.info(f"Read cell metadata from {cell_meta_path} with shape {cell_meta.shape}")

    # Convert proseg 3D shapes to vpt parquet format
    _proseg_3d_to_vpt_parquet(
        shapes_3d=shapes_3d,
        micron_per_z=micron_per_z,
        output_path= proseg_dir / reg_name / "proseg_polygons.parquet"
    )
    logger.info(f"Converted proseg 3D shapes to vpt parquet format at {proseg_dir / reg_name / 'proseg_polygons.parquet'}")
    
    # Sum signals 
    signals = _cli_sum_signals(
       root_dir = root_dir,
       output_dir = proseg_dir,
       region = reg_name,
       input_boundaries = "proseg_polygons.parquet",
       output_signals = "sum_signals.csv",
    )

    # Merge the signals with the cell metadata and save
    signals.index.name = 'EntityID'
    if (cell_meta.index.astype(int) == 0).any(): 
        signals.index = signals.index.astype(int) - 1
    cell_meta = cell_meta.merge(signals, left_on="cell", right_index=True).copy()
    cell_meta.to_csv(proseg_dir / reg_name / "cell_metadata_with_signals.csv", index=False)
    logger.info(f"Merged signals with cell metadata and saved to {proseg_dir / reg_name / 'cell_metadata_with_signals.csv'}")

@DeprecationWarning
def align_proseg_transcripts(
    exp_name: str,
    reg_name: str,
    seed_prefix_name: str,
    prefix_name: str,
    x: str = "x",
    y: str = "y",
    z: str = "global_z",
    cell_column: str = "cell_id",
    barcode_column: str = "barcode_id",
    gene_column: str = "gene",
    fov_column: str = "fov",
    cell_missing: str = "-1",
    min_jaccard: float = 0.4,
    min_prob: float = 0.5,
    filter_blank: bool = False,
    merged_transcripts_fname: str = "merged_transcripts.csv",
    merged_cell_by_gene_fname: str = "merged_cell_by_gene.csv",
    merged_cell_polygons_fname: str = "merged_cell_polygons.parquet",
    save_dir: str | pathlib.Path = None,
):
    """
    Align Proseg Transcripts to Seed Transcripts in order to filter out erroneous Proseg Segmentation Results.

    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    seed_prefix_name (str): Prefix for the seed transcripts.
    prefix_name (str): Prefix for the proseg transcripts.
    x, y, z: Column names for spatial coordinates.
    cell_column: Column name for cell identifiers in seed transcripts.
    barcode_column: Column name for barcode identifiers in seed transcripts.
    gene_column: Column name for gene identifiers in seed transcripts.
    fov_column: Column name for field of view identifiers in seed transcripts.
    cell_missing: Value indicating missing cells in the dataset.
    min_jaccard: Minimum Jaccard index to consider a match valid.
    min_prob: Minimum probability threshold for matching.
    filter_blank: Whether to filter out blank transcripts from the final cell by gene counts.

    merged_transcripts_fname: Filename for saving the merged transcripts.
    merged_cell_by_gene_fname: Filename for saving the cell by gene counts.
    merged_cell_polygons_fname: Filename for saving the merged cell polygons.

    Returns:
    None
    """
    from spida._utilities import _gen_keys
    from spida._constants import POINTS_KEY, SHAPES_KEY
    import spatialdata as sd
    import numpy as np
    import pandas as pd

    # KEYS
    SEED_KEYS = _gen_keys(seed_prefix_name, exp_name, reg_name)
    PROSEG_KEYS = _gen_keys(prefix_name, exp_name, reg_name)

    zarr_store = os.getenv("ZARR_STORAGE_PATH")
    zarr_path = f"{zarr_store}/{exp_name}/{reg_name}"

    # Read in the sdata object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sdata = sd.read_zarr(zarr_path)

    # Load the seed and the proseg data
    seed_transcripts = sdata[SEED_KEYS[POINTS_KEY]].copy().compute()
    seed_polygons = sdata[SEED_KEYS[SHAPES_KEY]].copy()
    proseg_transcripts = sdata[PROSEG_KEYS[POINTS_KEY]].copy().compute()
    proseg_polygons = sdata[PROSEG_KEYS[SHAPES_KEY]].copy()

    # Format seed transcripts
    transcripts = seed_transcripts[
        [x, y, z, cell_column, gene_column, barcode_column, fov_column]
    ]
    transcripts.loc[:, cell_column] = (
        transcripts.loc[:, cell_column].astype(str).replace(cell_missing, np.nan)
    )
    transcripts = transcripts.round(2)

    # Format proseg transcripts
    proseg = proseg_transcripts.query("probability > 0.5 & background == 0")[
        ["observed_x", "observed_y", "observed_z", "assignment", "fov"]
    ].rename({"observed_x": "x", "observed_y": "y", "observed_z": "z"}, axis=1)
    proseg = proseg.round(2)

    logger.info(
        f"{proseg.shape[0] / proseg_transcripts.shape[0] * 100:.2f}% of Proseg Transcripts with Probability > 0.5 and not Background"
    )
    logger.info(
        f"unique cellpose cells: {transcripts[cell_column].nunique()}; unique proseg cells {proseg['assignment'].nunique()}"
    )

    # Merge Transcripts on X and Y coordinates (drops transcripts in proseg without a match in seed transcripts)
    merged = transcripts.merge(proseg, on=[x, y], how="left", suffixes=("", "_proseg"))

    logger.info(
        f"num of aligned proseg transcripts: {(~pd.isna(merged['assignment'])).sum()}"
    )
    logger.info(
        f"percent of aligned proseg transcripts: {(~pd.isna(merged['assignment'])).sum() / proseg.shape[0] * 100:.2f}"
    )

    # Each row is a an overlapping seed, proseg cell combination.
    # Count is the number of transcripts that overlap between the two cells.
    mapping = (
        merged.query(f"assignment.notna() & {cell_column} != {cell_missing}")
        .groupby([cell_column, "assignment"])
        .size()
        .reset_index(name="count")
    )
    # Misleading due to the filtering of transcripts in Proseg cells
    logger.info(
        f"{mapping['count'].sum() / transcripts.shape[0] * 100:.3f}% Seed Transcripts within Aligned Cells "
    )
    logger.info(
        f"{mapping['count'].sum() / proseg.shape[0] * 100:.3f}% Proseg Transcripts within Aligned Cells "
    )

    # Generating the Jaccard Index for each cell pair -> The similarity measure used between the cells
    mapping["cell_count"] = mapping.groupby(cell_column)["count"].transform("sum")
    mapping["proseg_count"] = mapping.groupby("assignment")["count"].transform("sum")
    mapping["jaccard"] = mapping["count"] / (
        mapping["cell_count"] + mapping["proseg_count"] - mapping["count"]
    )

    # Only keeping the seed, proseg cell pairing with the highst jaccard index
    mapping = (
        mapping.sort_values("jaccard", ascending=False)
        .groupby("assignment")
        .first()
        .reset_index()
        .groupby(cell_column)
        .first()
        .reset_index()
    )

    # Filter out cell combinations with a low jaccard index (Proseg not aligned with Seed)
    mapping = mapping.query("jaccard >= @min_jaccard").copy()

    logger.info("Number of Mapped Seed Cells: %i" % mapping[cell_column].nunique())
    logger.info("Number of Assignments: %i" % mapping["assignment"].nunique())
    logger.info("Should Match")

    # Merging the mapping between seed and proseg cells with the overall merged transcript dataframe
    mapping.rename(columns={cell_column: "mapped_cell"}, inplace=True)
    merged = merged.merge(
        mapping[["assignment", "mapped_cell"]], on="assignment", how="left"
    )

    logger.info(
        "Number of Unmapped Seed Cells: %i; Number of Removed Proseg Cells: %i"
        % (
            merged.groupby(cell_column)["mapped_cell"]
            .apply(lambda x: x.isna().all())
            .sum(),
            merged.query("assignment.notna() & mapped_cell.isna()")[
                "assignment"
            ].nunique(),
        )
    )

    # This essentially combines the mapping and default in a way that all seed unmapped cells
    # remain in the dataset, but the seed mapped to proseg cells are merged with the proseg cells
    merged[cell_column] = merged["mapped_cell"].fillna(merged[cell_column])

    logger.info(
        "%.2f%% of transcripts in Final cells, %.2f%% of transcripts in only joint cells, %.2f%% in proseg cells"
        % (
            100 * (1 - merged[cell_column].isna().mean()),
            100 * (1 - merged["mapped_cell"].isna().mean()),
            100 * (1 - merged["assignment"].isna().mean()),
        )
    )

    proseg_to_seed_dict = mapping.set_index("assignment")["mapped_cell"].to_dict()
    _final_seeds = merged[cell_column].dropna().unique()
    logger.info("Number of Final Seeds: %i" % len(_final_seeds))
    _final_proseg = mapping["assignment"].unique()
    logger.info("Number of Final Proseg Cells: %i" % len(_final_proseg))

    # Get counts - Cell x Gene
    counts = merged.groupby(cell_column)[gene_column].value_counts().unstack().fillna(0)
    counts = counts.astype(int)
    counts.index = counts.index.astype(int).values
    counts.columns.name = None
    if filter_blank:
        logger.info("Filtering Blank Transcripts from cell by gene")
        counts = counts.loc[:, ~counts.columns.str.contains("Blank")]

    # Mapping the proseg polygons to the matching seed IDs
    proseg_polygons = proseg_polygons.loc[_final_proseg.astype(int).astype(str)]
    proseg_polygons["seed_id"] = (
        proseg_polygons["cell"].map(proseg_to_seed_dict).astype(str)
    )

    # Merging the two geometries together
    merge_polygons = seed_polygons.loc[_final_seeds]
    merge_polygons = merge_polygons.merge(
        proseg_polygons,
        how="left",
        left_index=True,
        right_on="seed_id",
        suffixes=("_seed", "_proseg"),
    )
    merge_polygons["is_proseg"] = ~merge_polygons["geometry_proseg"].isna()
    merge_polygons["geometry"] = merge_polygons["geometry_proseg"].combine_first(
        merge_polygons["geometry_seed"]
    )
    drop_cols = [
        "geometry_seed",
        "geometry_proseg",
        "seed_id",
        "cell",
        "area",
        "area_seed",
        "area_proseg",
    ]
    merge_polygons = (
        merge_polygons.set_index("EntityID")
        .drop(drop_cols, axis=1, errors="ignore")
        .set_geometry("geometry")
    )

    # Merging the transcripts
    merge_transcripts = merged.drop(
        ["assignment", "fov_proseg", "mapped_cell", "z"], axis=1
    )
    merge_transcripts[cell_column] = merge_transcripts[cell_column].fillna(cell_missing)
    merge_transcripts = merge_transcripts.rename(
        columns={x: "global_x", y: "global_y", z: "global_z"}
    )

    ### Saving the Outputs:
    logger.info(f"Saving Outputs to {save_dir}")
    counts.to_csv(save_dir / merged_cell_by_gene_fname)
    merge_polygons.to_file(save_dir / merged_cell_polygons_fname, driver="GeoJSON")
    merge_transcripts.to_csv(save_dir / merged_transcripts_fname)

@DeprecationWarning
def filt_to_ids(meta_path, geom_path, tz_path, cbg_path):
    """
    Filter metadata and geometry to only include cells with IDs in the provided list.
    Additionally rename the index of the merged_cell_by_gene
    """
    import pandas as pd
    import geopandas as gpd

    cbg = pd.read_csv(cbg_path, index_col=0)
    transcripts = pd.read_csv(tz_path, index_col=0)
    meta = pd.read_csv(meta_path)
    geom = gpd.read_parquet(geom_path)
    ids = transcripts["cell_id"].unique()
    meta_filt = meta[meta["EntityID"].isin(ids)]
    geom_filt = geom[geom["EntityID"].isin(ids)]
    cbg.index = cbg.index.astype(str)
    cbg.index.name = "cell"

    meta_filt.to_csv(meta_path, index=False)
    geom_filt.to_parquet(geom_path, index=False)
    cbg.to_csv(cbg_path)
    logger.info(
        f"Filtered metadata and geometry to {len(meta_filt)} cells with IDs in {len(ids)} unique IDs."
    )
    logger.info(
        f"Metadata shape: {meta_filt.shape}, Geometry shape: {geom_filt.shape}, CBG shape: {cbg.shape}"
    )
