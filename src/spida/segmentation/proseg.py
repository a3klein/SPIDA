import os
from dotenv import load_dotenv # type: ignore
load_dotenv()

import warnings
import pathlib
import subprocess

def _add_proseg_binary(): 

    """
    Add the path to the proseg binary to the PATH environment variable.
    This is necessary to ensure that the proseg command can be found when running the script.
    """
    # Need to add the rust path to the PATH variable
    os.environ['PATH'] += os.pathsep + os.getenv("RUST_BIN_PATH")

    #Test that proseg can be run: 
    ret = subprocess.run(["proseg", "-h"], capture_output=True, check=True)
    if ret.returncode != 0: 
        raise ValueError("Proseg failed to run. Check the installation and the PATH variable.")
    
    print("proseg binary added to PATH successfully.")

def _execute_cli_proseg(root_dir:str, output_dir:str, region:str, **proseg_params): 
    """
    Run the proseg command with the prespecified parameters using subprocess. 
    """

    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \
        --output-path {output_dir}/{region} \
        --output-expected-counts expected-counts.csv.gz \
        --output-cell-metadata cell-metadata.csv.gz \
        --output-transcript-metadata transcript-metadata.csv.gz \
        --output-cell-polygons cell-polygons.geojson.gz \
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \
        --output-union-cell-polygons union-cell-polygons.geojson.gz \
        --voxel-layers 4 \
        """
    for _key, _val in proseg_params.items():
        proseg_run_cmd += f"--{_key} {_val} \\\n"

    # proseg_run_cmd += f"--nthreads {os.cpu_count()} \\\n"
    print(f"Running proseg with command:\n{proseg_run_cmd}")

    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret
    

def run_proseg(root_dir:str, output_dir:str, region:str, **proseg_params):
    """
    Run the proseg command to process detected transcripts in a specified region.

    Parameters:
    root_dir (str): The root directory containing the detected transcripts.
    output_dir (str): The directory where the output files will be saved.
    region (str): The name of the region to process.
    """
    
    # Ensure the output directory exists
    pathlib.Path(f"{output_dir}/{region}").mkdir(parents=True, exist_ok=True)
    
    # Add proseg path to environment
    _add_proseg_binary()

    # Prepare the command to run proseg
    ret = _execute_cli_proseg(root_dir, output_dir, region, **proseg_params)

    print(ret.returncode)
    print(ret.stdout.decode('utf-8'))


def align_proseg_transcripts(
        exp_name: str, 
        reg_name: str, 
        seed_prefix_name:str, 
        prefix_name:str,
        x : str = "x",
        y : str = "y",
        z : str = "global_z",
        cell_column : str = "cell_id",
        barcode_column : str = "barcode_id",
        gene_column : str = "gene",
        fov_column : str = "fov",
        cell_missing : str = '-1',
        min_jaccard : float = 0.4,
        min_prob : float = 0.5,
        filter_blank : bool = False,
        merged_transcripts_fname:str = "merged_transcripts.csv",
        merged_cell_by_gene_fname:str = "merged_cell_by_gene.csv",
        merged_cell_polygons_fname:str = "merged_cell_polygons.parquet",
        save_dir: str | pathlib.Path = None
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
    transcripts = seed_transcripts[[x, y, z, cell_column, gene_column, barcode_column, fov_column]]
    transcripts.loc[:, cell_column] = transcripts.loc[:, cell_column].astype(str).replace(cell_missing, np.nan)
    transcripts = transcripts.round(2)

    # Format proseg transcripts
    proseg = proseg_transcripts.query("probability > 0.5 & background == 0")[['observed_x', 'observed_y', 'observed_z', 'assignment', 'fov']].rename(
        {'observed_x': 'x', 'observed_y': 'y', 'observed_z': 'z'}, axis=1
    )
    proseg = proseg.round(2)

    print(proseg.shape[0] / proseg_transcripts.shape[0] * 100, "% of Proseg Transcripts with Probability > 0.5 and not Background")
    print("unique cellpose cells: %i; unique proseg cells %i" %
       (transcripts[cell_column].nunique(), proseg['assignment'].nunique()))
    
    # Merge Transcripts on X and Y coordinates (drops transcripts in proseg without a match in seed transcripts)
    merged = transcripts.merge(proseg, on=[x, y], how="left", suffixes=("", "_proseg"))

    print("num of aligned proseg transcripts", (~pd.isna(merged['assignment'])).sum())
    print("percent of aligned proseg transcripts", (~pd.isna(merged['assignment'])).sum() / proseg.shape[0] * 100)

    # Each row is a an overlapping seed, proseg cell combination. 
    # Count is the number of transcripts that overlap between the two cells.
    mapping = (merged.query(f"assignment.notna() & {cell_column} != {cell_missing}").
            groupby([cell_column, "assignment"])
            .size()
            .reset_index(name="count")
    )
    # Misleading due to the filtering of transcripts in Proseg cells
    print("%.3f%% Seed Transcripts within Aligned Cells " % (mapping['count'].sum() / transcripts.shape[0] * 100))
    print("%.3f%% Proseg Transcripts within Aligned Cells " % (mapping['count'].sum() / proseg.shape[0] * 100))

    # Generating the Jaccard Index for each cell pair -> The similarity measure used between the cells
    mapping['cell_count'] = mapping.groupby(cell_column)['count'].transform('sum')
    mapping['proseg_count'] = mapping.groupby("assignment")['count'].transform('sum')
    mapping["jaccard"] = mapping["count"] / (mapping["cell_count"] + mapping["proseg_count"] - mapping["count"])

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
    
    print("Number of Mapped Seed Cells: %i" % mapping[cell_column].nunique())
    print("Number of Assignments: %i" % mapping['assignment'].nunique())
    print("Should Match")

    # Merging the mapping between seed and proseg cells with the overall merged transcript dataframe
    mapping.rename(columns={cell_column: "mapped_cell"}, inplace=True)
    merged = merged.merge(mapping[['assignment', 'mapped_cell']], on="assignment", how="left")

    print("Number of Unmapped Seed Cells: %i; Number of Removed Proseg Cells: %i" % 
      (merged.groupby(cell_column)['mapped_cell'].apply(lambda x : x.isna().all()).sum(),
        merged.query("assignment.notna() & mapped_cell.isna()")['assignment'].nunique()))

    # This essentially combines the mapping and default in a way that all seed unmapped cells 
    # remain in the dataset, but the seed mapped to proseg cells are merged with the proseg cells
    merged[cell_column] = merged["mapped_cell"].fillna(merged[cell_column])

    print("%.2f%% of transcripts in Final cells, %.2f%% of transcripts in only joint cells, %.2f%% in proseg cells" %
       (100 * ( 1 - merged[cell_column].isna().mean()),
        100 * (1 - merged['mapped_cell'].isna().mean()),
        100 * (1 - merged['assignment'].isna().mean())))
    
    proseg_to_seed_dict = mapping.set_index("assignment")["mapped_cell"].to_dict()
    _final_seeds = merged[cell_column].dropna().unique()
    print("Number of Final Seeds: %i" % len(_final_seeds))
    _final_proseg = mapping['assignment'].unique()
    print("Number of Final Proseg Cells: %i" % len(_final_proseg))

    # Get counts - Cell x Gene
    counts = merged.groupby(cell_column)[gene_column].value_counts().unstack().fillna(0)
    counts = counts.astype(int)
    counts.index = counts.index.astype(int).values
    counts.columns.name = None
    if filter_blank: 
        print("Filtering Blank Transcripts from cell by gene")
        counts = counts.loc[:, ~counts.columns.str.contains("Blank")]

    # Mapping the proseg polygons to the matching seed IDs
    proseg_polygons = proseg_polygons.loc[_final_proseg.astype(int).astype(str)]
    proseg_polygons["seed_id"] = proseg_polygons['cell'].map(proseg_to_seed_dict).astype(str)

    # Merging the two geometries together
    merge_polygons = seed_polygons.loc[_final_seeds]
    merge_polygons = merge_polygons.merge(proseg_polygons, how="left", left_index=True, right_on="seed_id", suffixes=("_seed", "_proseg"))
    merge_polygons['is_proseg'] = ~merge_polygons['geometry_proseg'].isna()
    merge_polygons['geometry'] = merge_polygons['geometry_proseg'].combine_first(merge_polygons['geometry_seed'])
    drop_cols = ["geometry_seed", "geometry_proseg", "seed_id", "cell", 'area', 'area_seed', 'area_proseg']
    merge_polygons = (merge_polygons
                      .set_index("EntityID")
                      .drop(drop_cols, axis=1, errors='ignore')
                      .set_geometry("geometry") )

    # Merging the transcripts 
    merge_transcripts = merged.drop(["assignment", "fov_proseg", "mapped_cell", "z"], axis=1)
    merge_transcripts[cell_column] = merge_transcripts[cell_column].fillna(cell_missing)

    ### Saving the Outputs: 
    print(f"Saving Outputs to {save_dir}")
    counts.to_csv(save_dir / merged_cell_by_gene_fname)
    merge_polygons.to_file(save_dir / merged_cell_polygons_fname, driver="GeoJSON")
    merge_transcripts.to_csv(save_dir / merged_transcripts_fname)
