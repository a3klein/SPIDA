# sdata Keys table keys
SHAPES_KEY = "shapes_key"
POINTS_KEY = "points_key"
TABLE_KEY = "table_key"
IMAGE_KEY = "image_key"

# For annotation with table on shapes
REGION_KEY = "cells_region"

# Constants for column names
X_LOC = "X_GLOBAL"
Y_LOC = "Y_GLOBAL"
Z_LOC = "Z_GLOBAL"
CELL_ID = "CELL_ID"
GENE_NAME = "GENE_NAME"
GENE_ID = "GENE_ID"
FOV_COL = "FOV"

CELL_X = "CENTER_X"
CELL_Y = "CENTER_Y"
CELL_FOV = "fov"
CELL_VOLUME = "volume"
GEOMETRY_COL = "geometry"
DAPI_COL = "DAPI_high_pass"

DEFAULT_PRESET = {
    "tz_cell_id": "cell",
    "meta_map": {
        "EntityID": CELL_ID,
        "fov": CELL_FOV,
        "center_x": CELL_X,
        "center_y": CELL_Y,
        "volume": CELL_VOLUME,
    },
    "tz_col_mapper": {
        "global_x": X_LOC,
        "global_y": Y_LOC,
        "global_z": Z_LOC,
        "cell_id": CELL_ID,
        "gene": GENE_NAME,
        "transcript_id": GENE_ID,
        "fov": FOV_COL,
    },
    "geom_depth_col": "ZIndex",
    "METADATA_CELL_KEY": "EntityID",
}

PROSEG_PRESET = {
    "tz_cell_id": None,
    "meta_map": {
        "cell": CELL_ID,
        "fov": CELL_FOV,
        "centroid_x": CELL_X,
        "centroid_y": CELL_Y,
        "volume": CELL_VOLUME,
    },
    "meta_cols": {
        "cell_col": "cell",
        "x_col": "centroid_x",
        "y_col": "centroid_y",
        "fov_col": "fov",
        "volume_col": "volume",
    },
    "tz_col_mapper": {
        "x": X_LOC,
        "y": Y_LOC,
        "z": Z_LOC,
        "assignment": CELL_ID,
        "gene": GENE_NAME,
        "transcript_id": GENE_ID,
        "fov": FOV_COL,
    },
    "geom_depth_col": "layer",
    "tz_col": {"x_col": "x", "y_col": "y", "gene_col": "gene"},
}


def rename_exp_salk(x):
    exp_n = x.split("_")[1]
    temp = exp_n.split("-")[-2]
    if len(temp) == 1:
        temp = exp_n.split("-")[-3]  # noqa: E701
    return temp


def rename_reg_salk(x):
    reg_n = x.split("_")[1]
    temp = reg_n.split("-")[0]
    if len(temp) == 3:
        temp = "".join(reg_n.split("-")[0:2])  # noqa: E701
    return temp


def rename_exp_ucsd(x):
    exp_n = x.split("_")[1].split("BICAN")[-1]
    return exp_n


def rename_reg_ucsd(x):
    reg_n = x.split("_")[1].split("Q0")[0]
    return reg_n

ren_to_exp_map = {
    "Ren35": "CAH",
    "Ren36": "PU",
    "Ren37": "GP",
    "Ren38": "CAB",
    "Ren39": "MGM1",
    "Ren40": "NAC",
    "Ren41": "SUBTH",
    "Ren42": "CAT",
}