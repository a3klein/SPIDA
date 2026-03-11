import os
from pathlib import Path
import warnings
import logging

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
# import spatialdata as sd
import geopandas as gpd
import shapely as shp

from spida.utilities.tiling import create_hexagonal_grid, hex_area_from_size


logger = logging.getLogger(__package__)

def _ensure_geodataframe(points, x_col: str = "x", y_col: str = "y") -> gpd.GeoDataFrame:
    if hasattr(points, "compute"):
        points = points.compute()
    if isinstance(points, gpd.GeoDataFrame):
        return points
    if isinstance(points, pd.DataFrame):
        return gpd.GeoDataFrame(points, geometry=gpd.points_from_xy(points[x_col], points[y_col]))
    raise TypeError("points must be a GeoDataFrame, DataFrame, or dask DataFrame")

def _hex_table_key(hex_size: float, hex_overlap: float) -> str:
    return f"adata_hex_s{hex_size}_o{hex_overlap}"

def _add_hex_metadata(adata_hex: ad.AnnData, hex_size: float, hex_overlap: float) -> None:
    adata_hex.uns.setdefault("hexgrid", {})
    adata_hex.uns["hexgrid"]["hex_size"] = hex_size
    adata_hex.uns["hexgrid"]["hex_overlap"] = hex_overlap
    adata_hex.uns["hexgrid"]["hex_area_um2"] = hex_area_from_size(hex_size)

def _apply_hex_filter(
    adata_hex: ad.AnnData,
    min_transcripts: int | None = None,
    min_density: float | None = None,
) -> ad.AnnData:
    if min_density is not None:
        mask = adata_hex.obs["density"] >= min_density
    elif min_transcripts is not None:
        mask = adata_hex.obs["tz_count"] >= min_transcripts
    else:
        return adata_hex
    return adata_hex[mask].copy()

def _serialize_obs_geometry(obs_df: pd.DataFrame, geometry_col: str = "geometry") -> pd.DataFrame:
    """
    Convert geometry objects to a serializable representation for AnnData .obs.
    """
    if isinstance(obs_df, gpd.GeoDataFrame):
        out = pd.DataFrame(obs_df.drop(columns=[geometry_col], errors="ignore")).copy()
        out["geometry_wkt"] = obs_df.geometry.to_wkt()
        return out

    out = obs_df.copy()
    if geometry_col in out.columns:
        # shapely objects in object dtype -> WKT
        out["geometry_wkt"] = out[geometry_col].apply(
            lambda g: g.wkt if g is not None and hasattr(g, "wkt") else None
        )
        out = out.drop(columns=[geometry_col], errors="ignore")
    return out

def _obs_to_grid_geodf(obs_df: pd.DataFrame) -> gpd.GeoDataFrame:
    """
    Reconstruct a GeoDataFrame from AnnData .obs when geometry is stored as WKT.
    """
    if "geometry_wkt" in obs_df.columns:
        geom = gpd.GeoSeries.from_wkt(obs_df["geometry_wkt"])
        base = obs_df.drop(columns=["geometry_wkt"]).copy()
        return gpd.GeoDataFrame(base, geometry=geom)
    if "geometry" in obs_df.columns:
        return gpd.GeoDataFrame(obs_df.copy(), geometry="geometry")
    return gpd.GeoDataFrame(obs_df.copy())

def overlay_hexes(
    points : gpd.GeoDataFrame | pd.DataFrame,
    x_col : str = "x",
    y_col : str = "y",
    hex_size : float = 30,
    hex_overlap : float = 0,
    gene_col: str = "gene",
    return_grid: bool = False,
): 
    """
    Overlay hexagonal grid on points and return an AnnData object with hexagons as observations and genes as variables.
    Parameters
    ----------
    points : gpd.GeoDataFrame or pd.DataFrame
        DataFrame containing the point coordinates.
    x_col : str, optional
        Name of the column containing x coordinates, by default 'x'.
    y_col : str, optional
        Name of the column containing y coordinates, by default 'y'.
    hex_size : float, optional
        Size of the hexagons (distance from center to any vertex), by default 30.
    hex_overlap : float, optional
        Amount of overlap between hexagons (0 to 1), by default 0.
    gene_col : str, optional
        Name of the column in points GeoDataFrame that contains gene information, by default 'gene'.
    return_grid : bool, optional
        If True, return the hexagonal grid as a GeoDataFrame.
    Returns
    -------
    ad.AnnData
        AnnData object with hexagons as observations and genes as variables, containing the count of points in each hexagon for each gene.
    pd.DataFrame
        DataFrame containing the hexagonal grid geometry and metadata.
    """

    points = _ensure_geodataframe(points, x_col=x_col, y_col=y_col)

    # Create hexagonal grid
    total_bounds = points.total_bounds  # [minx, miny, maxx, maxy]
    grid = create_hexagonal_grid(total_bounds, hex_size, hex_overlap)
    grid["hex_id"] = grid.index.astype(str)
    grid = grid.set_index("hex_id")
    
    # Spatial join to count points in each hexagon
    joint_grid = gpd.sjoin(grid, points, how="inner", predicate="contains")
    logger.info(f"Created hexagonal grid with {len(grid)} hexes.")

    # Aggregate counts of points in each hexagon for each gene
    temp_count = joint_grid[gene_col].to_frame().reset_index()
    temp_count["count"] = 1
    temp_count = temp_count.groupby(["hex_id", gene_col], observed=True)["count"].sum().reset_index()
    X = temp_count.pivot(index="hex_id", columns=gene_col, values="count").fillna(0)
    del joint_grid
    del temp_count

    # TODO: Do the removal of blank with a regex type string
    df_obs = grid.loc[grid.index.isin(X.index)].copy()
    X.index = X.index.astype(int)
    X = X.sort_index()
    X = X.loc[:, ~X.columns.str.contains("Blank-")]

    obs_serializable = _serialize_obs_geometry(df_obs)

    adata_hex = ad.AnnData(
        obs=obs_serializable,
        X = X, 
        var = pd.DataFrame(index=X.columns)
    )
    adata_hex.obs["tz_count"] = np.asarray(adata_hex.X.sum(axis=1)).ravel()
    adata_hex.obs["n_genes"] = np.asarray((adata_hex.X > 0).sum(axis=1)).ravel()
    adata_hex.var["cell_count"] = np.asarray(adata_hex.X.sum(axis=0)).ravel()

    hex_area = hex_area_from_size(hex_size)
    adata_hex.obs["density"] = adata_hex.obs["tz_count"] / hex_area
    _add_hex_metadata(adata_hex, hex_size=hex_size, hex_overlap=hex_overlap)

    if return_grid:
        return adata_hex, _obs_to_grid_geodf(adata_hex.obs)
    return adata_hex

def get_or_create_hex_adata(
    sdata,
    points_key: str,
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    table_key: str | None = None,
    force_recompute: bool = False,
):
    key = table_key or _hex_table_key(hex_size, hex_overlap)
    if not force_recompute and hasattr(sdata, "tables") and key in sdata.tables:
        return sdata.tables[key], key

    points = sdata[points_key]
    adata_hex = overlay_hexes(
        points,
        x_col=x_col,
        y_col=y_col,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        return_grid=False,
    )

    _add_hex_metadata(adata_hex, hex_size=hex_size, hex_overlap=hex_overlap)
    sdata.tables[key] = adata_hex
    sdata.write_element(key)
    return adata_hex, key

def transcript_qc(
    adata_hex : ad.AnnData,
    min_transcripts : int | None = None,
    min_density : float | None = None,
): 
    if "density" not in adata_hex.obs.columns:
        hex_area = adata_hex.uns.get("hexgrid", {}).get("hex_area_um2", None)
        if hex_area is None:
            raise ValueError("hex_area_um2 not found in adata_hex.uns['hexgrid']")
        adata_hex = adata_hex.copy()
        adata_hex.obs["density"] = adata_hex.obs["tz_count"] / hex_area

    adata_filt = _apply_hex_filter(
        adata_hex,
        min_transcripts=min_transcripts,
        min_density=min_density,
    )

    grid = _obs_to_grid_geodf(adata_hex.obs.copy())
    grid["filtered"] = True
    grid.loc[adata_filt.obs.index, "filtered"] = False
    return adata_filt, grid

def cluster_hexes(
    adata_hex : ad.AnnData,
    leiden_resolution : float = 1.0,
    min_cells : int = 10,
    min_genes : int = 100,
    n_top_genes : int = 200,    
    pca_comps : int = 50,
): 
    """
    Cluster precomputed hex-binned AnnData.
    """
    adata_hex = adata_hex.copy()

    orig_hexes = adata_hex.shape[0]
    orig_genes = adata_hex.shape[1]    
    sc.pp.filter_genes(adata_hex, min_cells=min_cells)
    sc.pp.filter_cells(adata_hex, min_genes=min_genes)
    logger.info(
        "Filtered hexes from %s to %s and genes from %s to %s",
        orig_hexes, adata_hex.shape[0], orig_genes, adata_hex.shape[1]
    )

    adata_hex.layers["counts"] = adata_hex.X.copy()
    sc.pp.normalize_total(adata_hex, target_sum=1e4)
    sc.pp.log1p(adata_hex)
    sc.pp.regress_out(adata_hex, ["tz_count"])
    sc.pp.scale(adata_hex, max_value=10)

    sc.pp.highly_variable_genes(adata_hex, n_top_genes=n_top_genes)
    sc.tl.pca(adata_hex, n_comps=pca_comps)
    sc.pp.neighbors(adata_hex)
    sc.tl.umap(adata_hex)
    sc.tl.leiden(adata_hex, flavor="igraph", n_iterations=5, resolution=leiden_resolution)

    sc.tl.rank_genes_groups(adata_hex, groupby="leiden", method="wilcoxon")
    return adata_hex

def run_transcript_qc(
    sdata,
    points_key: str,
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int | None = None,
    min_density: float | None = None,
    table_key: str | None = None,
    qc_shapes_key: str = "transcript_qc_shapes",
    persist_qc_shapes: bool = True,
    force_recompute: bool = False,
):
    adata_hex, table_key = get_or_create_hex_adata(
        sdata,
        points_key=points_key,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        table_key=table_key,
        force_recompute=force_recompute,
    )
    adata_qc, grid = transcript_qc(
        adata_hex,
        min_transcripts=min_transcripts,
        min_density=min_density,
    )

    if persist_qc_shapes and hasattr(sdata, "__setitem__"):
        try:
            # Override existing on-disk element when available.
            if hasattr(sdata, "elements_paths_on_disk") and hasattr(
                sdata, "delete_element_from_disk"
            ):
                shapes_path = f"shapes/{qc_shapes_key}"
                if shapes_path in sdata.elements_paths_on_disk():
                    sdata.delete_element_from_disk(qc_shapes_key)

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                from spatialdata.models import ShapesModel
            grid_shapes = ShapesModel.parse(grid)
            sdata[qc_shapes_key] = grid_shapes
            if hasattr(sdata, "write_element"):
                sdata.write_element(qc_shapes_key)
        except Exception as e:
            logger.warning(
                "Failed to persist transcript QC shapes '%s': %s",
                qc_shapes_key,
                e,
                exc_info=True,
            )

    return adata_qc, grid, table_key

def run_cluster_hexes(
    sdata,
    points_key: str,
    hex_size: float = 30,
    hex_overlap: float = 0,
    gene_col: str = "gene",
    x_col: str = "x",
    y_col: str = "y",
    min_transcripts: int | None = None,
    min_density: float | None = None,
    leiden_resolution: float = 1.0,
    min_cells: int = 10,
    min_genes: int = 100,
    n_top_genes: int = 200,
    pca_comps: int = 50,
    table_key: str | None = None,
    cluster_table_key: str | None = None,
    force_recompute: bool = False,
):
    adata_hex, table_key = get_or_create_hex_adata(
        sdata,
        points_key=points_key,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gene_col=gene_col,
        x_col=x_col,
        y_col=y_col,
        table_key=table_key,
        force_recompute=force_recompute,
    )

    adata_filt, _ = transcript_qc(
        adata_hex,
        min_transcripts=min_transcripts,
        min_density=min_density,
    )

    adata_clustered = cluster_hexes(
        adata_filt,
        leiden_resolution=leiden_resolution,
        min_cells=min_cells,
        min_genes=min_genes,
        n_top_genes=n_top_genes,
        pca_comps=pca_comps,
    )

    cluster_key = cluster_table_key or f"{_hex_table_key(hex_size, hex_overlap)}_clustered"
    
    if cluster_key in sdata.tables:
        warning_msg = f"Overwriting existing table with key '{cluster_key}' in sdata.tables."
        warnings.warn(warning_msg)
        logger.warning(warning_msg)
        sdata.delete_element_from_disk(cluster_key)  # Remove the old table from disk
    sdata.tables[cluster_key] = adata_clustered
    sdata.write_element(cluster_key)

    return adata_clustered, cluster_key
