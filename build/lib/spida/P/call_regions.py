# imports
import os
from pathlib import Path
import logging
import warnings

import numpy as np
import pandas as pd
import spatialdata as sd

import geopandas as gpd
from sklearn.mixture import GaussianMixture
import libpysal as lps
import networkx as nx

from spida.utilities.tiling import create_hexagonal_grid

logger = logging.getLogger(__package__)
warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

### TODO: Figure out the plotting mechanism
def call_regions(
    zarr_path : str | Path, 
    geoms_name : str,
    points_key : str,
    transfer_genes : str | list[str] = "BCAS1",
    save_geoms_path : str | None = None,
    dsc_comp_min_size = 5,
    hex_size = 50,
    hex_overlap = 0,
    gmm_ncomp : int | str = "auto",
    gmm_cov_type : str = "full",
    gene_agreement_thr = 0.25,
    top_n_comp: int = 1,
    gen_plots: bool = False,
    plot_save_path : str | Path = None,
    image_key : str = None,
):
    """
    Call regions based on specific gene expression in space using hexagonal tiling and GMM clustering
    zarr_path : str | Path
        Path to the zarr file containing the spatial transcriptomics data.
    geoms_name : str
        Name for the output geometries in the SpatialData object.
    points_key : str
        Key for the points model in the SpatialData object.
    transfer_genes : str | list[str], optional
        Gene or list of genes to use for region calling. Default is "BCAS1" for White Matter.
    save_geoms_path : str | None, optional
        Path to save the resulting geometries. If None, geometries are not saved. Default is None.
    dsc_comp_min_size : int, optional
        Minimum number of hexes in a disconnected component to consider as a region. Default is 5.
    hex_size : int, optional
        Size of the hexagons in the grid. Default is 50.
    hex_overlap : int, optional
        Overlap between hexagons in the grid. Default is 0.
    gmm_ncomp : int | str, optional
        Number of components for the Gaussian Mixture Model. If "auto", selects the best number based on BIC. Default is "auto".
    gmm_cov_type : str, optional
        Covariance type for the Gaussian Mixture Model. Default is "full".
    gene_agreement_thr : float, optional
        Portion of genes that must agree to call a hexagon as part of the region. Default is 0.25.
    top_n_comp: int, optional
        Number of top components to consider for region calling. Default is 1.
    gen_plots : bool, optional
        Whether to generate plots for the region calling results. Default is False.
    plot_save_path : str | Path, optional
        Path to save the generated plots. If None, plots are not saved. Default is None
    image_key : str, optional
        Key for the image model in the SpatialData object for plotting. Default is None.
    """

    if isinstance(transfer_genes, str):
        transfer_genes = [transfer_genes]

    run_goodness_of_fit = False
    if gmm_ncomp == "auto": 
        run_goodness_of_fit = True
    elif isinstance(gmm_ncomp, str) and (gmm_ncomp.isdigit()):
        gmm_ncomp = int(gmm_ncomp)
    else: 
        raise ValueError("gmm_ncomp must be an integer or 'auto'.")


    try: 
        sdata = sd.read_zarr(zarr_path)
    except Exception as e:
        raise RuntimeError(f"Failed to read zarr from path {zarr_path}: {e}")

    ret = call_region_worker(
        sdata=sdata,
        geoms_name=geoms_name,
        points_key=points_key,
        transfer_genes=transfer_genes,
        save_geoms_path=save_geoms_path,
        dsc_comp_min_size=dsc_comp_min_size,
        hex_size=hex_size,
        hex_overlap=hex_overlap,
        gmm_ncomp=gmm_ncomp,
        gmm_cov_type=gmm_cov_type,
        gene_agreement_thr=gene_agreement_thr,
        run_goodness_of_fit=run_goodness_of_fit,
        top_n_comp=top_n_comp,
        plot=gen_plots,
    )

    #TODO: import .pl workers for plotting if gen_plots is True
    if gen_plots and plot_save_path is not None:
        import matplotlib.pyplot as plt
        import spatialdata_plot as sdp # type: ignore
        plt.rcParams['figure.dpi'] = 150

        df_gene_dict, grid = ret

        fig = plt.figure(constrained_layout=True, figsize=(15, 5 * len(transfer_genes)))
        fig.suptitle("Region Calling by Gene", y=1.05, fontsize=16)
        subfigs = fig.subfigures(len(transfer_genes), 1, hspace=0.1)
        if len(transfer_genes) == 1:
            subfigs = [subfigs]
        for _gene, subfig in zip(transfer_genes, subfigs):
            df_gene = df_gene_dict[_gene]
            subfig.suptitle(_gene, y=1.05, fontsize=14)
            axes = subfig.subplots(1, 3)

            ax = axes[0]
            grid.plot(ax=ax, column=f"{_gene}_count", cmap="YlOrRd", edgecolor='k', alpha=0.7, linewidth=0.5).axis("off")
            ax.set_title("counts", fontsize=12)

            ax = axes[1]
            minv, maxv = df_gene[_gene + '_count'].min(), df_gene[_gene + '_count'].max()
            ax.hist(df_gene[_gene + '_count'], bins=30, range=(minv, maxv), color='lightgray', edgecolor='k', alpha=0.7, histtype='stepfilled')
            ax.hist(df_gene.loc[df_gene['predict'] == 1, _gene + '_count'], bins=30, range=(minv, maxv), color='green', edgecolor='k', alpha=0.7, histtype='stepfilled')
            # sns.histplot(ax=ax, data=df_gene, x=_gene + '_count', bins=30, hue="predict", palette="viridis", edgecolor='k')
            # sns.despine(ax=ax)
            ax.set_title(f'GMM Prediction', fontsize=12)
            
            ax = axes[2]
            grid.plot(ax=ax, column=f"{_gene}_predict", cmap="YlOrRd", edgecolor='k', alpha=0.7, linewidth=0.5).axis("off")
            ax.set_title("prediction", fontsize=12)
        plt.savefig(f"{plot_save_path}/region_calling_by_gene.png")
        plt.close()

        fig, axes = plt.subplots(ncols = 2, nrows = 1, figsize=(12,5))
        cs = "global"
        ch = "MBP"        
        (
            sdata.pl.render_images(image_key, channel=ch, cmap="gray")
            .pl.render_shapes(geoms_name, color="none", outline_color="red", outline_width=2, outline_alpha=1, fill_alpha=0)
            .pl.show(ax=axes[0], coordinate_systems=cs)
        )
        axes[0].set_title(f"WM Regions")

        ch = "DAPI"
        (
            sdata.pl.render_images(image_key, channel=ch, cmap="gray")
            .pl.render_shapes(geoms_name, color="none", outline_color="red", outline_width=2, outline_alpha=1, fill_alpha=0)
            .pl.show(ax=axes[1], coordinate_systems=cs)
        )
        axes[1].set_title(f"WM Regions")
        plt.savefig(f"{plot_save_path}/region_calling_on_image.png")
        plt.close()
        


# TODO: implement a second module based on spatial weights of the points with no hexes
def call_region_worker(
    sdata : sd.SpatialData,
    geoms_name : str,
    points_key : str,
    transfer_genes : str | list[str] = "BCAS1",
    save_geoms_path : str | None = None,
    dsc_comp_min_size = 5,
    hex_size = 50,
    hex_overlap = 0,
    gmm_ncomp : int | str = 3,
    gmm_cov_type : str = "full",
    gene_agreement_thr = 0.75,
    run_goodness_of_fit : bool = False,
    top_n_comp: int = 1,
    plot : bool = False,
): 
    """
    Call regions based on specific gene expression in space using hexagonal tiling and GMM clustering.

    sdata : sd.SpatialData
        SpatialData object containing the spatial transcriptomics data.
    geoms_name : str
        Name for the output geometries in the SpatialData object.
    points_key : str
        Key for the points model in the SpatialData object.
    transfer_genes : str | list[str], optional
        Gene or list of genes to use for region calling. Default is "BCAS1" for White Matter.
    save_geoms_path : str | None, optional
        Path to save the resulting geometries. If None, geometries are not saved. Default is None.
    hex_size : int, optional
        Size of the hexagons in the grid. Default is 50.
    hex_overlap : int, optional
        Overlap between hexagons in the grid. Default is 0.
    gmm_ncomp : int, optional
        Number of components for the Gaussian Mixture Model. Default is 3.
    top_n_comp: int, optional
        Number of top components to consider for region calling. Default is 1.
    gmm_cov_type : str, optional
        Covariance type for the Gaussian Mixture Model. Default is "full".
    gene_agreement_thr : float, optional
        Portion of genes that must agree to call a hexagon as part of the region. Default is 0.75.
    dsc_comp_min_size : int, optional
        Minimum number of hexes in a disconnected component to consider as a region. Default is 5.
    """
    # TODO: Figure out how the plotting mechanism will work here. 

    # Get the points from the spatialdata object and make sure genes are categorical
    fts = sdata[points_key].compute()
    fts = fts.reset_index()
    fts['gene'] = fts['gene'].astype("category")

    # Load the points object into a geopandas dataframe and subset only to the genes of interest
    gdf = gpd.GeoDataFrame(fts, geometry=gpd.points_from_xy(fts['x'], fts['y'])) 
    gdf = gdf.loc[gdf['gene'].isin(transfer_genes)].copy()
    print(f"Loaded points for region calling. Total points: {len(gdf)}")

    # Get the Distance Band Weights     
    kd = lps.cg.KDTree(gdf.geometry.apply(lambda geom: (geom.x, geom.y)).tolist())
    wnndb = lps.weights.DistanceBand(kd, threshold=100, p=1, binary=False, alpha=-1, ids=gdf.index.tolist())

    # For tz counts in each hex (module 1)
    for _gene in transfer_genes: 
        gdf[_gene] = (gdf['gene'] == _gene).astype(int)

    # Create the grid
    total_bounds = gdf.total_bounds  # (minx, miny, maxx, maxy)
    grid = create_hexagonal_grid(total_bounds, hex_size, overlap=hex_overlap)

    # Join the grid for the accumulated counts
    grid['hex_id'] = grid.index.astype(str)
    grid = grid.set_index("hex_id")
    joint_grid = gpd.sjoin(grid, gdf, how="inner", predicate="contains")
    print(f"Created hexagonal grid with {len(grid)} hexes.")

    # Get the number of total tz's per hex (from the transfer list)
    grid['cell_count'] = joint_grid.groupby('hex_id').size()
    # Get the number of gene-specific tz's per hex
    for _gene in transfer_genes:
        grid[f'{_gene}_count'] = joint_grid.groupby('hex_id')[_gene].sum()

    if plot:
        df_gene_dict = {}
    for i, _gene in enumerate(transfer_genes):
        print(f"Fitting GMM for {_gene}...")
        df_gene = grid[[_gene + '_count']].dropna().copy() # Get the gene count / hex names 
        # Choosing an ncomp best on BIC if set to auto
        if run_goodness_of_fit:
            lowest_bic = np.inf
            best_ncomp = 1
            n_components_range = range(1, 4)
            for n_components in n_components_range:
                gmm = GaussianMixture(n_components=n_components, random_state=0, covariance_type=gmm_cov_type)
                gmm.fit(df_gene[[_gene + '_count']].values)
                bic = gmm.bic(df_gene[[_gene + '_count']].values)
                logger.info(f"GMM BIC for {_gene} with {n_components} components: {bic}")
                if bic < lowest_bic:
                    lowest_bic = bic
                    best_ncomp = n_components
            gmm_ncomp = best_ncomp
            print(f"Selected {gmm_ncomp} components for {_gene} based on BIC.")
        # fitting gmm
        gmm = GaussianMixture(n_components=gmm_ncomp, random_state=0, covariance_type=gmm_cov_type).fit(df_gene[[_gene + '_count']].values)
        gene_prediction = gmm.predict(grid[[_gene + '_count']].dropna().values)
        df_gene['predict'] = gene_prediction
        pred_vals = df_gene.groupby('predict').mean()[_gene + "_count"].nlargest(top_n_comp).index.tolist()
        # pred_val = df_gene.groupby("predict").mean().idxmax(axis=0)[0]
        df_gene['predict'] = (df_gene['predict'].isin(pred_vals)).astype(int)
        grid.loc[df_gene.index, _gene + '_predict'] = df_gene['predict']
        if plot:
            df_gene_dict[_gene] = df_gene
    
    # Accumulating the hex specific predictions per gene 
    hex_ids = {}
    for _gene in transfer_genes: 
        hexes = set(grid[grid[_gene + "_predict"] == 1].index)
        for _h in hexes: 
            hex_ids[_h] = 1 if _h not in hex_ids else hex_ids[_h] + 1

    # Getting the final hexes based on agreement threshold
    df_hids = pd.DataFrame.from_dict(hex_ids, orient="index")
    df_hids = df_hids / len(transfer_genes)
    chosen_cells = df_hids[df_hids[0] >= gene_agreement_thr].index
    chosen_cells = grid.loc[list(chosen_cells)]
    print(f"Chosen Hexes: {len(chosen_cells)}")

    # Turning hexes into a CCP to aggregate them 
    chosen_cells = chosen_cells.reset_index()
    W = lps.weights.Queen.from_dataframe(chosen_cells)
    G = W.to_networkx()
    connected_components = list(nx.connected_components(G))
    disconnected_comp = [comp for comp in connected_components]

    # aggregating all connected components into geometries 
    chosen_cells['comp'] = -1
    geoms = []
    for i, disc in enumerate(disconnected_comp): 
        if len(disc) < dsc_comp_min_size:
            continue
        chosen_cells.loc[list(disc), "comp"] = i
        temp = chosen_cells[chosen_cells.index.isin(disc)]
        geom = temp.union_all().convex_hull
        geoms.append(geom)

    # Combine all geometries into a single GeoDataFrame
    gdf_geoms = gpd.GeoDataFrame(geometry=geoms)
    union = gdf_geoms.union_all()
    gdf_geoms = gpd.GeoDataFrame(geometry=[union])

    # Save Geoms if save_geoms is set
    gdf_geoms = gdf_geoms.explode()
    if save_geoms_path is not None:
        gdf_geoms.to_file(save_geoms_path, driver="GPKG")


    print(f"Final geometries created: {len(gdf_geoms)}. Saving to spatialdata")

    if geoms_name in sdata.shapes.keys():
        logger.warning(f"Overwriting existing geometries at {geoms_name} in SpatialData.")
        sdata.delete_element_from_disk(geoms_name)  # Remove the old table from disk

    sdata[geoms_name] = sd.models.ShapesModel().parse(gdf_geoms)
    sd.transformations.set_transformation(
        sdata[geoms_name],
        sd.transformations.get_transformation(sdata[points_key], to_coordinate_system="pixel"),
        to_coordinate_system="pixel"
    )

    sd.transformations.set_transformation(
        sdata[geoms_name],
        sd.transformations.get_transformation(sdata[points_key], to_coordinate_system="global"),
        to_coordinate_system="global"
    )
    sdata.write_element(geoms_name, overwrite=True)

    if plot: 
        return df_gene_dict, grid
    return None