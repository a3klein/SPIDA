import json

import polars as pl
import numpy as np
import anndata as ad

from spida._utilities import _validate_adata
from spida._constants import (CELL_ID, CELL_X, CELL_Y, CELL_FOV, CELL_VOLUME,
                            PROSEG_PRESET, DEFAULT_PRESET) # type: ignore

class Filter():

    def __init__(self, adata:ad.AnnData,
                 exp_name:str, 
                 reg_name:str,
                 seg_name:str,
                 donor_name:str=None,
                 PRESET:dict=DEFAULT_PRESET):  
        self.adata = adata
        self.exp_name = exp_name
        self.reg_name = reg_name
        self.seg_name = seg_name
        self.donor_name = donor_name
        self.PRESET = PRESET

    def load_info(self, df:pl.DataFrame) -> pl.DataFrame:
        """
        Load information into a Polars DataFrame and add metadata columns.
        """

        df = df.with_columns(
                experiment = pl.lit(self.exp_name),
                region = pl.lit(self.reg_name),
                segmentation = pl.lit(self.seg_name),
                donor = pl.lit(self.donor_name),
            )
        
        df = df.with_columns(pl.concat_str([pl.col(CELL_ID),
                                            pl.col("experiment"),
                                            pl.col("region"),
                                            pl.col("segmentation"),
                                            ], separator=".").alias("Index"))
        return df

    def load_metadata(self) -> pl.DataFrame: 
        """
        Load metadata from an AnnData object and return it as a Polars DataFrame.
        """
        try: 
            df_obs = pl.DataFrame(self.adata.obs).rename(self.PRESET['meta_map'], strict=True)
        except pl.exceptions.ColumnNotFoundError as e: 
            df_obs = pl.DataFrame(self.adata.obs)
        
        df_meta = df_obs.select([CELL_ID, CELL_X, CELL_Y, CELL_FOV, CELL_VOLUME])
        df_meta = self.load_info(df_meta)
        return df_meta 


    def load_cbg(self) -> pl.DataFrame: 
        """
        Load cell-by-gene data from an AnnData object and return it as a Polars DataFrame.
        """        
        df_cbg = pl.DataFrame(self.adata.X.copy()).sum_horizontal().to_frame()
        cell_ids = self.adata.obs.index.str.split('.', n=1, expand=True).get_level_values(0).to_frame(name="index")
        df_cbg = df_cbg.with_columns(pl.DataFrame(cell_ids)) # Need to use CELL_ID not index as the new index
        # 2 Indexes renaming because capitalization was not consistent in the past
        df_cbg = df_cbg.rename({"index" : CELL_ID, "Index" : CELL_ID, "sum" : "nCount_RNA"}, strict=False )
        df_cbg = df_cbg.with_columns(nFeature_RNA=(pl.DataFrame(self.adata.X.copy()) != 0).sum_horizontal())
        df_cbg = df_cbg.with_columns(nBlank=pl.DataFrame(self.adata.copy().obsm['blank']).sum_horizontal())
        df_cbg = self.load_info(df_cbg)
        return df_cbg
    
    def get_features(self) -> pl.DataFrame: 
        """
        Get features by joining cell-by-gene data with metadata.
        """
        df_join = self.load_cbg().join(self.load_metadata(), on="Index", how="inner", coalesce=True)
        # Adding the count / volume ratio per cell 
        df_join = df_join.with_columns(pl.col("nCount_RNA").truediv(pl.col("volume")).alias("nCount_RNA_per_Volume"))

        return df_join.select(['Index', CELL_ID, "experiment", "region", "segmentation", "donor", 
                               CELL_X, CELL_Y, CELL_FOV, CELL_VOLUME,
                               'nCount_RNA', 'nFeature_RNA', 'nBlank', 'nCount_RNA_per_Volume'])
    

    def filter_cells(self, cutoffs:dict) -> pl.DataFrame: 
        """
        Filter cells based on specified cutoffs for various features.
        """

        df_feature = self.get_features()
        # Getting the distribution based feature cutoffs: 
        volume_min=cutoffs['volume_min']
        volume_max=cutoffs['volume_max']
        n_count_min=cutoffs['n_count_min']
        n_count_max=cutoffs['n_count_max']
        n_gene_min=cutoffs['n_gene_min']
        n_gene_max=cutoffs['n_gene_max']
        n_blank_min=cutoffs['n_blank_min']
        n_blank_max=cutoffs['n_blank_max']
        ncpv_min_quantile=cutoffs['ncpv_min_quantile']
        ncpv_max_quantile=cutoffs['ncpv_max_quantile']

        vol_max_mult = cutoffs.get('vol_max_mult', 3)
        max_vol_by_median = cutoffs.get('max_vol_by_median', False)
        count_max_mult = cutoffs.get('count_max_mult', 3)
        max_count_by_median = cutoffs.get('max_count_by_median', False)

        if max_vol_by_median: 
            volume_max = max(df_feature['volume'].median() * vol_max_mult, volume_max)
        if max_count_by_median: 
            n_count_max = max(df_feature['nCount_RNA'].median() * count_max_mult, n_count_max)
        n_count_per_volume = df_feature['nCount_RNA_per_Volume'].to_numpy()
        qMax = np.quantile(n_count_per_volume, ncpv_max_quantile)
        qMin = np.quantile(n_count_per_volume, ncpv_min_quantile)
        n_count_per_volume_min = qMin
        n_count_per_volume_max = qMax

        # Store cutoffs in the AnnData object for reference
        self.adata.uns['cutoffs']['volume_max'] = volume_max  
        self.adata.uns['cutoffs']['n_count_per_volume_min'] = n_count_per_volume_min
        self.adata.uns['cutoffs']['n_count_per_volume_max'] = n_count_per_volume_max
        
        # Apply the filtering 
        judge = df_feature.filter(
            pl.col("volume").is_between(volume_min, volume_max) & 
            pl.col("nCount_RNA").is_between(n_count_min, n_count_max) &
            pl.col("nFeature_RNA").is_between(n_gene_min, n_gene_max) &
            pl.col("nBlank").is_between(n_blank_min, n_blank_max) &
            pl.col("nCount_RNA_per_Volume").is_between(n_count_per_volume_min, n_count_per_volume_max)
            ).select(CELL_ID)
        
        
        # df_tt.select(pl.col(CELL_ID).is_in(judge[CELL_ID].implode()))
        df_feature = df_feature.with_columns(pl.col(CELL_ID).is_in(judge[CELL_ID].implode()).alias("pass_qc"))

        return df_feature
    


def run_filtering(adata:ad.AnnData, 
                  exp_name:str,
                  reg_name:str,
                  prefix_name:str,
                  donor_name=None,
                  cutoffs_path:str="/ceph/cephatlas/aklein/bican/reference/filtering_cutoffs.json"): 
    """
    Run the filtering process on the spatialdata object.
    This function initializes the Filter class and returns the filtered features.

    Parameters:
    sdata (sd.SpatialData): The spatialdata object containing the data to be filtered.
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    prefix_name (str): Name of the segmentation.
    donor_name (str, optional): Name of the donor. Defaults to None. 
    """

    ## Defining the filtering params
    seg_fam_map = {"default" : "default",
                   "proseg" : "proseg", 
                   "cellpose_nuclei" : "default",
                   "cellposeSAM" : "default",
                   "proseg_nuclei" : "proseg",
                   "proseg_SAM" : "proseg",
                   "mesmer" : "default"}
    seg_fam = seg_fam_map.get(prefix_name, "default")
    PRESET = PROSEG_PRESET if seg_fam == "proseg" else DEFAULT_PRESET # PRESET for column names
    
    # Cutoffs for the filtering 
    with open(cutoffs_path, 'r') as f: 
        cutoffs = json.load(f)
    adata.uns['cutoffs'] = cutoffs  # Store cutoffs in the AnnData object for reference

    # Initialize the filter class
    filter_instance = Filter(adata, exp_name, reg_name, prefix_name, donor_name=donor_name, PRESET=PRESET)

    # Perform the filtering 
    filtered_features = filter_instance.filter_cells(cutoffs)

    adata.obs = filtered_features.to_pandas().set_index("Index").join(adata.obs.drop(columns=['region', 'fov', 'volume']), on=CELL_ID, how='inner')
    adata = _validate_adata(adata)  # Validate the AnnData object
    
    return adata.copy()
