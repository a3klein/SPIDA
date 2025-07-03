### Integration / Annotations main file 
import warnings
import fire # type: ignore
from dotenv import load_dotenv  # type: ignore
load_dotenv()


class I_cli(): 
    """
    Command line interface for integration and annotation tasks in the BICAN project.
    This class provides methods to backup AnnData objects, run ALLCools integration, 
    setup MapMyCells integration, and run MapMyCells annotation for specific regions or entire experiments.
    
    Methods:
    - backup_adata_region: Backup function for AnnData objects for a specific region.
    - backup_adata_experiment: Backup function for AnnData objects for an entire experiment.
    - allcools_integration_region: Run ALLCools integration on a given experiment and region.
    - allcools_integration_experiment: Run ALLCools integration for an entire experiment.
    - mmc_setup: Setup function for MapMyCells integration.
    - mmc_annotation_region: Run MapMyCells annotation on a given experiment and region.
    - mmc_annotation_experiment: Run MapMyCells annotation for an entire experiment.
    - moscot_integration: Placeholder for MOSCOT integration (not implemented yet).
    """
    
    def backup_adata_region(exp_name:str, reg_name:str, prefix_name:str, adata_path:str=None):
        """
        Backup function for AnnData objects. Due to env incompatibilities with 
        SpatialData and some of the annotation / integration libraries. 

        Parameters:
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        """

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I._utils import backup_adata_region as func
        func(exp_name, reg_name, prefix_name, adata_path)

    def backup_adata_experiment(exp_name:str, prefix_name:str, adata_path:str=None):
        """
        Backup function for AnnData objects for an entire experiment.
        
        Parameters:
        exp_name (str): Name of the experiment.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        adata_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        """
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I._utils import backup_adata_experiment as func
        func(exp_name, prefix_name, adata_path)

    
    def allcools_integration_region(exp_name:str,
                                    reg_name:str, 
                                    prefix_name:str, 
                                    ref_path:str,
                                    anndata_store_path:str=None, 
                                    annotations_store_path:str=None,
                                    **kwargs):
        """
        Run ALLCools integration on a given experiment and region.
        ENVIRONMENT = "preprocessing"

        Parameters:
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        ref_path (str): Path to the reference RNA AnnData object .
        anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
        to None.
        **kwargs: Additional keyword arguments for ALLCools integration.
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I.allcools import allcools_integration_region as func
        func(exp_name, reg_name, prefix_name, ref_path, anndata_store_path, annotations_store_path, **kwargs)

    def allcools_integration_experiment(exp_name:str,
                                        prefix_name:str, 
                                        ref_path:str,
                                        anndata_store_path:str=None, 
                                        annotations_store_path:str=None,
                                        **kwargs):
        """
        Run ALLCools integration for an entire experiment.

        Parameters:
        exp_name (str): Name of the experiment.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        ref_path (str): Path to the reference RNA AnnData object .
        anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults
        to None.
        **kwargs: Additional keyword arguments for ALLCools integration.
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I.allcools import allcools_integration_experiment as func
        func(exp_name, prefix_name, ref_path, anndata_store_path, annotations_store_path, **kwargs)


    def mmc_setup(ref_path:str,
                  heirarchy_list:list,
                  BRAIN_REGION:str,
                  CODEBOOK:str,
                  codebook_path:str=None,
                  mmc_store_path:str=None,
                  ref_norm:str="log2CPM", 
                  **kwargs
                ): 
        """
        Setup function for MapMyCells integration.
        
        Parameters:
        ref_path (str): Path to the reference data.
        heirarchy_list (list): List of hierarchy levels for the annotation.
        BRAIN_REGION (str): Brain region for the annotation.
        CODEBOOK (str): Codebook for the annotation.
        codebook_path (str, optional): Path to the codebook file. Defaults to None.
        mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
        ref_norm (str, optional): Normalization method for the reference AnnData.X. Defaults to "log2CPM".
        **kwargs: Additional keyword arguments.

        
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I.mmc import mmc_setup as func
        func(ref_path, BRAIN_REGION, CODEBOOK, codebook_path, heirarchy_list, mmc_store_path, ref_norm, **kwargs)

    def mmc_annotation_region(
                            exp_name:str, 
                            reg_name:str, 
                            prefix_name:str, 
                            BRAIN_REGION:str, 
                            CODEBOOK:str, 
                            mmc_store_path:str=None, 
                            anndata_store_path:str=None, 
                            annotations_store_path:str=None,
                            **kwargs
                            ):
        """
        Run MapMyCells annotation on a given experiment and region.
        
        Parameters: 
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        BRAIN_REGION (str): Brain region for the annotation.
        CODEBOOK (str): Codebook for the annotation.
        mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
        anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
        **kwargs: Additional keyword arguments.
        """
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I.mmc import mmc_annotation_region as func
        func(exp_name, reg_name, prefix_name, BRAIN_REGION, CODEBOOK, mmc_store_path,
              anndata_store_path, annotations_store_path, **kwargs)
        
    def mmc_annotation_experiment(exp_name:str, 
                                  prefix_name:str, 
                                  BRAIN_REGION:str, 
                                  CODEBOOK:str, 
                                  mmc_store_path:str=None, 
                                  anndata_store_path:str=None, 
                                  annotations_store_path:str=None,
                                  **kwargs
                                  ):
        """
        Run MapMyCells annotation for an entire experiment.
        
        Parameters:
        exp_name (str): Name of the experiment.
        prefix_name (str): Prefix for the keys in the spatialdata object.
        BRAIN_REGION (str): Brain region for the annotation.
        CODEBOOK (str): Codebook for the annotation.
        mmc_store_path (str, optional): Path to the store of MapMyCells markers. Defaults to None.
        anndata_store_path (str, optional): Path to the store of AnnData objects. Defaults to None.
        annotations_store_path (str, optional): Path to the store of annotation specific files. Defaults to None.
        **kwargs: Additional keyword arguments.
        """

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from spida.I.mmc import mmc_annotation_experiment as func
        func(exp_name, prefix_name, BRAIN_REGION, CODEBOOK, mmc_store_path,
              anndata_store_path, annotations_store_path, **kwargs)
    
    def moscot_integration(): 
        raise NotImplementedError("MOSCOT integration is not implemented yet.")



if __name__ == "__main__":
    fire.Fire(I_cli)
    # fire.Fire({"backup_adata_region" : backup_adata_region,
    #            "backup_adata_experiment": backup_adata_experiment,
    #             "allcools_integration_region" : allcools_integration_region, 
    #             "allcools_integration_experiment": allcools_experiment,
    #            "setup_mmc": mmc_setup,
    #            "mmc_annotation_region": mmc_annotation_region,
    #            "mmc_annotation_experiment": mmc_annotation_experiment,
    #            "moscot_integration": moscot_integration})