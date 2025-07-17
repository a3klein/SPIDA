# import os
# import sys
# import fire # type: ignore
# import warnings

# from dotenv import load_dotenv # type: ignore
# load_dotenv()

# # TODO: What if I defined a logger at the top of each cli module file? then all downstream functions use the same inherited logger? 

# class io_cli(): 
#     """
#     A Class for handling IO operations in SPIDA.
#     This class provides methods to ingest new spatial data, load segmentations, and plot results.

#     Methods: 
#     - ingest_region: Ingest a specific region of an experiment into a spatialdata object.
#     - ingest_all: Ingest all regions of an experiment into spatialdata objects.
#     - load_segmentation_region: Load segmentation data for a specific region into a spatialdata object.
#     - load_segmentation_all: Load segmentation data for all regions of an experiment into spatialdata objects.
#     """

#     def ingest_region(self, 
#                       exp_name:str, 
#                       reg_name:str,
#                       type:str="merscope",  
#                       prefix_name:str="default", 
#                       source:str="machine", 
#                       plot:bool=False):
#         """
#         Entry point for ingesting spatial data / segmentation outputs into spatialdata objects.

#         Parameters:
#         exp_name (str): Name of the experiment.
#         reg_name (str): Name of the region.
#         type (str): Type of the data to ingest (default is "merscope").
#         prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").

#         """
#         with warnings.catch_warnings():
#             warnings.filterwarnings("ignore")
#             from spida.S.io.main import ingest_region as func
#         func(exp_name, reg_name, type=type, prefix_name=prefix_name, source=source, plot=plot)

#     def ingest_all(self, 
#                    exp_name:str, 
#                    type:str="merscope", 
#                    prefix_name:str="default",
#                    source:str="machine",
#                    plot:bool=False):
#         """
#         Ingest all regions of an experiment into spatialdata objects.

#         Parameters:
#         exp_name (str): Name of the experiment.
#         type (str): Type of the data to ingest (default is "merscope").
#         prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
#         source (str): Source of the data (default is "machine").
#         """

#         with warnings.catch_warnings():
#             warnings.filterwarnings("ignore")
#             from spida.S.io.main import ingest_all as func
#         func(exp_name, type=type, prefix_name=prefix_name, source=source, plot=plot)


#     def load_segmentation_region(self, 
#                                  exp_name:str, 
#                                  reg_name:str, 
#                                  seg_dir:str,
#                                  type:str="vpt", 
#                                  prefix_name:str="default", 
#                                  plot:bool=False,
#                                  **load_kwargs
#                                 ):
#         """
#         Load segmentation data into spatialdata objects.

#         Parameters:
#         exp_name (str): Name of the experiment.
#         reg_name (str): Name of the region.
#         type (str): Type of the segmentation data to load (default is "vpt").
#         prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
#         """

#         with warnings.catch_warnings():
#             warnings.filterwarnings("ignore")
#             from spida.S.io.main import load_segmentation_region as func
#         func(exp_name, reg_name, seg_dir=seg_dir, type=type, prefix_name=prefix_name, plot=plot, **load_kwargs)

    
#     def load_segmentation_all(self, 
#                               exp_name:str,
#                             seg_dir:str,
#                             type:str="vpt", 
#                             prefix_name:str="default",
#                             plot:bool=False,
#                             ):
#         """
#         Load segmentation data for all regions of an experiment into spatialdata objects.

#         Parameters:
#         exp_name (str): Name of the experiment.
#         seg_dir (str): Directory containing the segmentation data.
#         type (str): Type of the segmentation data to load (default is "vpt").
#         prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
#         """

#         with warnings.catch_warnings():
#             warnings.filterwarnings("ignore")
#             from spida.S.io.main import load_segmentation_all as func
#         func(exp_name, seg_dir=seg_dir, type=type, prefix_name=prefix_name, plot=plot)


# if __name__ == "__main__": 
#     fire.Fire(io_cli)
