
from importlib.metadata import version

from . import io, P, I, segmentation, __main__ 
__all__ = ["io", "P", "I", "segmentation", "__main__"]
# __version__ = version("spida")


# import os 
# import sys
# from pathlib import Path
# import fire

# spida_path = "/ceph/cephatlas/aklein/spida/src"
# sys.path.append(spida_path)
# from spida._templates import *


# def gen_run(type:str,
#             exp_name:str, 
#             prefix_name:str, 
#             seg_dir:Path="",
#             seg_type:str="",
#             cutoffs_path:Path=None, 
#             plot_load:bool=False,
#             plot_filter:bool=False,
#             plot_setup:bool=False,
#             ):
#     """
#     Generate and run the main function for filtering cells in a region.
#     """

#     if type == "start": 
#         def_template.format(
#             EXPERIMENT=exp_name,
#             PREFIX=prefix_name,
#             LOAD_PLOT=plot_load,
#             PLOT_FILTERING=plot_filter,
#             PLOT_ADATA=plot_setup,
#             CUTOFF_JSON_PATH=cutoffs_path if cutoffs_path else "",
#         )
#     elif type == "load_seg": 
#         load_seg_template.format(
#             EXPERIMENT=exp_name,
#             SEG_DIR=seg_dir,
#             SEG_LOAD_PLOT=plot_load,
#             SEG_TYPE=seg_type,
#             PREFIX=prefix_name,
#             PLOT_FILTERING=plot_filter,
#             CUTOFF_JSON_PATH=cutoffs_path if cutoffs_path else "",
#             PLOT_ADATA=plot_setup,
#         )
    
# if __name__ == "__main__": 
#     fire.Fire({"run_template" : gen_run})

