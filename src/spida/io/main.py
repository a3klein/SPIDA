import sys
import glob
import pathlib
import fire # type: ignore
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import spatialdata as sd

from spida._utilities import _gen_keys
from spida._constants import *
from spida.io.ingest_exp import read_merscope, load_vpt_segmentation, load_proseg_segmentation
from spida.pl import plot_images, plot_shapes, plot_points, plot_overlap

from matplotlib.backends.backend_pdf import PdfPages

# Setting Logging for this module 
import logging
logging.basicConfig(filename="/ceph/cephatlas/aklein/spida/tests/spida.log", level=logging.INFO)


class io_cli(): 
    """
    A Class for handling IO operations in SPIDA.
    This class provides methods to ingest new spatial data, load segmentations, and plot results.

    Methods: 
    - ingest_region: Ingest a specific region of an experiment into a spatialdata object.
    - ingest_all: Ingest all regions of an experiment into spatialdata objects.
    - load_segmentation_region: Load segmentation data for a specific region into a spatialdata object.
    - load_segmentation_all: Load segmentation data for all regions of an experiment into spatialdata objects.
    """

    def ingest_region(self, 
                      exp_name:str, 
                      reg_name:str,
                      type:str="merscope",  
                      prefix_name:str="default", 
                      source:str="machine", 
                      plot:bool=False):
        """
        Entry point for ingesting spatial data / segmentation outputs into spatialdata objects.

        Parameters:
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        type (str): Type of the data to ingest (default is "merscope").
        prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").

        """

        print("INGESTING REGION; EXPERIMENT %s, REGION %s " %(exp_name, reg_name) )
        KEYS = _gen_keys(prefix_name, exp_name, reg_name)

        if type == "merscope" and source == "machine":
            root_path = "/ceph/cephatlas/merscope_data/processed"
            processed_path = "/data/aklein/bican_zarr"
            input_path = f"{root_path}/{exp_name}/out/{reg_name}"
            zarr_path = f"{processed_path}/{exp_name}/{reg_name}"

            sdata = read_merscope(input_path, zarr_path, exp_name=exp_name, reg_name=reg_name, prefix_name=prefix_name)

            image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
            image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())

            print(sdata.tables.keys())
            if plot: 
                image_path = f"/ceph/cephatlas/aklein/bican/images/{exp_name}/default/{reg_name}/pixi-ing.pdf"
                pathlib.Path(image_path).parent.mkdir(parents=True, exist_ok=True)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    with PdfPages(image_path) as pdf:
                        plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                        plot_shapes(sdata, KEYS[SHAPES_KEY], table_name=KEYS[TABLE_KEY], cs="pixel", pdf_file=pdf)
                        plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                        plot_overlap(sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[POINTS_KEY], KEYS[TABLE_KEY], image_scale_keys, cs="pixel", pdf_file=pdf)


    def ingest_all(self, 
                   exp_name:str, 
                   type:str="merscope", 
                   prefix_name:str="default",
                   source:str="machine",
                   plot:bool=False):
        """
        Ingest all regions of an experiment into spatialdata objects.

        Parameters:
        exp_name (str): Name of the experiment.
        type (str): Type of the data to ingest (default is "merscope").
        prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
        source (str): Source of the data (default is "machine").
        """

        if type == "merscope" and source == "machine":
            root_path = "/ceph/cephatlas/merscope_data/processed"
            input_path = f"{root_path}/{exp_name}/out"
            regions = glob.glob(f"{input_path}/region_*")
            for reg in regions: 
                self.ingest_region(exp_name, reg.split("/")[-1], type=type, prefix_name=prefix_name, source=source, plot=plot)



    def load_segmentation_region(self, 
                                 exp_name:str, 
                                 reg_name:str, 
                                 seg_dir:str,
                                 type:str="vpt", 
                                 prefix_name:str="default", 
                                 plot:bool=False
                                ):
        """
        Load segmentation data into spatialdata objects.

        Parameters:
        exp_name (str): Name of the experiment.
        reg_name (str): Name of the region.
        type (str): Type of the segmentation data to load (default is "vpt").
        prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
        """

        print("LOADING SEGMENTATION; EXPERIMENT %s, REGION %s, SEGMENTATION %s" %(exp_name, reg_name, type) )


        zarr_path = f"/data/aklein/bican_zarr/{exp_name}/{reg_name}"
        print(zarr_path)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            sdata = sd.read_zarr(zarr_path)

        if type == "vpt":
            sdata = load_vpt_segmentation(sdata, exp_name, reg_name, vpt_path=seg_dir, prefix_name=prefix_name)
        elif type == "proseg":
            sdata = load_proseg_segmentation(sdata, exp_name, reg_name, proseg_path=seg_dir, prefix_name=prefix_name)

        if plot: 
            # plot params
            KEYS = _gen_keys(prefix_name, exp_name, reg_name)
            image_channels = sd.models.get_channel_names(sdata[KEYS[IMAGE_KEY]])
            image_scale_keys = list(sdata[KEYS[IMAGE_KEY]].keys())
            
            # define plot pdf
            image_path = f"/ceph/cephatlas/aklein/bican/images/{exp_name}/{prefix_name}/{reg_name}/pixi-load.pdf"
            pathlib.Path(image_path).parent.mkdir(parents=True, exist_ok=True)
            
            # plot
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                with PdfPages(image_path) as pdf:
                    plot_images(sdata, KEYS[IMAGE_KEY], image_scale_keys, image_channels, cs="global", pdf_file=pdf)
                    plot_shapes(sdata, KEYS[SHAPES_KEY], table_name=KEYS[TABLE_KEY], cs="pixel", pdf_file=pdf)
                    plot_points(sdata, KEYS[POINTS_KEY], KEYS[TABLE_KEY], cs="pixel", cmap="tab10", pdf_file=pdf)
                    plot_overlap(sdata, KEYS[IMAGE_KEY], KEYS[SHAPES_KEY], KEYS[POINTS_KEY], KEYS[TABLE_KEY], image_scale_keys, cs="pixel", pdf_file=pdf)


    def load_segmentation_all(self, 
                              exp_name:str,
                            seg_dir:str,
                            type:str="vpt", 
                            prefix_name:str="default",
                            plot:bool=False,
                            ):
        """
        Load segmentation data for all regions of an experiment into spatialdata objects.

        Parameters:
        exp_name (str): Name of the experiment.
        seg_dir (str): Directory containing the segmentation data.
        type (str): Type of the segmentation data to load (default is "vpt").
        prefix_name (str): Prefix for the keys in the spatialdata object (default is "default").
        """

        regions = glob.glob(f"{seg_dir}/region_*")
        for reg_dir in regions: 
            reg_name = reg_dir.split("/")[-1]
            self.load_segmentation_region(exp_name, reg_name, reg_dir, type=type, prefix_name=prefix_name, plot=plot)

if __name__ == "__main__": 
    fire.Fire(io_cli)
    # fire.Fire({"ingest_region": ingest_region,
    #            "ingest_all": ingest_asll,
    #            "load_segmentation_region" : load_segmentation_region,
    #            "load_segmentation_all" : load_segmentation_all,
    # })
