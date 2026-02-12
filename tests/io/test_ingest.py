# import os
# import sys
# import warnings

# with warnings.catch_warnings():
#     warnings.filterwarnings("ignore")
#     import spatialdata as sd

# spida_path = "/ceph/cephatlas/aklein/spida/src"
# sys.path.append(spida_path)
# from spida.io.ingest_exp import read_merscope
# from spida.plotting.ingest import plot_images, plot_shapes, plot_points, plot_overlap

# from matplotlib.backends.backend_pdf import PdfPages


# # parameters
# EXPERIMENT = "202503281042_BICAN4x1A9-E-01_VMSC31910"
# REGION = "region_UCI-4723-A9-E1"

# prefix_name = "default"

# IMAGE_KEY = f"{prefix_name}_{EXPERIMENT}_{REGION}_z3"
# POINTS_KEY = f"{prefix_name}_{EXPERIMENT}_{REGION}_transcripts"
# SHAPES_KEY = f"{prefix_name}_{EXPERIMENT}_{REGION}_polygons"

# processed_path = "/data/aklein/bican_zarr"
# root_path = "/ceph/cephatlas/merscope_data/processed"

# input_path = f"{root_path}/{EXPERIMENT}/out/{REGION}"
# zarr_path = f"{processed_path}/{EXPERIMENT}/{REGION}"


# image_path = f"/ceph/cephatlas/aklein/bican/images/{EXPERIMENT}/pixi-test.pdf"

# print("Reading Merscope data...")
# print("EXPERIMENT:", EXPERIMENT)
# print("REGION:", REGION)
# print("zarr_path:", zarr_path)


# sdata = read_merscope(
#     input_path, zarr_path, exp_name=EXPERIMENT, reg_name=REGION, prefix_name=prefix_name
# )

# image_channels = sd.models.get_channel_names(sdata[IMAGE_KEY])
# image_scale_keys = list(sdata[IMAGE_KEY].keys())

# with PdfPages(image_path) as pdf:
#     plot_images(
#         sdata, IMAGE_KEY, image_scale_keys, image_channels, cs="global", pdf_file=pdf
#     )
#     plot_shapes(sdata, SHAPES_KEY, table_name="table", cs="pixel", pdf_file=pdf)
#     plot_points(sdata, POINTS_KEY, cs="pixel", cmap="tab10", pdf_file=pdf)
#     plot_overlap(
#         sdata,
#         IMAGE_KEY,
#         SHAPES_KEY,
#         POINTS_KEY,
#         image_scale_keys,
#         cs="pixel",
#         pdf_file=pdf,
#     )
