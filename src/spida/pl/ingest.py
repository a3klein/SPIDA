import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
plt.rcParams['axes.facecolor'] = 'black'

import spatialdata_plot


 
def plot_images(sdata, IMAGE_KEY, image_scale_keys, image_channels, cs:str="global", pdf_file:str=None): 
    max_int = sdata[IMAGE_KEY][image_scale_keys[-1]]['image'].max(['x', 'y']).compute().to_dataframe().to_dict()['image']
    min_int = sdata[IMAGE_KEY][image_scale_keys[-1]]['image'].min(['x', 'y']).compute().to_dataframe().to_dict()['image']


    fig, axes = plt.subplots(1, len(image_channels), figsize=(20, 8))

    for i, channel in enumerate(image_channels):
        norm = Normalize(vmin=min_int[channel], vmax=max_int[channel]*0.5)
        sdata.pl.render_images(IMAGE_KEY, channel=channel, cmap="grey", norm=norm).pl.show(ax=axes[i], title=channel, coordinate_systems=cs)

    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches='tight')
    else: 
        fig.show()

def plot_shapes(sdata, SHAPES_KEY, table_name:str="table", cs:str="pixel", pdf_file:str=None):
    fig, axes = plt.subplots(1, 3, figsize=(20, 8))

    sdata.pl.render_shapes(SHAPES_KEY, table_name=table_name, method="matplotlib").pl.show(ax=axes[0], title="all shapes", coordinate_systems=cs)
    sdata.pl.render_shapes(SHAPES_KEY, table_name=table_name, scale=2, method="matplotlib", color="transcript_count", cmap="viridis").pl.show(ax=axes[1], title="Transcript Count", coordinate_systems=cs)
    sdata.pl.render_shapes(SHAPES_KEY, table_name=table_name, scale=2, method="matplotlib", color="volume", cmap="viridis").pl.show(ax=axes[2], title="Volume", coordinate_systems=cs)

    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches='tight')
    else: 
        fig.show()


def plot_points(sdata, POINTS_KEY, TABLE_KEY, cs:str="pixel", cmap:str="tab10", pdf_file:str=None): 

    fig, axes = plt.subplots(1, 3, figsize=(20, 8))
    for i in range(3): 
        gene_ex = sdata[TABLE_KEY].var.index[i]
        sdata.pl.render_points(POINTS_KEY, groups=gene_ex, color="gene", cmap=cmap).pl.show(ax=axes[i], title=gene_ex, coordinate_systems=cs)
    
    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches='tight')
    else: 
        fig.show()

def plot_overlap(sdata, IMAGE_KEY, SHAPES_KEY, POINTS_KEY, TABLE_KEY, image_scale_keys, cs:str="pixel", pdf_file:str=None):
    """
    Plot the overlap of shapes and points on the image.
    """
    gene_ex = sdata[TABLE_KEY].var.index[0]
    
    # Render images with shapes and points
    fig, axes = plt.subplots(1, 3, figsize=(20, 8))

    sdata.pl.render_images(IMAGE_KEY, scale=image_scale_keys[-1]).pl.render_shapes(SHAPES_KEY).pl.show(ax=axes[0], title="Shapes on Image", coordinate_systems=cs)
    sdata.pl.render_images(IMAGE_KEY, scale=image_scale_keys[-1]).pl.render_points(POINTS_KEY, groups=gene_ex, color="gene", cmap="tab10").pl.show(ax=axes[1], title="Points on Image", coordinate_systems=cs)
    sdata.pl.render_shapes(SHAPES_KEY).pl.render_points(POINTS_KEY, groups=gene_ex, color="gene", cmap="tab10").pl.show(ax=axes[2], title="Points on Shapes", coordinate_systems=cs)

    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches='tight')
    else: 
        fig.show()
    
