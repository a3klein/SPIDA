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


    ncols = 4
    nrows = int(len(image_channels) / ncols) + (len(image_channels) % ncols > 0)
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows*5), dpi=200)
    axes = axes.flatten()

    for i, channel in enumerate(image_channels):
        norm = Normalize(vmin=min_int[channel], vmax=max_int[channel]*0.5)
        sdata.pl.render_images(IMAGE_KEY, channel=channel, cmap="grey", norm=norm).pl.show(ax=axes[i], title=channel, coordinate_systems=cs, colorbar=False)
        
    plt.tight_layout()

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



def plot_seg_load(sdata, image_key, shapes_key, cs="global", pdf_file=None): 
    from ..utilities.boxes import generate_multiple_random_boxes, visualize_boxes

    # Subsetting for faster cropping
    sub_sdata = sdata.subset([image_key, shapes_key])
    
    # Generating image metadata
    image_scale_keys = list(sub_sdata[image_key].keys())
    max_int = sub_sdata[image_key][image_scale_keys[-1]]['image'].max(['x', 'y']).compute().to_dataframe().to_dict()['image']
    min_int = sub_sdata[image_key][image_scale_keys[-1]]['image'].min(['x', 'y']).compute().to_dataframe().to_dict()['image']
    norm = Normalize(vmin=min_int["DAPI"], vmax=max_int["DAPI"]*0.5)

    # generating boxes
    c, frame_height, frame_width = sub_sdata[image_key][image_scale_keys[0]]['image'].shape
    boxes = generate_multiple_random_boxes(frame_width, frame_height, 1500, 1500, num_boxes=20, avoid_overlap=True)
    
    # plotting
    fig, axes = plt.subplots(2, 3, figsize=(20, 10), dpi=200)
    axes = axes.flatten()
    # First plot is the boxes overlaid on the DAPI image
    pp = sub_sdata.pl.render_images(image_key, channel="DAPI", cmap="grey", norm=norm).pl.show(ax=axes[0], title="DAPI", coordinate_systems=cs, colorbar=False)
    visualize_boxes(frame_width=frame_width, frame_height=frame_height, boxes=boxes, ax=axes[0], show_plot=False)

    # Iterating for 5 example fovs
    plot_counter = 1
    for b_cont, _b in enumerate(boxes):
        if plot_counter == 6: 
            continue
        xmin, ymin, width, height = _b
        xmax = xmin + width
        ymax = ymin + height
        
        try: 
            cropped_sdata = sub_sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=[xmin, ymin],
                max_coordinate=[xmax, ymax],
                target_coordinate_system=cs,
            )
        except Exception as e:
            print(f"Error cropping data: {e}")
            continue
        if shapes_key not in cropped_sdata:
            print(f"Shapes key {shapes_key} not found in cropped_sdata")
            continue
        
        ax = axes[plot_counter]
        plot_counter += 1
        pp = cropped_sdata.pl.render_images(image_key, channel="DAPI", cmap="grey", norm=norm)
        pp = pp.pl.render_shapes(shapes_key, fill_alpha=0, cmap="jet", outline_alpha=1, outline_color="red")
        pp.pl.show(ax=ax, title="Shapes", coordinate_systems=cs, colorbar=False)
        ax.axis('off')
        ax.set_title("Shapes - " + str(b_cont+1))
    
    if pdf_file: 
        pdf_file.savefig(fig, bbox_inches='tight')
    else: 
        fig.show()

    
