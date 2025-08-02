from .cluster_plots import plot_scatter
from .ingest import (
    plot_images,
    plot_shapes,
    plot_points,
    plot_overlap,
    plot_seg_load,
    plot_example_doublets
)
from .P_plots import (
    plot_filtering,
    plot_setup,
    plot_doublets,
    plot_resolvi,
    plot_dataset
)
from ._scatteplot import plot_categorical, plot_continuous
from ._spatial import plot_spatial_continuous, plot_spatial_categorical
from .qc_plots import plot_violin_QC
