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

from .I_plots import (
    plot_allcools_joint_embeddings,
    plot_allcools_spatial_annot
)
from ._scatteplot import plot_categorical, plot_continuous, categorical_scatter, continuous_scatter
from ._spatial import plot_spatial_continuous, plot_spatial_categorical
from .qc_plots import plot_violin_QC
from .palettes import add_colors
# from .palettes import register_spatial_cmaps
from .transcript_plots import plot_hex_qc, plot_hex_clusters
