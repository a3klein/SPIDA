"""
MapMyCells annotation: MMCConfig, MMCPreprocessor, MMCAnnotator.

AnnData-centric: nothing here imports spatialdata.
"""

import json
import glob
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

logger = logging.getLogger(__name__)


# ============================================================================
# Configuration Management
# ============================================================================

@dataclass
class MMCConfig:
    """MapMyCells run settings plus the two MMC-owned paths."""

    mmc_store_path: Path
    gene_panel_path: Path | None = None

    n_cpu: int = 8
    bootstrap_factor: float = 0.8
    bootstrap_iterations: int = 100
    rng_seed: int = 13
    ref_norm: str = "log2CPM"
    query_norm: str = "log2CPM"

    def __post_init__(self):
        self.mmc_store_path = Path(self.mmc_store_path)
        if self.gene_panel_path is not None:
            self.gene_panel_path = Path(self.gene_panel_path)

    @classmethod
    def from_spida(cls, config, **settings) -> "MMCConfig":
        """Paths from a SpidaConfig (``mmc_store_path``, ``gene_panel_path``); settings as given."""
        return cls(
            mmc_store_path=settings.pop("mmc_store_path", None) or config.require("mmc_store_path"),
            gene_panel_path=settings.pop("gene_panel_path", None) or config.gene_panel_path,
            **settings,
        )

    def get_mmc_store_dir(self, identifier: str) -> Path:
        """Get or create the MMC store directory for a given identifier."""
        store_dir = self.mmc_store_path / identifier
        store_dir.mkdir(parents=True, exist_ok=True)
        return store_dir
    
    def get_marker_dir(self, identifier: str) -> Path:
        """Get marker directory."""
        marker_dir = self.get_mmc_store_dir(identifier) / "marker_dir"
        marker_dir.mkdir(parents=True, exist_ok=True)
        return marker_dir
    
    def get_precomp_stats_path(
        self, identifier: str, ref_name: str
    ) -> Path:
        """Get precomputed stats file path."""
        store_dir = self.get_mmc_store_dir(identifier)
        return store_dir / f"{ref_name}_precomputed_stats.h5"
    
    def get_query_path(self, identifier: str) -> Path:
        """Get query AnnData path."""
        store_dir = self.get_mmc_store_dir(identifier)
        return store_dir / "query.h5ad"
    
    def _gene_panel(self) -> Path:
        if self.gene_panel_path is None:
            raise FileNotFoundError("no gene_panel_path; pass the gene files or set GENE_PANEL_PATH in the config")
        return self.gene_panel_path

    def get_gene_names_path(self, identifier: str) -> Path:
        """Find gene_names.txt file for identifier."""
        self._gene_panel()
        matches = glob.glob(str(self.gene_panel_path / f"*{identifier}*gene_names.txt"))
        if not matches:
            raise FileNotFoundError(
                f"No gene_names.txt found matching '{identifier}' in {self.gene_panel_path}"
            )
        return Path(matches[0])
    
    def get_gene_name_mapping_path(self, identifier: str) -> Path:
        """Find gene_name_mapping.json file for identifier."""
        self._gene_panel()
        matches = glob.glob(str(self.gene_panel_path / f"*{identifier}*gene_name_mapping.json"))
        if not matches:
            raise FileNotFoundError(
                f"No gene_name_mapping.json found matching '{identifier}' in {self.gene_panel_path}"
            )
        return Path(matches[0])
    
    def validate(self) -> None:
        """Configured paths must exist; gene_panel_path is optional."""
        paths = {"mmc_store_path": self.mmc_store_path, "gene_panel_path": self.gene_panel_path}
        for name, path in paths.items():
            logger.info(f"{name}: {path}")
            if path is not None and not path.exists():
                raise FileNotFoundError(f"{name} does not exist: {path}")


# ============================================================================
# Preprocessing Component
# ============================================================================

class MMCPreprocessor:
    """Handles all reference preprocessing for MapMyCells annotation."""
    
    def __init__(self, config: MMCConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def load_gene_name_mapping(self, identifier: str) -> Dict[str, str]:
        """
        Load gene name mapping from JSON file.
        
        Parameters:
        -----------
        identifier : str
            Reference identifier (e.g., 'mouse_motor_cortex')
            
        Returns:
        --------
        dict mapping reference gene names to MERFISH gene names
        """
        try:
            mapping_path = self.config.get_gene_name_mapping_path(identifier)
            with open(mapping_path, 'r') as f:
                mapping = json.load(f)
            self.logger.info(f"Loaded gene name mapping from {mapping_path}")
            self.logger.info(f"Mapping contains {len(mapping)} genes")
            return mapping
        except FileNotFoundError:
            self.logger.warning(f"No gene name mapping found for {identifier}, using identity mapping")
            return {}
    
    def apply_gene_name_mapping(
        self,
        ref_adata: ad.AnnData,
        mapping: Dict[str, str]
    ) -> ad.AnnData:
        """
        Apply gene name mapping to reference AnnData.
        
        Maps reference gene names to MERFISH gene names using the provided mapping.
        If a gene is not in the mapping, it retains its original name.
        
        Parameters:
        -----------
        ref_adata : ad.AnnData
            Reference AnnData with unmapped gene names
        mapping : dict
            Dictionary mapping reference gene names to MERFISH gene names
            
        Returns:
        --------
        ad.AnnData with mapped gene names
        """
        if not mapping:
            self.logger.info("No gene name mapping provided, using original names")
            return ref_adata
        
        # Create copy to avoid modifying original
        ref_adata = ref_adata.copy()
        
        # Map gene names
        original_names = ref_adata.var_names.tolist()
        mapped_names = [mapping.get(name, name) for name in original_names]
        
        ref_adata.var_names = mapped_names
        
        unmapped_count = sum(1 for orig, mapped in zip(original_names, mapped_names) if orig == mapped)
        mapped_count = len(original_names) - unmapped_count
        
        self.logger.info(f"Applied gene name mapping: {mapped_count} genes mapped, {unmapped_count} unchanged")
        
        return ref_adata
    
    def load_gene_names(self, identifier: str) -> List[str]:
        """
        Load gene names from gene_names.txt file.
        
        Parameters:
        -----------
        identifier : str
            Reference identifier
            
        Returns:
        --------
        list of gene names to use
        """
        try:
            gene_names_path = self.config.get_gene_names_path(identifier)
            with open(gene_names_path, 'r') as f:
                gene_names = [line.strip() for line in f if line.strip()]
            self.logger.info(f"Loaded {len(gene_names)} gene names from {gene_names_path}")
            return gene_names
        except FileNotFoundError:
            raise FileNotFoundError(
                f"No gene_names.txt found for identifier '{identifier}' in {self.config.gene_panel_path}"
            )
    
    def create_empty_query_adata(
        self,
        ref_path: str | Path,
        gene_list: List[str],
        output_path: str | Path,
    ) -> None:
        """
        Create an empty AnnData with reference genes intersected with query genes.
        
        This is used by cell_type_mapper to identify which genes to use.
        """
        self.logger.info(f"Creating empty query AnnData at {output_path}")
        
        ref_adata = ad.read_h5ad(ref_path, backed="r")
        var_index = ref_adata.var.index
        
        # Intersect with available genes
        query_index = var_index.intersection(pd.Index(gene_list))
        self.logger.info(f"Gene intersection: {len(query_index)} / {len(gene_list)} genes found")
        
        # Create minimal AnnData
        query = ad.AnnData(
            X=sparse.csr_matrix(np.zeros((1, len(query_index)))),
            var=pd.DataFrame(index=query_index),
        )
        
        query.write_h5ad(output_path)
        self.logger.info(f"Created query AnnData with {len(query_index)} genes")
    
    def precompute_reference_stats(
        self,
        ref_path: str | Path,
        hierarchy_list: List[str],
        identifier: str,
    ) -> Path:
        """Precompute statistics for the reference AnnData."""
        from cell_type_mapper.cli.precompute_stats_scrattch import (  # type: ignore
            PrecomputationScrattchRunner,
        )
        
        ref_name = Path(ref_path).stem
        output_path = self.config.get_precomp_stats_path(identifier, ref_name)
        
        self.logger.info(f"Precomputing reference stats -> {output_path}")
        
        precomputation_config = {
            "hierarchy": hierarchy_list,
            "h5ad_path": str(ref_path),
            "output_path": str(output_path),
            "n_processors": self.config.n_cpu,
            "normalization": self.config.ref_norm,
            "clobber": True,
        }
        self.logger.info(f"Precomputation config: {precomputation_config}")
        
        runner = PrecomputationScrattchRunner(args=[], input_data=precomputation_config)
        runner.run()
        
        self.logger.info("Precomputation complete")
        return output_path
    
    def compute_reference_markers(
        self,
        precomp_path: str | Path,
        identifier: str,
        query_path: str | Path = None,
        n_valid: int = 10,
    ) -> Path:
        """Compute reference markers for the hierarchy."""
        from cell_type_mapper.cli.reference_markers import (  # type: ignore
            ReferenceMarkerRunner,
        )
        
        marker_dir = self.config.get_marker_dir(identifier)
        
        self.logger.info(f"Computing reference markers -> {marker_dir}")
        
        reference_config = {
            "precomputed_path_list": [str(precomp_path)],
            "n_valid": n_valid,
            "output_dir": str(marker_dir),
            "clobber": True,
            "query_path": str(query_path) if query_path else None,
        }
        self.logger.info(f"Reference marker config: {reference_config}")
        
        runner = ReferenceMarkerRunner(args=[], input_data=reference_config)
        runner.run()
        
        markers_path = marker_dir / "reference_markers.h5"
        self.logger.info(f"Reference markers saved to {markers_path}")
        return markers_path
    
    def compute_query_markers(
        self,
        ref_markers_path: str | Path,
        identifier: str,
        query_path: str | Path,
        n_per_utility: int = 10,
    ) -> Path:
        """Downsample markers to those present in query."""
        from cell_type_mapper.cli.query_markers import (  # type: ignore
            QueryMarkerRunner,
        )
        
        marker_dir = self.config.get_marker_dir(identifier)
        output_path = marker_dir / "calc_markers.json"
        
        self.logger.info(f"Computing query markers -> {output_path}")
        
        query_config = {
            "output_path": str(output_path),
            "reference_marker_path_list": [str(ref_markers_path)],
            "n_per_utility": n_per_utility,
            "n_processors": self.config.n_cpu,
            "query_path": str(query_path),
        }
        self.logger.info(f"Query marker config: {query_config}")

        runner = QueryMarkerRunner(args=[], input_data=query_config)
        runner.run()
        
        self.logger.info(f"Query markers saved to {output_path}")
        return output_path
    
    def setup_annotation_pipeline(
        self,
        ref_path: str | Path,
        identifier: str,
        hierarchy_list: List[str],
        gene_names_path: str | Path = None,
        gene_name_mapping_path: str | Path = None,
        **preprocessing_params,
    ) -> Dict[str, Path]:
        """
        Complete setup: creates query AnnData, precomputes stats, and computes markers.
        
        Parameters:
        -----------
        ref_path : str or Path
            Path to reference AnnData
        identifier : str
            Unique identifier for this reference/codebook combination
        hierarchy_list : List[str]
            Hierarchy levels for annotation
        gene_names_path : str or Path, optional
            Path to gene_names.txt file. If not provided, will search using identifier
        gene_name_mapping_path : str or Path, optional
            Path to gene_name_mapping.json file. If not provided, will search using identifier
        **preprocessing_params
            Additional params (n_valid, n_per_utility, etc.)
        
        Returns:
        --------
        dict with paths to precomputed files
        """
        self.logger.info("Starting annotation pipeline setup")
        self.logger.info(f"  Identifier: {identifier}")
        self.logger.info(f"  Hierarchy: {hierarchy_list}")
        
        # Load gene name mapping
        if gene_name_mapping_path:
            with open(gene_name_mapping_path, 'r') as f:
                mapping = json.load(f)
            self.logger.info(f"Loaded gene name mapping from {gene_name_mapping_path}")
        else:
            mapping = self.load_gene_name_mapping(identifier)
        
        # Load and apply gene name mapping to reference
        ref_adata = ad.read_h5ad(ref_path)
        ref_adata = self.apply_gene_name_mapping(ref_adata, mapping)
        
        # Create temporary reference file with mapped names
        temp_ref_path = self.config.get_mmc_store_dir(identifier) / "temp_ref_mapped.h5ad"
        ref_adata.write_h5ad(temp_ref_path)
        self.logger.info(f"Saved reference with mapped gene names to {temp_ref_path}")
        
        # Load gene names from gene_names.txt
        if gene_names_path:
            with open(gene_names_path, 'r') as f:
                gene_list = [line.strip() for line in f if line.strip()]
            self.logger.info(f"Loaded {len(gene_list)} gene names from {gene_names_path}")
        else:
            gene_list = self.load_gene_names(identifier)
        
        # Create query AnnData
        query_path = self.config.get_query_path(identifier)
        self.create_empty_query_adata(str(temp_ref_path), gene_list, query_path)
        
        # Precompute reference stats
        precomp_stats_path = self.precompute_reference_stats(
            str(temp_ref_path), hierarchy_list, identifier
        )
        
        # Compute reference markers
        n_valid = preprocessing_params.get("n_valid", 10)
        ref_markers_path = self.compute_reference_markers(
            precomp_stats_path, identifier, query_path, n_valid
        )
        
        # Compute query markers
        n_per_utility = preprocessing_params.get("n_per_utility", 10)
        query_markers_path = self.compute_query_markers(
            ref_markers_path, identifier, query_path, n_per_utility
        )
        
        result = {
            "precomp_stats_path": precomp_stats_path,
            "ref_markers_path": ref_markers_path,
            "query_markers_path": query_markers_path,
            "query_adata_path": query_path,
        }
        
        self.logger.info("Annotation pipeline setup complete")
        for _key, _path in result.items():
            self.logger.info(f"  {_key}: {_path}")
        
        return result


# ============================================================================
# Annotation Component
# ============================================================================

class MMCAnnotator:
    """Handles query annotation with MapMyCells."""
    
    def __init__(self, config: MMCConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def annotate(
        self,
        query_adata_path: str | Path,
        identifier: str,
        output_dir: str | Path,
        **annotation_params,
    ) -> Tuple[Path, Path]:
        """Run MapMyCells annotation on query data."""
        from cell_type_mapper.cli.from_specified_markers import (  # type: ignore
            FromSpecifiedMarkersRunner,
        )
        
        # Get paths
        marker_dir = self.config.get_marker_dir(identifier)
        markers_path = marker_dir / "calc_markers.json"
        mmc_store_dir = self.config.get_mmc_store_dir(identifier)
        precomp_path = list(mmc_store_dir.glob("*_precomputed_stats.h5"))[0]
        # ref_name = list(marker_dir.glob("*.h5"))[0] # .stem.replace("_reference_markers", "")
        # precomp_path = self.config.get_precomp_stats_path(brain_region, codebook, ref_name)
        
        # Output paths
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        extended_results_path = output_dir / "extended_results.json"
        csv_results_path = output_dir / "csv_results.csv"
        
        self.logger.info(f"Running annotation for {query_adata_path}")
        self.logger.info(f"  Output: {output_dir}")
        
        # Annotation parameters
        bootstrap_factor = annotation_params.get("bootstrap_factor", self.config.bootstrap_factor)
        bootstrap_iterations = annotation_params.get("bootstrap_iterations", self.config.bootstrap_iterations)
        rng_seed = annotation_params.get("rng_seed", self.config.rng_seed)
        n_cpu = annotation_params.get("n_cpu", self.config.n_cpu)
        
        config = {
            "precomputed_stats": {"path": str(precomp_path)},
            "query_markers": {"serialized_lookup": str(markers_path)},
            "type_assignment": {
                "bootstrap_factor": bootstrap_factor,
                "bootstrap_iteration": bootstrap_iterations,
                "rng_seed": rng_seed,
                "n_processors": n_cpu,
                "normalization": self.config.query_norm,
            },
            "query_path": str(query_adata_path),
            "extended_result_path": str(extended_results_path),
            "csv_result_path": str(csv_results_path),
        }
        
        self.logger.info(f"Annotation config: {config}")
        
        runner = FromSpecifiedMarkersRunner(args=[], input_data=config)
        runner.run()
        
        self.logger.info(f"Annotation complete. Results at {output_dir}")
        return extended_results_path, csv_results_path
    
    def transfer_labels(
        self,
        adata: ad.AnnData,
        csv_results_path: str | Path,
        annot_name: str = "mmc",
    ) -> ad.AnnData:
        """Transfer labels to ``{annot_name}_{level}`` / ``{annot_name}_{level}_transfer_score``."""
        self.logger.info(f"Transferring labels from {csv_results_path}")
        
        results_df = pd.read_csv(csv_results_path, comment="#", index_col=0)
        
        # Extract hierarchy levels from column names
        hierarchy_list = []
        for col in results_df.columns:
            if col.endswith("_bootstrapping_probability"):
                hierarchy_list.append(col.split("_bootstrapping_probability")[0])
        
        self.logger.info(f"Found hierarchy levels: {hierarchy_list}")
        
        #TODO: a more elegant solution to align the index types maybe? 
        adata.obs.index = adata.obs.index.astype(str)
        results_df.index = results_df.index.astype(str)

        # Add columns to obs
        for h in hierarchy_list:
            score_col = f"{h}_bootstrapping_probability"
            label_col = f"{h}_name"
            
            adata.obs[f"{annot_name}_{h}_transfer_score"] = results_df[score_col]
            adata.obs[f"{annot_name}_{h}"] = results_df[label_col]

        # Log what was added
        added_cols = [f"{annot_name}_{h}" for h in hierarchy_list] + [f"{annot_name}_{h}_transfer_score" for h in hierarchy_list]
        self.logger.info(f"Added columns to adata.obs: {added_cols}")
        
        return adata


# ============================================================================
# High-Level Public API (New Architecture)
# ============================================================================

def setup_annotation_pipeline(
    ref_path: str | Path,
    identifier: str,
    hierarchy_list: List[str],
    config: MMCConfig,
    gene_names_path: str | Path = None,
    gene_name_mapping_path: str | Path = None,
    **preprocessing_params,
) -> Dict[str, Path]:
    """
    Complete setup for MapMyCells annotation.
    This creates all necessary reference data and markers for annotation.
    Should be run once per reference/identifier combination.

    Parameters:
        ref_path: Path to reference single-cell RNA-seq AnnData.
        identifier: Unique identifier for this reference (e.g. 'mouse_motor_cortex').
        hierarchy_list: Hierarchy levels for annotation.
        config: MMCConfig (MMC store, gene panel dir, run settings).
        gene_names_path: gene_names.txt; if not provided, searched in config.gene_panel_path by identifier.
        gene_name_mapping_path: gene_name_mapping.json; if not provided, searched the same way.
        **preprocessing_params: Additional params (n_valid, n_per_utility, etc.).

    Returns:
        dict with paths to precomputed files
    """
    config.validate()
    preprocessor = MMCPreprocessor(config)
    return preprocessor.setup_annotation_pipeline(
        ref_path, identifier, hierarchy_list, gene_names_path, gene_name_mapping_path, **preprocessing_params
    )


def annotate_adata(
    query_adata: ad.AnnData,
    identifier: str,
    config: MMCConfig,
    query_path: str | Path,
    output_dir: str | Path,
    filter_annot: bool = False,
    plot_path: str | Path | None = None,
    palette_path: str | Path | None = None,
    **annotation_params,
) -> ad.AnnData:
    """
    Annotate a query AnnData with MapMyCells; labels go to ``mmc_<identifier>_<level>`` columns.
    MapMyCells reads its query from disk, so the AnnData is written to ``query_path``,
    annotated, and written back there.

    Parameters:
        query_adata: Query data to annotate.
        identifier: Reference identifier.
        config: MMCConfig (MMC store, run settings).
        query_path: AnnData file the query is written to and annotated in.
        output_dir: Directory for the raw MapMyCells results.
        filter_annot: Whether to apply filtering to annotations based on scores and counts.
        plot_path: pdf for the annotation plots (None: no plots).
        palette_path: Path to color palette for plotting.
        **annotation_params: Additional params (bootstrap_factor, rng_seed, filter thresholds, etc.).

    Returns:
        Annotated AnnData with MMC labels
    """
    config.validate()
    annot = f"mmc_{identifier}"
    query_path = Path(query_path)
    query_path.parent.mkdir(parents=True, exist_ok=True)
    query_adata.write_h5ad(query_path)

    annotator = MMCAnnotator(config)
    _, csv_path = annotator.annotate(query_path, identifier, output_dir, **annotation_params)
    query_adata = annotator.transfer_labels(query_adata, csv_path, annot_name=annot)

    hierarchy_list = [
        col[len(annot) + 1:-len("_transfer_score")]
        for col in query_adata.obs.columns
        if col.startswith(f"{annot}_") and col.endswith("_transfer_score")
    ]
    if filter_annot:
        query_adata.obs[f'{annot}_total_filt'] = filt_adata_mmc(query_adata, hierarchy_list, annot_name=annot, score_ext="transfer_score", **annotation_params)
    for col in hierarchy_list:
        query_adata.obs[f'{annot}_{col}'] = query_adata.obs[f'{annot}_{col}'].astype(str).astype("category")
        query_adata.obs[f'{annot}_{col}_transfer_score'] = query_adata.obs[f'{annot}_{col}_transfer_score'].astype(float).fillna(0)
    query_adata.write_h5ad(query_path)

    if plot_path is not None:
        logger.info("Generating annotation plots")
        plot_path = Path(plot_path)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        from spida.pl.annot_plots import plot_annot
        plot_annot(
            adata = query_adata,
            hierarchy_list = hierarchy_list,
            annot_name = annot,
            palette_path = palette_path,
            pdf_path = plot_path,
            plot_filt = filter_annot,
            filt_column = "total_filt" if filter_annot else None,
            score_ext = "transfer_score",
            umap_col = "base_umap",
            single_clust_plot_thr = 40,
        )
    return query_adata


def annotate_region(
    dataset,
    prefix_name: str,
    suffix: str,
    identifier: str,
    config: MMCConfig,
    adata_file: str | Path | None = None,
    refresh_table: bool = False,
    filter_annot: bool = False,
    plot: bool = False,
    palette_path: str | Path = None,
    plot_name: str = "mmc_annotation_plots",
    **annotation_params,
) -> ad.AnnData:
    """
    Annotate a specific experiment region: segmentation ``prefix_name`` (table suffix ``suffix``).
    Labels are added to the region's AnnData file (or ``adata_file``); MapMyCells results go to
    ``<annotations>/<table>/mmc/<identifier>``, figures to ``<figures>/<prefix>/<plot_name>_<identifier>.pdf``.

    Parameters:
        dataset: The region's locations (Dataset).
        prefix_name: Segmentation name, e.g. cellpose_cell.
        suffix: Table suffix, e.g. _filt.
        identifier: Reference identifier.
        config: MMCConfig (MMC store, run settings).
        adata_file: Annotate this AnnData file instead; labels are written back to it.
        refresh_table: Re-read the table from the zarr store even if its AnnData file exists.
        filter_annot: Whether to apply filtering to annotations based on scores and counts.
        plot: Whether to generate annotation plots.
        palette_path: Path to color palette for plotting.
        **annotation_params: Additional params for annotation.

    Returns:
        Annotated AnnData with MMC labels
    """
    from spida.I._io import load_table
    from spida.io.layout import table_key

    output_dir = dataset.annotation_dir(table_key(prefix_name, suffix), "mmc", identifier)
    plot_path = dataset.figures_dir(prefix_name) / f"{plot_name}_{identifier}.pdf" if plot else None
    adata, query_path = load_table(dataset, prefix_name, suffix, adata_file=adata_file, refresh=refresh_table)
    return annotate_adata(
        adata,
        identifier,
        config,
        query_path,
        output_dir,
        filter_annot=filter_annot,
        plot_path=plot_path,
        palette_path=palette_path,
        **annotation_params,
    )


def annotate_experiment(
    spida_config,
    exp_name: str,
    prefix_name: str,
    suffix: str,
    identifier: str,
    config: MMCConfig,
    **kwargs,
) -> int:
    """
    Annotate all regions in an experiment (``annotate_region`` on each region in the zarr store).

    Parameters:
        spida_config: SpidaConfig (store roots).
        exp_name: Experiment name.
        config: MMCConfig (MMC store, run settings).
        **kwargs: Passed to annotate_region.

    Returns:
        0 on success
    """
    from spida.io.layout import TreeLayout

    layout = TreeLayout(spida_config)
    for reg in layout.units(exp_name):
        logger.info(f"Processing region: {reg}")
        annotate_region(layout.dataset(exp_name, reg), prefix_name, suffix, identifier, config, **kwargs)
    return 0


def filt_adata_mmc(
    adata : ad.AnnData,
    hierarchy_list,
    annot_name : str = "mmc",
    score_ext : str | None = "transfer_score",
    **kwargs
):
    num_l = len(hierarchy_list)
    df_annot = adata.obs[[f"{annot_name}_{_cc}" for _cc in hierarchy_list] + [f"{annot_name}_{_cc}_{score_ext}" for _cc in hierarchy_list]].copy()
    # TODO: Change logic here
    prev_cc = None
    for i, _cc in enumerate(hierarchy_list):
        l_thr = kwargs.get(f"l{i+1}_thr", kwargs.get(f"{_cc}_thr", 0.5))
        l_count_thr = kwargs.get(f"l{i+1}_count_thr", kwargs.get(f"{_cc}_count_thr", 0))
        df_annot[f"{annot_name}_{_cc}_filt"] = False
        if prev_cc is None: 
            df_annot.loc[df_annot[f"{annot_name}_{_cc}_{score_ext}"] > l_thr, f"{annot_name}_{_cc}_filt"] = True
        else: 
            df_annot.loc[
                (df_annot[f"{annot_name}_{_cc}_{score_ext}"] > l_thr) & 
                (df_annot[f"{annot_name}_{prev_cc}_filt"]),
                f"{annot_name}_{_cc}_filt"
            ] = True
        prev_cc = _cc
        vc = df_annot.loc[df_annot[f"{annot_name}_{_cc}_filt"] == True][f"{annot_name}_{_cc}"].value_counts()
        remove = vc.loc[vc < l_count_thr].index
        df_annot.loc[df_annot[f"{annot_name}_{_cc}"].isin(remove), f"{annot_name}_{_cc}_filt"] = False
    df_annot['total_filt'] = df_annot[[f"{annot_name}_{_cc}_filt" for _cc in hierarchy_list]].all(axis=1)
    return df_annot['total_filt']
