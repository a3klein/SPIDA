import os
from dotenv import load_dotenv  # type: ignore
from pathlib import Path
import glob
import logging

load_dotenv()
logger = logging.getLogger(__package__)


# The controller that calls the specified segmentation algorithms
def run_segmentation(
    type: str,
    exp_name: str,
    reg_name: str,
    input_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    root_dir: str | Path | None = None,  
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    Run an implemented segmentation algorithm on a given region

    Parameters:
    type (str): Type of segmentation to run (e.g., "proseg", "vpt", "cellpose", "mesmer").
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    input_dir (Path): Directory containing the input data.
    output_dir (Path): Directory to save the output data.
    root_dir (Path): Root directory for the experiment data, used for reading necessary files for segmentation (sometimes same as input_dir).
    root_path (Path, optional): Root path for the data (default is None, uses environment variable).
    segmentation_store (Path, optional): Default Path to store the segmentation results (default is None, uses environment variable).
    vpt_bin_path (Path, optional): Path to the VPT binary (default is None, uses environment variable).
    rust_bin_path (Path, optional): Path to the Rust binary for Proseg (default is None, uses environment variable).
    **kwargs: Additional keyword arguments to pass to the segmentation function.
    """

    if root_path is None: 
        processed_root_path = os.getenv("PROCESSED_ROOT_PATH")
    else:
        processed_root_path = root_path
    if input_dir is None:
        input_dir = f"{processed_root_path}/{exp_name}/out"
    if root_dir is None: 
        root_dir = f"{processed_root_path}/{exp_name}/out"
    if output_dir is None:
        if segmentation_store is None:
            seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        else:
            seg_out_path = segmentation_store
        output_dir = f"{seg_out_path}/{exp_name}/{type}"

    logger.info(
        f"Running segmentation of type {type} on region {reg_name} in experiment {exp_name}."
    )
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    if type == "proseg":
        from .backends.proseg import run_proseg, add_signals_meta_to_proseg
        
        run_proseg(input_dir, output_dir, reg_name, rust_bin_path=rust_bin_path, **kwargs["kwargs"])
        add_signals_meta_to_proseg(root_dir, output_dir, reg_name, vpt_bin_path=vpt_bin_path, **kwargs["kwargs"])
    elif type == "vpt":
        from .vpt import run_vpt

        config_path = kwargs.get(
            "config_path",
            "/ceph/cephatlas/aklein/vpt/config_files/cellpose_nuclei_Z3.json",
        )
        run_vpt(input_dir, output_dir, reg_name, config_path=config_path, vpt_bin_path=vpt_bin_path, **kwargs["kwargs"])
    elif type == "cellpose":
        from .backends.cellpose import run_cellpose
        from .vpt import seg_to_vpt

        # TODO: change other segmentation methods to take in the direct image dir
        is_3d = run_cellpose(input_dir, output_dir, reg_name, **kwargs["kwargs"])
        input_dir_vpt = Path(input_dir).parents[1]
        seg_to_vpt(input_dir_vpt, output_dir, reg_name, vpt_bin_path=vpt_bin_path, is_3d=is_3d, **kwargs["kwargs"])
    elif type == "mesmer":
        from .backends.mesmer import run_mesmer
        from .vpt import seg_to_vpt

        run_mesmer(input_dir, output_dir, reg_name, **kwargs["kwargs"])
        seg_to_vpt(input_dir, output_dir, reg_name, vpt_bin_path=vpt_bin_path, **kwargs["kwargs"])

    else:
        raise ValueError(f"Unknown segmentation type: {type}")


def segment_experiment(
    type: str,
    exp_name: str,
    input_dir: Path = None,
    output_dir: Path = None,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    zarr_store: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs
):
    """
    Run segmentation for all regions in an experiment.
    Parameters:
    type (str): Type of segmentation to run (e.g., "proseg", "vpt").
    exp_name (str): Name of the experiment.
    input_dir (Path, optional): Directory containing the input data. Defaults to None.
    output_dir (Path, optional): Directory to save the output data. Defaults to None.
    **kwargs: Additional keyword arguments to pass to the segmentation function.
    """

    if zarr_store is not None:
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    exp_path = Path(f"{zarr_store}/{exp_name}")
    regions = glob.glob(f"{exp_path}/region_*")
    for reg in regions:
        reg_name = reg.split("/")[-1]
        run_segmentation(
            type=type,
            exp_name=exp_name,
            reg_name=reg_name,
            input_dir=input_dir,
            output_dir=output_dir,
            root_path=root_path,
            segmentation_store=segmentation_store,
            vpt_bin_path=vpt_bin_path,
            rust_bin_path=rust_bin_path,
            **kwargs,
        )


def vpt_on_segmentation(
    exp_name: str,
    reg_name: str,
    input_dir: Path = None,
    output_dir: Path = None,
    vpt_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    Convert segmentation results to VPT format.
    Parameters:
    exp_name (str): Name of the experiment.
    reg_name (str): Name of the region.
    input_dir (Path, optional): Directory containing the input data. Defaults to None.
    output_dir (Path, optional): Directory to save the output data. Defaults to None.
    **kwargs: Additional keyword arguments to pass to the conversion function.
    """

    from .vpt import seg_to_vpt

    if input_dir is None:
        processed_root_path = os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"
    if output_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        output_dir = f"{seg_out_path}/{exp_name}/vpt"

    seg_to_vpt(input_dir, output_dir, reg_name, vpt_bin_path=vpt_bin_path, **kwargs)


#TODO: work on this
def align_geometries(
    exp_name,
    reg_name,
    prefix1 : str = "cellpose_nuc",
    prefix2 : str = "cellpose_cell",
    output_dir : str | None = None,
    geometry_mode : str = "larger", 
    cell_id : str = "EntityID",
    coordinate_system : str = "global",
    min_intersection_area : float = 0.0,
    out_dir_name : str = "align",
    zarr_store : str | None = None,
    segmentation_store: str | None = None,
    root_path : str | None = None,
    vpt_bin_path : str | None = None
): 
    from .vpt import seg_to_vpt
    from .align import align_segmentations

    if zarr_store is None: 
        zarr_store = os.getenv("ZARR_STORAGE_PATH")
    if root_path is None: 
        root_path = os.getenv("PROCESSED_ROOT_PATH")
    input_dir = f"{root_path}/{exp_name}/out"
    if segmentation_store is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        segmentation_store = Path(f"{seg_out_path}/{exp_name}/{out_dir_name}")
    elif exp_name not in str(segmentation_store):
        segmentation_store = Path(f"{segmentation_store}/{exp_name}/{out_dir_name}")
    # if isinstance(segmentation_store, str):
    #     segmentation_store = Path(segmentation_store)

    zarr_path = Path(zarr_store) / exp_name / reg_name


    ret_dict = align_segmentations(
        zarr_path=zarr_path,
        exp_name=exp_name,
        reg_name=reg_name,
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir = output_dir,
        segmentation_store = segmentation_store,
        geometry_mode = geometry_mode,
        cell_id = cell_id,
        min_intersection_area = min_intersection_area,
        coordinate_system = coordinate_system,
    )

    logger.info(f"Converting aligned geometries to VPT format for region {reg_name} in experiment {exp_name}.")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Segmentation store: {segmentation_store}")
    logger.info(f"Input boundaries file: {Path(ret_dict['output_path']).name}")

    input_boundaries = Path(ret_dict['output_path']).name
    seg_to_vpt(input_dir, segmentation_store, reg_name, vpt_bin_path, input_boundaries=input_boundaries)


# @depreceated in the current version of proseg (v > 3.0)
def align_proseg(
    exp_name: str,
    reg_name: str,
    seed_prefix_name: str = "default",
    prefix_name: str = "proseg",
    out_prefix_name: str = "proseg_aligned",
    input_dir: str | Path = None,
    seg_dir: str | Path = None,
    x: str = "x",
    y: str = "y",
    z: str = "global_z",
    cell_column: str = "cell_id",
    barcode_column: str = "barcode_id",
    gene_column: str = "gene",
    fov_column: str = "fov",
    cell_missing: int = "-1",
    min_jaccard: float = 0.4,
    min_prob: float = 0.5,
    filter_blank: bool = False,
    cell_metadata_fname: str = "merged_cell_metadata.csv",
    cell_by_gene_fname: str = "merged_cell_by_gene.csv",
    detected_transcripts_fname: str = "merged_transcript_metadata.csv",
    cell_polygons_fname: str = "merged_cell_polygons.geojson",
    vpt_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    Align Proseg transcripts to seed transcripts.
    """

    from .backends.proseg import align_proseg_transcripts
    from .vpt import generate_metadata  # , seg_to_vpt
    # from spida.S.io import load_segmentation_region

    if input_dir is None:
        processed_root_path = os.getenv("PROCESSED_ROOT_PATH")
        input_dir = f"{processed_root_path}/{exp_name}/out"

    if seg_dir is None:
        seg_out_path = os.getenv("SEGMENTATION_OUT_PATH")
        seg_dir = Path(f"{seg_out_path}/{exp_name}/{prefix_name}")
    if isinstance(seg_dir, str):
        seg_dir = Path(seg_dir)

    logger.info(
        f"Aligning Proseg transcripts for region {reg_name} in experiment {exp_name}."
    )
    # Aligning the proseg transcripts to the seed transcripts
    align_proseg_transcripts(
        exp_name=exp_name,
        reg_name=reg_name,
        seed_prefix_name=seed_prefix_name,
        prefix_name=prefix_name,
        x=x,
        y=y,
        z=z,
        cell_column=cell_column,
        barcode_column=barcode_column,
        gene_column=gene_column,
        fov_column=fov_column,
        cell_missing=cell_missing,
        min_jaccard=min_jaccard,
        min_prob=min_prob,
        filter_blank=filter_blank,
        merged_transcripts_fname=detected_transcripts_fname,
        merged_cell_by_gene_fname=cell_by_gene_fname,
        merged_cell_polygons_fname=cell_polygons_fname,
        save_dir=seg_dir / reg_name,
    )

    logger.info(
        f"Generating metadata for aligned transcripts in region {reg_name} of experiment {exp_name}."
    )
    # Generating the metadata for the aligned transcripts
    generate_metadata(
        root_dir=input_dir,
        seg_out_dir=seg_dir,
        region=reg_name,
        input_boundaries=cell_polygons_fname,
        output_boundaries="merged_converted_boundaries.parquet",
        input_transcripts=detected_transcripts_fname,
        input_entity_by_gene=cell_by_gene_fname,
        output_metadata=cell_metadata_fname,
        output_signals="merged_sum_signals.csv",
        output_entity_by_gene="cell_by_gene.csv",
        output_transcripts="detected_transcripts.csv",
        vpt_bin_path=vpt_bin_path,
    )

# ===========================================================================
# Stage-2 orchestrators: segment (backend env) + process (preprocessing env).
# These supersede the bundled `run_segmentation` above (kept only for the
# legacy batch/align helpers). The CLI exposes `segment-region` and
# `process-segmentation-region`; `run-segmentation-region` is deprecated.
# ===========================================================================
def segment_region(
    method: str,
    exp_name: str,
    reg_name: str,
    version: str | None = None,
    *,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs,
):
    """Run only the segmentation *backend* for a region (produces raw boundaries).

    Runs in the backend's env (``spec.env``): cellpose/mesmer -> "cellpose",
    proseg -> "preprocessing". Output goes to
    ``{SEGMENTATION_OUT_PATH}/{exp}/{method}/{region}``.
    """
    from .backends import get_spec

    spec = get_spec(method, version)
    spec.require_env("segment")

    processed_root = root_path or os.getenv("PROCESSED_ROOT_PATH")
    seg_out = segmentation_store or os.getenv("SEGMENTATION_OUT_PATH")
    input_dir = f"{processed_root}/{exp_name}/out"
    output_dir = f"{seg_out}/{exp_name}/{method}"
    logger.info("segment_region: %s v%s on %s/%s", method, spec.version, exp_name, reg_name)

    if method == "cellpose":
        # cellpose/mesmer read mosaic stain tiffs from the region's images dir;
        # proseg reads {out}/{region}/detected_transcripts.csv (so it takes {out}).
        from .backends.cellpose import run_cellpose
        run_cellpose(f"{input_dir}/{reg_name}/images", output_dir, reg_name, **kwargs)
    elif method == "mesmer":
        from .backends.mesmer import run_mesmer
        run_mesmer(f"{input_dir}/{reg_name}/images", output_dir, reg_name, **kwargs)
    elif method == "proseg":
        from .backends.proseg import run_proseg
        run_proseg(input_dir, output_dir, reg_name, rust_bin_path=rust_bin_path, **kwargs)
    else:
        raise ValueError(f"segment_region: unknown method {method!r}")


def process_segmentation_region(
    method: str,
    exp_name: str,
    reg_name: str,
    version: str | None = None,
    *,
    backend: str = "native",
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    micron_per_z: float = 1.5,
    n_z_planes: int = 7,
    n_jobs: int = 7,
    vpt_bin_path: str | Path | None = None,
    **kwargs,
):
    """Normalize + post-process a segmenter's raw output into the segmentation schema.

    Runs in the ``preprocessing`` env. ``backend="native"`` (default) uses the
    pure-Python path (ingest + the steps in ``spec.needs``); ``backend="vpt"``
    falls back to the VPT binary path (``seg_to_vpt``). Writes the standardized
    segmentation schema files to ``{SEGMENTATION_OUT_PATH}/{exp}/{method}/{region}``.
    """
    import pandas as pd
    from .backends import get_spec
    from .backends.base import (
        SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_CELL_METADATA,
        SCHEMA_TRANSCRIPTS, SCHEMA_SUM_SIGNALS,
    )

    spec = get_spec(method, version)
    spec.require_env("process")

    processed_root = root_path or os.getenv("PROCESSED_ROOT_PATH")
    seg_out = segmentation_store or os.getenv("SEGMENTATION_OUT_PATH")
    raw_region = Path(processed_root) / exp_name / "out" / reg_name
    images_dir = raw_region / "images"
    m2m = images_dir / "micron_to_mosaic_pixel_transform.csv"
    seg_region = Path(seg_out) / exp_name / method / reg_name

    logger.info("process_segmentation_region: %s v%s (%s backend) needs=%s",
                method, spec.version, backend, spec.needs)

    if backend == "vpt":
        from .vpt import seg_to_vpt
        seg_to_vpt(str(raw_region.parents[0]), str(seg_region.parents[0]),
                   reg_name, vpt_bin_path=vpt_bin_path, **kwargs)
        return

    from .ingest import ingest_polygons
    from .segmentation_utils import (
        partition_transcripts, derive_entity_metadata, sum_signals,
    )

    pixel = method in ("cellpose", "mesmer")
    cir_boundaries = seg_region / SCHEMA_BOUNDARIES

    # 1. raw boundaries -> segmentation schema boundaries (micron space)
    ingest_polygons(
        spec, seg_region / spec.boundaries_file, cir_boundaries,
        micron_to_mosaic=str(m2m) if pixel else None,
        micron_per_z=micron_per_z, n_z_planes=n_z_planes,
    )

    # 2. cell-by-gene (only if the method doesn't provide counts itself)
    if "partition_transcripts" in spec.needs:
        partition_transcripts(
            cir_boundaries, raw_region / spec.raw_transcripts_file,
            output_entity_by_gene=seg_region / SCHEMA_CELL_BY_GENE,
            output_transcripts=seg_region / SCHEMA_TRANSCRIPTS,
        )

    # 3. per-cell geometric metadata (derived) or native (proseg)
    if "derive_entity_metadata" in spec.needs:
        meta = derive_entity_metadata(cir_boundaries, seg_region / SCHEMA_CELL_BY_GENE)
    else:
        native = pd.read_csv(seg_region / spec.metadata_file)
        meta = native.set_index(spec.columns.cell_id)
        meta.index.name = "EntityID"      # proseg cell id == segmentation schema EntityID

    # 4. image intensity (always) and merge into cell_metadata
    signals = sum_signals(cir_boundaries, images_dir, m2m,
                          output_csv=seg_region / SCHEMA_SUM_SIGNALS, n_jobs=n_jobs)
    merged = meta.join(signals)
    merged.to_csv(seg_region / SCHEMA_CELL_METADATA)
    logger.info("process_segmentation_region: wrote segmentation schema outputs to %s", seg_region)


def _deprecated_run_segmentation_message(method: str = "cellpose",
                                         exp_name: str = "EXP",
                                         reg_name: str = "REGION") -> str:
    """Migration guidance for the retired bundled `run-segmentation-region`."""
    from .backends import get_spec
    try:
        env = get_spec(method).env
    except Exception:
        env = "cellpose"
    return (
        "`run-segmentation-region` is deprecated and no longer runs.\n"
        "Segmentation and its post-processing are now two steps in two envs "
        "(and output filenames changed, e.g. `boundaries_micron.parquet`):\n\n"
        f"  pixi run -e {env} python -m spida.S.cli segment-region "
        f"{method} {exp_name} {reg_name}\n"
        f"  pixi run -e preprocessing python -m spida.S.cli process-segmentation-region "
        f"{method} {exp_name} {reg_name}\n\n"
        "Use `--backend vpt` on process-segmentation-region for the legacy VPT path."
    )
