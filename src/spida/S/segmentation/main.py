import os
from dotenv import load_dotenv  # type: ignore
from pathlib import Path
import logging

load_dotenv()
logger = logging.getLogger(__package__)


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
    """Reconcile two segmentations (e.g. nuclear ``cellpose_nuc`` vs whole-cell
    ``cellpose_cell``) into one boundary set by spatial overlap, then convert to
    the segmentation schema via VPT.

    .. warning::
        NOT YET MIGRATED to the segmentation-schema redesign and may be buggy.
        Specifically: the ``seg_to_vpt`` call below does not thread
        ``is_3d``/``spacing_z``, so it is effectively **2D-only** (3D aligned
        geometries would get default z-spacing on the 2D path), and the outputs
        still use the pre-redesign ``align``/nuc-cell conventions. Kept as a
        feature but treat results with caution until migrated.

    Parameters:
    exp_name (str): The name of the experiment.
    reg_name (str): The name of the region.
    prefix1 (str): Shapes-layer prefix for the first segmentation (default "cellpose_nuc").
    prefix2 (str): Shapes-layer prefix for the second segmentation (default "cellpose_cell").
    output_dir (str | None): Directory for the aligned boundary output (default is None).
    geometry_mode (str): Which polygon to keep per matched pair: one of "larger", "prefix1", "prefix2", "intersection" (default "larger").
    cell_id (str): The cell-id column name (default "EntityID").
    coordinate_system (str): The spatialdata coordinate system to align in (default "global").
    min_intersection_area (float): Minimum intersection area (pixels) required for a match (default 0.0).
    out_dir_name (str): Sub-directory name under the segmentation store for outputs (default "align").
    zarr_store (str | None): Path to the zarr storage (default is None, uses env ZARR_STORAGE_PATH).
    segmentation_store (str | None): Segmentation output root (default is None, uses env SEGMENTATION_OUT_PATH).
    root_path (str | None): Raw MERSCOPE root (default is None, uses env PROCESSED_ROOT_PATH).
    vpt_bin_path (str | None): Path to the VPT binary (default is None, uses env VPT_BIN_PATH).
    """
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


# ===========================================================================
# Stage-2 orchestrators: segment (backend env) + process (preprocessing env).
# These replace the old bundled single-env controller. The CLI exposes
# `segment-region` and `process-segmentation-region`; the legacy
# `run-segmentation-region` remains only as a deprecation redirect.
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
                   reg_name, vpt_bin_path=vpt_bin_path,
                   is_3d=(n_z_planes > 1), spacing_z=micron_per_z, **kwargs)
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
