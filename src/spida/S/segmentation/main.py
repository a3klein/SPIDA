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
    *,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs,
):
    """Run only the segmentation *backend* for a region (produces raw boundaries).

    Runs in the backend's env (``spec.env``): cellpose -> "cellpose", mesmer ->
    "deepcell", proseg -> "preprocessing". Output goes to
    ``{SEGMENTATION_OUT_PATH}/{exp}/{method}/{region}``.
    """
    from .backends import get_spec

    spec = get_spec(method)
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


def process_custom_segmentation(
    boundaries_path: str | Path,
    output_dir: str | Path,
    *,
    boundaries_space: str = "micron",
    micron_to_mosaic_path: str | Path | None = None,
    z_spacing: float = 1.5,
    n_z_planes: int = 1,
    cell_id_col: str | None = None,
    boundary_z_col: str | None = None,
    entity_type: str = "cell",
    transcripts_path: str | Path | None = None,
    transcript_z_col: str = "global_z",
    transcript_z_in_microns: bool = False,
    gene_col: str = "gene",
    barcode_col: str = "barcode_id",
    transcript_x_col: str = "global_x",
    transcript_y_col: str = "global_y",
    images_dir: str | Path | None = None,
    segmentation_z_index: int | None = None,
    metadata_path: str | Path | None = None,
    metadata_cell_id_col: str | None = None,
    n_jobs: int = 7,
) -> None:
    """Process a *user-provided* segmentation into the segmentation schema (native, no VPT).

    The spec-free sibling of :func:`process_segmentation_region`: instead of a registry
    backend, the caller describes the layout of their own files and this composes the same
    native building blocks (``ingest_custom_polygons`` -> ``partition_transcripts`` ->
    ``derive_entity_metadata`` -> ``sum_signals``) into the standardized schema files
    (``boundaries_micron.parquet``, ``cell_by_gene.csv``, ``detected_transcripts.csv``,
    ``cell_metadata.csv``, ``sum_signals.csv``) that ``load_segmentation`` reads directly.
    Steps are auto-selected by what you provide (transcripts -> partition; images ->
    sum-signals; metadata -> merged onto the derived metadata).

    Z-indexing (two independent conventions, matching VPT): **transcripts** are placed on a
    z-plane via ``transcript_z_col`` (integer plane, or micron ``ZLevel`` when
    ``transcript_z_in_microns``, converted with ``z_spacing``); **images** are indexed by the
    integer ``ZIndex`` in their ``mosaic_{stain}_z{N}.tif`` filenames. All z-planes share the
    single ``z_spacing``. When ``n_z_planes == 1`` the input is 2D: each polygon is expanded
    into a cylinder across the transcript (else image) planes for partitioning, while
    sum-signals is computed on a *single* plane only (``segmentation_z_index``, default the
    middle image plane) — the plane the segmentation was run on, not the whole stack. When
    ``n_z_planes > 1`` the boundary and transcript plane counts must both equal ``n_z_planes``
    (a mismatch raises); an image plane-count mismatch only warns.

    Parameters:
    boundaries_path (str | Path): Path to the user's polygons file (``.parquet`` or ``.geojson``/``.geojson.gz``).
    output_dir (str | Path): Directory to write the segmentation-schema output files into.
    boundaries_space (str): Coordinate space of the input polygons: ``"micron"`` or ``"pixel"`` (default "micron").
    micron_to_mosaic_path (str | Path | None): Path to the ``micron_to_mosaic_pixel_transform.csv`` file;
        required when ``boundaries_space="pixel"`` and/or when ``images_dir`` is given (sum-signals maps micron->pixel). Default is None.
    z_spacing (float): Micron thickness per z-plane; must be consistent between boundaries and transcripts (default 1.5).
    n_z_planes (int): Number of z-planes of the segmentation; 1 => 2D (cylinder expansion) (default 1).
    cell_id_col (str | None): Column naming the caller's cell identifier in the boundaries; groups rows into cells and is 
        preserved in the output. Required for 3D input and for merging ``metadata_path`` (default is None).
    boundary_z_col (str | None): Integer z-plane column in the boundaries for 3D input; leave None for 2D (default is None).
    entity_type (str): Value for the schema ``Type`` column (default "cell").
    transcripts_path (str | Path | None): Detected transcripts (``.csv`` or ``.parquet``); if given, drives ``partition_transcripts`` (default is None).
    transcript_z_col (str): Column giving each transcript's z (default "global_z").
    transcript_z_in_microns (bool): If True, ``transcript_z_col`` is micron ``ZLevel`` and is converted to an integer plane via ``z_spacing`` (default False).
    gene_col (str): Transcript gene-name column (default "gene").
    barcode_col (str): Transcript barcode-id column (default "barcode_id").
    transcript_x_col (str): Transcript x-coordinate column (default "global_x").
    transcript_y_col (str): Transcript y-coordinate column (default "global_y").
    images_dir (str | Path | None): Directory of ``mosaic_{stain}_z{N}.tif`` stain images; if given, drives ``sum_signals`` (default is None).
    segmentation_z_index (int | None): For 2D input, the single ``ZIndex`` plane the segmentation was run on, used for sum-signals; None => the middle image plane (default is None).
    metadata_path (str | Path | None): Optional user cell-metadata CSV; its columns are merged onto the derived metadata on the cell-id field. Requires ``cell_id_col`` (default is None).
    metadata_cell_id_col (str | None): Cell-id column in ``metadata_path`` if it differs from ``cell_id_col`` (default is None => use ``cell_id_col``).
    n_jobs (int): Parallel workers for sum-signals (default 7).
    """
    import re
    import glob
    import numpy as np
    import pandas as pd
    from .backends.base import (
        SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_CELL_METADATA,
        SCHEMA_TRANSCRIPTS, SCHEMA_SUM_SIGNALS,
    )
    from .ingest import ingest_custom_polygons
    from .segmentation_utils import (
        partition_transcripts, derive_entity_metadata, sum_signals,
        _ENTITY_COL, _ZINDEX_COL,
    )

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    boundaries_out = out / SCHEMA_BOUNDARIES
    if images_dir is not None and micron_to_mosaic_path is None:
        raise ValueError("images_dir given but micron_to_mosaic_path is None; sum-signals needs "
                         "the micron->mosaic transform to map polygons into pixel space.")
    logger.info("process_custom_segmentation: %s -> %s (n_z_planes=%d, space=%s)",
                boundaries_path, out, n_z_planes, boundaries_space)

    def _read_z(path, col):
        s = (pd.read_csv(path, usecols=[col])[col] if str(path).endswith(".csv")
             else pd.read_parquet(path, columns=[col])[col])
        if transcript_z_in_microns:
            return (np.round(s / z_spacing).astype("int64") - 1)
        return s.astype("int64")

    # discover data z-planes (transcripts by transcript_z_col, images by filename ZIndex)
    tx_planes = None
    if transcripts_path is not None:
        tx_planes = sorted(pd.unique(_read_z(transcripts_path, transcript_z_col)).tolist())
        logger.info("process_custom_segmentation: transcripts span z-planes %s", tx_planes)
    img_planes = None
    if images_dir is not None:
        rex = re.compile(r"mosaic_[\w-]+_z(\d+)\.tif$")
        img_planes = sorted({int(m.group(1)) for f in glob.glob(f"{images_dir}/mosaic_*_z*.tif")
                             if (m := rex.search(f))})
        logger.info("process_custom_segmentation: images span z-planes %s", img_planes)

    # resolve target planes + validate plane matching
    if n_z_planes == 1:
        if boundary_z_col is not None:
            raise ValueError("n_z_planes=1 (2D) but boundary_z_col is set; drop boundary_z_col for "
                             "2D input, or set n_z_planes to the boundary plane count for 3D.")
        target_planes = tx_planes or img_planes or [0]
        logger.info("process_custom_segmentation: 2D input -> cylinder across planes %s", target_planes)
    else:
        if boundary_z_col is None:
            raise ValueError(f"n_z_planes={n_z_planes} (3D) requires boundary_z_col naming the "
                             "z-plane column in the boundaries.")
        if tx_planes is not None and len(tx_planes) != n_z_planes:
            raise ValueError(f"transcripts span {len(tx_planes)} z-planes but n_z_planes={n_z_planes}; "
                             "boundary and transcript planes must match.")
        if img_planes is not None and len(img_planes) != n_z_planes:
            logger.warning("process_custom_segmentation: images span %d z-planes but n_z_planes=%d; "
                           "sum-signals will use only overlapping planes.", len(img_planes), n_z_planes)
        target_planes = [0]

    # boundaries -> segmentation schema
    boundaries = ingest_custom_polygons(
        boundaries_path, boundaries_out,
        target_planes=target_planes,
        boundaries_space=boundaries_space, micron_to_mosaic=micron_to_mosaic_path,
        z_spacing=z_spacing, cell_id_col=cell_id_col,
        boundary_z_col=(boundary_z_col if n_z_planes > 1 else None),
        entity_type=entity_type,
    )
    if n_z_planes > 1:
        got = int(boundaries[_ZINDEX_COL].nunique())
        if got != n_z_planes:
            raise ValueError(f"boundaries span {got} z-planes but n_z_planes={n_z_planes}; "
                             "boundary and transcript planes must match.")

    # partition transcripts (optional)
    tmp_files = []
    if transcripts_path is not None:
        tpath, z_arg = transcripts_path, transcript_z_col
        if transcript_z_in_microns:
            tpath = out / "_tx_zindexed.parquet"
            tdf = (pd.read_csv(transcripts_path) if str(transcripts_path).endswith(".csv")
                   else pd.read_parquet(transcripts_path))
            tdf["_ZIndex"] = np.round(tdf[transcript_z_col] / z_spacing).astype("int64") - 1
            tdf.to_parquet(tpath)
            z_arg = "_ZIndex"
            tmp_files.append(tpath)
            logger.info("process_custom_segmentation: converted transcript micron z -> integer plane")
        partition_transcripts(
            boundaries_out, tpath,
            output_entity_by_gene=out / SCHEMA_CELL_BY_GENE,
            output_transcripts=out / SCHEMA_TRANSCRIPTS,
            gene_col=gene_col, barcode_col=barcode_col,
            x_col=transcript_x_col, y_col=transcript_y_col, z_col=z_arg,
        )

    # cell metadata; merge user metadata if provided
    cbg_path = out / SCHEMA_CELL_BY_GENE
    meta = derive_entity_metadata(boundaries_out, cbg_path if cbg_path.exists() else None)
    if metadata_path is not None:
        if cell_id_col is None:
            raise ValueError("metadata_path requires cell_id_col so user metadata can be merged "
                             "onto the assigned EntityIDs.")
        mcid = metadata_cell_id_col or cell_id_col
        # EntityID <-> cell-id mapping is carried in the boundaries file
        xmap = (boundaries[[_ENTITY_COL, cell_id_col]].drop_duplicates()
                .rename(columns={cell_id_col: "__cellkey__"}))
        user_meta = pd.read_csv(metadata_path)
        user_meta = (user_meta.merge(xmap, left_on=mcid, right_on="__cellkey__", how="inner")
                     .set_index(_ENTITY_COL).drop(columns="__cellkey__"))
        meta = meta.join(user_meta, how="left", rsuffix="_user")
        logger.info("process_custom_segmentation: merged %d user-metadata column(s) on %s",
                    user_meta.shape[1], mcid)

    # sum-signals (optional) + merge into cell metadata
    if images_dir is not None:
        if n_z_planes == 1 and len(target_planes) > 1:
            pool = img_planes or target_planes
            seg_z = (segmentation_z_index if segmentation_z_index is not None
                     else int(pool[len(pool) // 2]))
            sig_bounds = out / "_signals_plane.parquet"
            one = boundaries[boundaries[_ZINDEX_COL] == boundaries[_ZINDEX_COL].min()].copy()
            one[_ZINDEX_COL] = seg_z
            one["ZLevel"] = (seg_z + 1) * z_spacing
            one.to_parquet(sig_bounds)
            tmp_files.append(sig_bounds)
            logger.info("process_custom_segmentation: 2D sum-signals on single plane z=%d", seg_z)
        else:
            sig_bounds = boundaries_out
        signals = sum_signals(sig_bounds, images_dir, micron_to_mosaic_path,
                              output_csv=out / SCHEMA_SUM_SIGNALS, n_jobs=n_jobs)
        meta = meta.join(signals)

    meta.to_csv(out / SCHEMA_CELL_METADATA)
    for f in tmp_files:
        Path(f).unlink(missing_ok=True)
    logger.info("process_custom_segmentation: wrote segmentation schema outputs to %s", out)


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
