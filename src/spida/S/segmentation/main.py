from pathlib import Path
import logging

logger = logging.getLogger(__package__)


def align_geometries(
    dataset,
    prefix1: str = "cellpose_nuc",
    prefix2: str = "cellpose_cell",
    geometry_mode: str = "larger",
    cell_id: str = "EntityID",
    coordinate_system: str = "global",
    min_intersection_area: float = 0.0,
    out_dir_name: str = "align",
    micron_per_z: float = 1.5,
    n_jobs: int = 7,
):
    """Reconcile two segmentations (e.g. nuclear ``cellpose_nuc`` vs whole-cell
    ``cellpose_cell``) into one boundary set by spatial overlap, then build the
    segmentation schema from it natively (``process_custom_segmentation``, no VPT).

    Output: ``{seg}/{out_dir_name}`` — the aligned polygons parquet plus the schema files.
    2D: each aligned polygon is one plane, expanded across the transcript planes for
    partitioning, as for any 2D custom segmentation.

    Parameters:
    dataset (Dataset): The dataset (zarr + raw + seg locations).
    prefix1 (str): First segmentation (default "cellpose_nuc").
    prefix2 (str): Second segmentation (default "cellpose_cell").
    geometry_mode (str): Which polygon to keep per matched pair: "larger", "prefix1", "prefix2", "intersection" (default "larger").
    cell_id (str): The cell-id column name (default "EntityID").
    coordinate_system (str): Coordinate system to align in (default "global" = mosaic pixels for MERSCOPE).
    min_intersection_area (float): Minimum intersection area (pixels) required for a match (default 0.0).
    out_dir_name (str): Label of the aligned segmentation under the dataset's seg dir (default "align").
    micron_per_z (float): Micron thickness per z-plane (default 1.5).
    n_jobs (int): sum-signals workers (default 7).
    """
    from spida._constants import SHAPES_KEY
    from .align import align_segmentations

    raw = dataset.require("raw")
    out_dir = dataset.seg_dir(out_dir_name)
    ret = align_segmentations(
        zarr_path=dataset.require("zarr"),
        shape_key1=dataset.keys(prefix1)[SHAPES_KEY],
        shape_key2=dataset.keys(prefix2)[SHAPES_KEY],
        prefix1=prefix1,
        prefix2=prefix2,
        output_dir=out_dir,
        geometry_mode=geometry_mode,
        cell_id=cell_id,
        min_intersection_area=min_intersection_area,
        coordinate_system=coordinate_system,
    )

    transcripts = raw / "detected_transcripts.parquet"
    if not transcripts.exists():
        transcripts = raw / "detected_transcripts.csv"
    logger.info("align_geometries: %d aligned cells -> schema in %s", ret["n_aligned"], out_dir)
    process_custom_segmentation(
        ret["output_path"],
        out_dir,
        boundaries_space="pixel" if coordinate_system == "global" else "micron",
        micron_to_mosaic_path=raw / "images" / "micron_to_mosaic_pixel_transform.csv",
        z_spacing=micron_per_z,
        n_z_planes=1,
        cell_id_col=f"{prefix1}_id",
        transcripts_path=transcripts,
        images_dir=raw / "images",
        n_jobs=n_jobs,
    )


# ===========================================================================
# Stage-2 orchestrators: segment (backend env) + process (preprocessing env).
# These replace the old bundled single-env controller. The CLI exposes
# `segment-region` and `process-segmentation-region`; the legacy
# `run-segmentation-region` remains only as a deprecation redirect.
# ===========================================================================

#: Filename that ``P/tissue_mask.py`` writes into a region's images dir (its``QC_MASK_TIF_NAME``).
TISSUE_MASK_FILENAME = "tissue_mask.tif"

#: Methods whose backend accepts ``tissue_mask_path``. this set
#: only records which backends have been wired to consume it so far.
_TISSUE_MASK_METHODS = frozenset({"cellpose"})


def resolve_tissue_mask_path(
    explicit: str | Path | None, images_dir: str | Path
) -> Path | None:
    """Locate a region's tissue mask, or ``None`` when it has none.

    Explicit and auto-discovered paths are treated differently on purpose: a
    mask the caller *asked* for must exist, whereas absence of the
    conventional file simply means this region has not had ``tissue-mask`` run
    and should segment ungated.
    """
    if explicit is not None:
        path = Path(explicit).expanduser()
        if path.is_dir():
            path = path / TISSUE_MASK_FILENAME
        if not path.is_file():
            raise FileNotFoundError(
                f"tissue_mask_path={explicit} does not exist (looked for "
                f"{path}). Run `tissue-mask` for this region first, or omit "
                f"--tissue_mask_path to segment without a tissue gate."
            )
        logger.info("tissue mask (explicit): %s", path)
        return path

    path = Path(images_dir) / TISSUE_MASK_FILENAME
    if path.is_file():
        logger.info("tissue mask (auto-discovered): %s", path)
        return path
    logger.info("no tissue mask at %s; segmenting the full extent", path)
    return None


def segment_region(
    method: str,
    dataset,
    *,
    prefix_name: str | None = None,
    rust_bin_path: str | Path | None = None,
    tissue_mask_path: str | Path | None = None,
    **kwargs,
):
    """Run only the segmentation *backend* for a region (produces raw boundaries).

    Runs in the backend's env (``spec.env``): cellpose -> "cellpose", mesmer ->
    "deepcell", proseg -> "preprocessing". Output goes to ``dataset.seg_dir(prefix_name)``
    (``prefix_name`` defaults to the method; set it to e.g. ``cellpose_nuc``/``cellpose_cell``).

    Parameters:
    dataset (Dataset): The dataset (raw region dir, segmentation dir).
    prefix_name (str | None): Segmentation name / output-directory leaf; defaults to the
        method name. This is the same label used as the ``prefix_name`` at load time
        (default is None).
    tissue_mask_path (str | Path | None): Boolean tissue mask restricting where
        segmentation runs, as written by ``tissue-mask``. Default None auto-discovers
        ``<raw images dir>/tissue_mask.tif`` and proceeds ungated if it is absent;
        an explicit path that does not exist raises.
    """
    from .backends import get_spec

    spec = get_spec(method)
    spec.require_env("segment")

    label = prefix_name or method
    raw_region = dataset.require("raw")
    reg_name = raw_region.name
    input_dir = str(raw_region.parent)
    images_dir = str(dataset.raw_images())
    output_dir = str(dataset.seg_dir(label))
    logger.info("segment_region: %s v%s on %s -> %s",
                method, spec.version, raw_region, output_dir)

    if tissue_mask_path is not None and method not in _TISSUE_MASK_METHODS:
        raise ValueError(
            f"tissue_mask_path was given but the {method!r} backend does not "
            f"consume it yet (wired for: {', '.join(sorted(_TISSUE_MASK_METHODS))})."
        )

    if method == "cellpose":
        # cellpose/mesmer read mosaic stain tiffs from the region's images dir;
        # proseg reads {out}/{region}/detected_transcripts.csv (so it takes {out}).
        from .backends.cellpose import run_cellpose
        run_cellpose(images_dir, output_dir, reg_name,
                     tissue_mask_path=resolve_tissue_mask_path(
                         tissue_mask_path, images_dir),
                     **kwargs)
    elif method == "mesmer":
        from .backends.mesmer import run_mesmer
        run_mesmer(images_dir, output_dir, reg_name, **kwargs)
    elif method == "proseg":
        from .backends.proseg import run_proseg
        run_proseg(input_dir, output_dir, reg_name, rust_bin_path=rust_bin_path, **kwargs)
    else:
        raise ValueError(f"segment_region: unknown method {method!r}")


def process_segmentation_region(
    method: str,
    dataset,
    version: str | None = None,
    *,
    prefix_name: str | None = None,
    backend: str = "native",
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
    segmentation schema files to ``dataset.seg_dir(prefix_name)``.

    Parameters:
    dataset (Dataset): The dataset (raw region dir, segmentation dir).
    prefix_name (str | None): Segmentation name / directory leaf; must match what was
        passed to ``segment-region``. Defaults to the method name (default is None).
    """
    import pandas as pd
    from .backends import get_spec
    from .backends.base import (
        SCHEMA_BOUNDARIES, SCHEMA_CELL_BY_GENE, SCHEMA_CELL_METADATA,
        SCHEMA_TRANSCRIPTS, SCHEMA_SUM_SIGNALS,
    )

    spec = get_spec(method, version)
    spec.require_env("process")

    label = prefix_name or method
    raw_region = dataset.require("raw")
    images_dir = dataset.raw_images()
    m2m = images_dir / "micron_to_mosaic_pixel_transform.csv"
    # where segment-region wrote the raw output (current layout, legacy fallback)
    seg_region = dataset.seg_dir(label, must_exist=True)

    logger.info("process_segmentation_region: %s v%s (%s backend) needs=%s -> %s",
                method, spec.version, backend, spec.needs, seg_region)

    if backend == "vpt":
        from .vpt import seg_to_vpt
        seg_to_vpt(str(raw_region), str(seg_region), vpt_bin_path=vpt_bin_path,
                   is_3d=(n_z_planes > 1), spacing_z=micron_per_z, **kwargs)
        return

    from .ingest import ingest_polygons
    from .postprocess import (
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
    from .postprocess import (
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


