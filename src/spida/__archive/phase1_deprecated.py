"""Archived 2026-09-28 (IO redesign phase 1, Amit): everything that carried a deprecation warning. Not importable as-is."""

# ==== S/segmentation/backends/proseg.py :: _execute_cli_proseg_v2 ====
@DeprecationWarning
def _execute_cli_proseg_v2(root_dir: str, output_dir: str, region: str, **proseg_params):
    """
    Run the proseg command with the prespecified parameters using subprocess.
    """

    # The proseg command to run
    proseg_run_cmd = f"""
        proseg --merscope \
        {root_dir}/{region}/detected_transcripts.csv \n
        --output-path {output_dir} \n
        --output-expected-counts expected-counts.csv.gz \n
        --output-cell-metadata cell-metadata.csv.gz \n
        --output-transcript-metadata transcript-metadata.csv.gz \n
        --output-cell-polygons cell-polygons.geojson.gz \n
        --output-cell-polygon-layers cell-polygons-layers.geojson.gz \n
        --output-union-cell-polygons union-cell-polygons.geojson.gz \n
        """
    for _key, _val in proseg_params.items():
        if _key not in proseg_args:
            logger.info("skipping unknown proseg parameter: %s" % _key)
            continue
        if isinstance(_val, bool):
            if _val:
                logger.info(f"Adding boolean proseg parameter: --{_key}")
                proseg_run_cmd += f"--{_key} \n"
        else:
            logger.info(f"Adding proseg parameter: --{_key} {_val}")
            proseg_run_cmd += f"--{_key} {_val} \n"

    # proseg_run_cmd += f"--nthreads {os.cpu_count()} \\\n"
    logger.info(f"Running proseg with command:\n{proseg_run_cmd}")

    # run command and report results
    ret = subprocess.run(proseg_run_cmd.split(), capture_output=True, check=True)
    return ret

# ==== S/segmentation/backends/proseg.py :: filt_to_ids ====
@DeprecationWarning
def filt_to_ids(meta_path, geom_path, tz_path, cbg_path):
    """
    Filter metadata and geometry to only include cells with IDs in the provided list.
    Additionally rename the index of the merged_cell_by_gene
    """
    import pandas as pd
    import geopandas as gpd

    cbg = pd.read_csv(cbg_path, index_col=0)
    transcripts = pd.read_csv(tz_path, index_col=0)
    meta = pd.read_csv(meta_path)
    geom = gpd.read_parquet(geom_path)
    ids = transcripts["cell_id"].unique()
    meta_filt = meta[meta["EntityID"].isin(ids)]
    geom_filt = geom[geom["EntityID"].isin(ids)]
    cbg.index = cbg.index.astype(str)
    cbg.index.name = "cell"

    meta_filt.to_csv(meta_path, index=False)
    geom_filt.to_parquet(geom_path, index=False)
    cbg.to_csv(cbg_path)
    logger.info(
        f"Filtered metadata and geometry to {len(meta_filt)} cells with IDs in {len(ids)} unique IDs."
    )
    logger.info(
        f"Metadata shape: {meta_filt.shape}, Geometry shape: {geom_filt.shape}, CBG shape: {cbg.shape}"
    )

# ==== S/cli.py :: filt_to_ids ====
@cli.command(
    name="filt-to-ids",
    aliases=["filt_to_ids", "filter_to_ids", "filter-to-ids"],
    cls=RichCommand,
    help="['filt_to_ids', 'filter_to_ids', 'filter-to-ids'] - Filter metadata and geometry to only include cells with IDs in the provided list.",
)
@click.argument("meta_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("geom_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("tz_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("cbg_path", type=click.Path(exists=True, dir_okay=False))
@click.pass_context
def filt_to_ids(
    ctx,
    meta_path: str | Path,
    geom_path: str | Path,
    tz_path: str | Path,
    cbg_path: str | Path,
):
    """
    Filter metadata and geometry to only include cells with IDs in the provided list.
    Additionally rename the index of the merged_cell_by_gene
    """
    from .segmentation.backends.proseg import filt_to_ids as func
    func(
        meta_path=meta_path,
        geom_path=geom_path,
        tz_path=tz_path,
        cbg_path=cbg_path,
    )

# ==== S/cli.py :: run_segmentation ====
@cli.command(
    name="run-segmentation-region",
    aliases=["run", "run_segmentation_region"],
    cls=RichCommand,
    help="['run', 'run_segmentation_region'], Run segmentation of kind TYPE for a given EXP_NAME and REG_NAME.",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.argument("seg_type", type=str)
@click.argument("exp_name", type=str)
@click.argument("reg_name", type=str)
@click.option("input_dir", "--input_dir", default=None, type=click.Path(), help="Directory containing the input data (default: None, uses environment variable)")
@click.option("output_dir", "--output_dir", default=None, type=click.Path(), help="Directory to save the output data (default: None, uses environment variable)")
@click.option("root_path", "--root_path", default=None, type=click.Path(), help="Root path for the data (default: None, uses environment variable)")
@click.option("segmentation_store", "--segmentation_store", default=None, type=click.Path(), help="Path to store the segmentation results (default: None, uses environment variable)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=None, type=click.Path(), help="Path to the VPT binary (default: None, uses environment variable)")
@click.option("rust_bin_path", "--rust_bin_path", default=None, type=click.Path(), help="Path to the Rust binary (default: None, uses environment variable)")
@click.pass_context
def run_segmentation(
    ctx,
    seg_type: str,
    exp_name: str,
    reg_name: str,
    input_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    root_path: str | Path | None = None,
    segmentation_store: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    rust_bin_path: str | Path | None = None,
    **kwargs,
):
    """
    DEPRECATED. Segmentation and its post-processing are now separate commands
    (see the printed migration guidance).
    """
    from .segmentation.main import _deprecated_run_segmentation_message

    raise click.ClickException(
        _deprecated_run_segmentation_message(seg_type, exp_name, reg_name)
    )

# ==== S/segmentation/main.py :: _deprecated_run_segmentation_message ====
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

# ==== utilities/degs.py :: call_degs_by_celltype ====
def call_degs_by_celltype(*args, **kwargs) -> dict[str, pd.DataFrame]:
    """Deprecated. Use call_degs_scanpy instead."""
    warnings.warn(
        "call_degs_by_celltype is deprecated; use call_degs_scanpy instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return call_degs_scanpy(*args, **kwargs)

# ==== S/segmentation/backends/cellpose_config.py :: CellposeConfig.from_kwargs ====
    @classmethod
    def from_kwargs(cls, **kwargs) -> "CellposeConfig":
        """As the base, but first folds the legacy mode flags into ``seg_mode``.

        ``do_3D`` and ``project_3d_to_2d`` are not simple renames — they select
        a mode — so they are translated here rather than via ``ALIASES``.
        """
        legacy = {}
        for old in ("do_3D", "project_3d_to_2d"):
            if old in kwargs:
                legacy[old] = kwargs.pop(old)

        if legacy and "seg_mode" in kwargs:
            raise ConfigError(
                f"seg_mode cannot be combined with the deprecated flag(s) "
                f"{', '.join(sorted(legacy))}; use seg_mode alone."
            )
        if legacy:
            if legacy.get("project_3d_to_2d"):
                kwargs["seg_mode"] = "2d"
                kwargs.setdefault("z_reduce", "mip")
            elif legacy.get("do_3D"):
                kwargs["seg_mode"] = "3d_true"
            else:
                kwargs["seg_mode"] = "3d_stitched"
            names = ", ".join(sorted(legacy))
            # Same signalling as ALIASES so every deprecated name behaves alike.
            warnings.warn(
                f"CellposeConfig: {names} is deprecated; use "
                f"seg_mode={kwargs['seg_mode']!r}.",
                DeprecationWarning,
                stacklevel=2,
            )
            logger.warning(
                "CellposeConfig: %s is deprecated; interpreted as seg_mode=%r",
                names, kwargs["seg_mode"],
            )

        return super().from_kwargs(**kwargs)  # type: ignore[return-value]

# ==== S/segmentation/backends/cellpose_config.py :: CellposeConfig.ALIASES ====
    ALIASES: ClassVar[dict[str, str]] = {
        "model_name": "model",
        "tile_size": "poly_tile_size",
        "overlap": "poly_tile_overlap",
        "clahe_clipLimit": "clahe_clip_limit",
        "clahe_tileGridSize": "clahe_tile_grid_size",
    }
