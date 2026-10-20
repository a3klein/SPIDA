"""Archived in phase 2 of the IO redesign (2026-10-17): replaced by ``sync-annotations``
(``spida.I._io.sync_annotations``), which writes only obs and runs without spatialdata.
These rewrote the whole zarr table from the AnnData file. Not imported anywhere.
"""
# flake8: noqa
# Originally in spida/I/cli.py.

@cli.command(name="backup-adata", cls=RichCommand)
@click.argument('prefix_name', type=click.STRING)
@click.option('--suffix', type=click.STRING, default='', help='Table suffix. Defaults to none.')
@click.option('--adata_file', type=click.Path(exists=True, path_type=str), default=None, help="AnnData to write (default: the region's file in the AnnData store).")
@dataset_options
@_options(_STORE_ROOTS)
@dataset_overrides
@table_options
@click.pass_context
def backup_adata_region(
    ctx,
    prefix_name: str,
    suffix: str = "",
    adata_file: str | None = None,
    exp: str | None = None,
    reg: str | None = None,
    zarr_path: str | None = None,
    raw_path: str | None = None,
    seg_path: str | None = None,
    figures_path: str | None = None,
    anndata_path: str | None = None,
    annotations_path: str | None = None,
    **roots,
):
    """
    Write the AnnData file of segmentation PREFIX_NAME back into the region's spatialdata zarr store.
    Can serve to bring annotations made outside spatialdata (i.e. mapmycells) back into the store.
    Needs an env with spatialdata.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        suffix: Table suffix, e.g. _filt.\n
        adata_file: AnnData to write instead of the one in the AnnData store.\n
    """
    from spida.io.layout import table_key
    from spida.I._io import read_adata, backup_table

    d = _dataset(ctx, exp, reg, zarr_path, raw_path, seg_path, figures_path, anndata_path, annotations_path, **roots)
    key = table_key(prefix_name, suffix)
    backup_table(d, read_adata(adata_file or d.table_h5ad(key)), key)


@cli.command(name='backup-adata-experiment', cls=RichCommand)
@click.argument('prefix_name', type=str)
@click.option("--exp", "exp_name", required=True, type=str, help="Experiment name.")
@click.option('--suffix', type=click.STRING, default='', help='Table suffix. Defaults to none.')
@_options(_STORE_ROOTS)
@click.pass_context
def backup_adata_experiment(
    ctx,
    prefix_name: str,
    exp_name: str,
    suffix: str = "",
    **roots,
):
    """
    Backup the AnnData files of segmentation PREFIX_NAME for every region of an entire experiment --exp.

    Parameters:\n
        prefix_name: Segmentation name, e.g. cellpose_cell.\n
        exp_name: Name of the experiment.\n
    """
    from spida.io.layout import TreeLayout, table_key
    from spida.I._io import read_adata, backup_table

    layout = TreeLayout(run_config(ctx, **roots))
    key = table_key(prefix_name, suffix)
    for reg in layout.units(exp_name):
        d = layout.dataset(exp_name, reg)
        backup_table(d, read_adata(d.table_h5ad(key)), key)


