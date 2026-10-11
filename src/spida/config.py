import os
import json
import click
from rich_click import RichCommand, RichGroup, echo as rich_echo # type: ignore
from dataclasses import dataclass
from pathlib import Path
from spida.utilities.script_utils import parse_click_kwargs


# config key -> field name
CONFIG_KEY_TO_FIELD = {
    "RUST_BIN_PATH" : "rust_bin_path",
    "VPT_BIN_PATH": "vpt_bin_path",
    "DECONWOLF_CONFIG": "deconwolf_config_file",
    "ZARR_STORAGE_PATH": "zarr_store",
    "PROCESSED_ROOT_PATH": "root_path",
    "SEGMENTATION_OUT_PATH": "segmentation_store",
    "ANNDATA_STORE_PATH": "anndata_store",
    "ANNOTATION_STORE_PATH": "annotation_store",
    "ANNOTATIONS_STORE_PATH": "annotation_store",
    "IMAGE_STORE_PATH": "image_store",
    "CUTOFFS_PATH": "cutoffs_path",
    "SPIDA_SITE_DIR": "site_dir",
    "MMC_DIR": "mmc_store_path",
    "GENE_PANEL_PATH": "gene_panel_path",
}

DEFAULT_CONFIG = {
    "RUST_BIN_PATH": "/path/to/rust/bin",
    "VPT_BIN_PATH": "/path/to/vpt/bin",
    "DECONWOLF_CONFIG": "/path/to/deconwolf/config/file",
    "ZARR_STORAGE_PATH": "/path/to/zarr/storage",
    "PROCESSED_ROOT_PATH": "/path/to/processed/root",
    "SEGMENTATION_OUT_PATH": "/path/to/segmentation/output",
    "ANNDATA_STORE_PATH": "/path/to/anndata/store",
    "ANNOTATION_STORE_PATH": "/path/to/annotation/store",
    "IMAGE_STORE_PATH": "/path/to/image/store",
    "CUTOFFS_PATH": "/path/to/cutoffs.json",
    "SPIDA_SITE_DIR": "/path/to/spida/site",
    # Optional.
    "MMC_DIR": "",
    "GENE_PANEL_PATH": "",
    # Optional. Unset -> CELLPOSE_LOCAL_MODELS_PATH, else cellpose's ~/.cellpose/models
    # (resolved in backends/cellpose_config.resolve_model_dir).
    "CELLPOSE_MODEL_PATH": "",
}

CONFIG_DESCRIPTIONS = {
    "RUST_BIN_PATH": "Path to the Rust binary used by SPIDA.",
    "VPT_BIN_PATH": "Path to the VPT binary used by SPIDA.",
    "DECONWOLF_CONFIG": "Path to the Deconwolf configuration file.",
    "ZARR_STORAGE_PATH": "Root directory for Zarr stores.",
    "PROCESSED_ROOT_PATH": "Root directory for processed outputs.",
    "SEGMENTATION_OUT_PATH": "Directory for segmentation outputs.",
    "ANNDATA_STORE_PATH": "Directory for AnnData (.h5ad) storage.",
    "ANNOTATION_STORE_PATH": "Directory for annotation artifacts/labels.",
    "IMAGE_STORE_PATH": "Directory for image outputs/exports.",
    "CUTOFFS_PATH": "Path to QC/cutoffs configuration file.",
    "SPIDA_SITE_DIR": "Path to the SPIDA website directory.",
    "MMC_DIR": "Optional. Root of the MapMyCells reference stores.",
    "GENE_PANEL_PATH": "Optional. Directory of gene-panel name mappings used by MapMyCells.",
    "CELLPOSE_MODEL_PATH": (
        "Optional. Directory holding cellpose model weights (downloads land here too). "
        "Defaults to CELLPOSE_LOCAL_MODELS_PATH, else ~/.cellpose/models. Override "
        "per-run with `segment-region cellpose ... --model_dir=<path>`."
    ),
}

def _normalize_config_keys(config: dict) -> dict:
    """Normalize known alias keys into canonical names."""
    if "ANNOTATIONS_STORE_PATH" in config and "ANNOTATION_STORE_PATH" not in config:
        config["ANNOTATION_STORE_PATH"] = config["ANNOTATIONS_STORE_PATH"]
    return config


def _config_from_env(defaults: dict) -> dict:
    """Overlay defaults with environment variables if present."""
    resolved = defaults.copy()
    for key in defaults:
        env_val = os.getenv(key)
        if env_val is not None:
            resolved[key] = env_val
    # allow alternative env var spelling for annotations store
    env_alt = os.getenv("ANNOTATIONS_STORE_PATH")
    if env_alt is not None:
        resolved["ANNOTATION_STORE_PATH"] = env_alt
    return resolved


def resolve_config(
    config_path: str | Path | None = None,
    defaults: dict | None = None,
) -> dict:
    """
    Resolve configuration with precedence:
    config file > environment variables > defaults (placeholders, i.e. unset).

    CLI options are applied on top per command, by ``spida.io.options.run_config``.
    """
    base_defaults = defaults or DEFAULT_CONFIG
    resolved = _config_from_env(base_defaults)

    if config_path:
        cfg = load_config(config_path)
        cfg = _normalize_config_keys(cfg)
        resolved.update(cfg)

    return resolved


def _default_for(key: str) -> str:
    """Return a default value using env > DEFAULT_CONFIG."""
    resolved = resolve_config(defaults=DEFAULT_CONFIG)
    return resolved.get(key, DEFAULT_CONFIG.get(key, ""))


def field_for(key: str) -> str:
    return CONFIG_KEY_TO_FIELD.get(key, key.lower())


@dataclass(frozen=True)
class SpidaConfig:
    """Resolved global settings (locations and binaries) for one run. Unset paths are None."""
    zarr_store: Path | None = None
    root_path: Path | None = None
    segmentation_store: Path | None = None
    anndata_store: Path | None = None
    annotation_store: Path | None = None
    image_store: Path | None = None
    cutoffs_path: Path | None = None
    site_dir: Path | None = None
    mmc_store_path: Path | None = None
    gene_panel_path: Path | None = None
    rust_bin_path: Path | None = None
    vpt_bin_path: Path | None = None
    deconwolf_config_file: Path | None = None
    cellpose_model_path: Path | None = None

    def require(self, name: str) -> Path:
        """Return a path setting, or raise naming the config key that sets it."""
        value = getattr(self, name)
        if value is None:
            key = next(k for k in DEFAULT_CONFIG if field_for(k) == name)
            raise ValueError(f"{name} is not configured; set {key} in the --config JSON or the environment")
        return value


def _as_path(value) -> Path | None:
    if value is None or str(value).strip() == "" or str(value).startswith("/path/to/"):
        return None
    return Path(os.path.expanduser(str(value)))


def load(config_path: str | Path | None = None, resolved: dict | None = None, **overrides) -> SpidaConfig:
    """
    SpidaConfig from ``resolve_config`` (config file > environment > defaults), then overrides.

    ``resolved`` reuses an already-resolved dict. Overrides use field names; None values are ignored.
    """
    resolved = resolve_config(config_path=config_path) if resolved is None else resolved
    values = {field_for(k): _as_path(resolved.get(k)) for k in DEFAULT_CONFIG}
    values.update({k: _as_path(v) for k, v in overrides.items() if v is not None and k in values})
    return SpidaConfig(**values)


@click.group(cls=RichGroup, help="Configuration management commands for SPIDA.")
@click.pass_context
def config_manager(ctx):
    ctx.ensure_object(dict)
    pass
    

@click.command(cls=RichCommand, help="Write a JSON config file from given values, else the environment.")
@click.option("config_store_path", "--config_store_path", default="spida_config.json", prompt=True, type=click.Path(), help="JSON file to write (default: spida_config.json)")
@click.option("rust_bin_path", "--rust_bin_path", default=lambda: _default_for("RUST_BIN_PATH"), show_default="RUST_BIN_PATH from the environment", prompt=True, type=click.Path(), help="Path to the Rust binary (default: RUST_BIN_PATH from the environment)")
@click.option("vpt_bin_path", "--vpt_bin_path", default=lambda: _default_for("VPT_BIN_PATH"), show_default="VPT_BIN_PATH from the environment", prompt=True, type=click.Path(), help="Path to the VPT binary (default: VPT_BIN_PATH from the environment)")
@click.option("deconwolf_config_file", "--deconwolf_config_file", default=lambda: _default_for("DECONWOLF_CONFIG"), show_default="DECONWOLF_CONFIG from the environment", prompt=True, type=click.Path(), help="Path to the Deconwolf configuration file (default: DECONWOLF_CONFIG from the environment)")
@click.option("zarr_storage_path", "--zarr_storage_path", default=lambda: _default_for("ZARR_STORAGE_PATH"), show_default="ZARR_STORAGE_PATH from the environment", prompt=True, type=click.Path(), help="Path to the Zarr storage (default: ZARR_STORAGE_PATH from the environment)")
@click.option("root_path", "--root_path", default=lambda: _default_for("PROCESSED_ROOT_PATH"), show_default="PROCESSED_ROOT_PATH from the environment", prompt=True, type=click.Path(), help="Root path for processed data (default: PROCESSED_ROOT_PATH from the environment)")
@click.option("segmentation_output_dir", "--segmentation_output_dir", default=lambda: _default_for("SEGMENTATION_OUT_PATH"), show_default="SEGMENTATION_OUT_PATH from the environment", prompt=True, type=click.Path(), help="Directory for segmentation outputs (default: SEGMENTATION_OUT_PATH from the environment)")
@click.option("anndata_store_dir", "--anndata_store_dir", default=lambda: _default_for("ANNDATA_STORE_PATH"), show_default="ANNDATA_STORE_PATH from the environment", prompt=True, type=click.Path(), help="Directory for AnnData storage (default: ANNDATA_STORE_PATH from the environment)")
@click.option("annotation_store_dir", "--annotation_store_dir", default=lambda: _default_for("ANNOTATION_STORE_PATH"), show_default="ANNOTATION_STORE_PATH from the environment", prompt=True, type=click.Path(), help="Directory for annotation storage (default: ANNOTATION_STORE_PATH from the environment)")
@click.option("image_store_path", "--image_store_path", default=lambda: _default_for("IMAGE_STORE_PATH"), show_default="IMAGE_STORE_PATH from the environment", prompt=True, type=click.Path(), help="Path to the image storage (default: IMAGE_STORE_PATH from the environment)")
@click.option("cutoffs_path", "--cutoffs_path", default=lambda: _default_for("CUTOFFS_PATH"), show_default="CUTOFFS_PATH from the environment", prompt=True, type=click.Path(), help="Path to the cutoffs file (default: CUTOFFS_PATH from the environment)")
@click.option("spida_site_dir", "--spida_site_dir", default=lambda: _default_for("SPIDA_SITE_DIR"), show_default="SPIDA_SITE_DIR from the environment", prompt=True, type=click.Path(), help="Path to the SPIDA website directory (default: SPIDA_SITE_DIR from the environment)")
@click.option("overwrite", "--overwrite", is_flag=True, default=False, prompt=True, help="Overwrite existing configuration file if it exists.")
def setup_config(
    config_store_path: str | Path = "spida_config.json",
    rust_bin_path: str | Path | None = None,
    vpt_bin_path: str | Path | None = None,
    deconwolf_config_file: str | Path | None = None,
    zarr_storage_path: str | Path | None = None,
    root_path: str | Path | None = None,
    segmentation_output_dir: str | Path | None = None,
    anndata_store_dir: str | Path | None = None,
    annotation_store_dir: str | Path | None = None,
    image_store_path: str | Path | None = None,
    cutoffs_path: str | Path | None = None,
    spida_site_dir: str | Path | None = None,
    overwrite: bool = False,
):
    """Write a project JSON config (pass it to any SPIDA command with --config)."""
    config = resolve_config(defaults=DEFAULT_CONFIG, config_path=None)
    given = {"RUST_BIN_PATH": rust_bin_path, "VPT_BIN_PATH": vpt_bin_path, "DECONWOLF_CONFIG": deconwolf_config_file, "ZARR_STORAGE_PATH": zarr_storage_path, "PROCESSED_ROOT_PATH": root_path, "SEGMENTATION_OUT_PATH": segmentation_output_dir, "ANNDATA_STORE_PATH": anndata_store_dir, "ANNOTATION_STORE_PATH": annotation_store_dir, "IMAGE_STORE_PATH": image_store_path, "CUTOFFS_PATH": cutoffs_path, "SPIDA_SITE_DIR": spida_site_dir}
    config.update({k: v for k, v in given.items() if v})

    config_path = Path(config_store_path)
    if config_path.suffix != ".json":
        config_path = config_path.with_suffix(".json")
    if config_path.exists() and not overwrite:
        click.echo(f"Configuration file already exists at {config_path.resolve()}. Use --overwrite to replace it.")
        return
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with config_path.open("w") as f:
        json.dump(config, f, indent=4)
    click.echo(f"Configuration file created at {config_path.resolve()}")

def load_config(config_path: str | Path) -> dict:
    """The JSON config file at ``config_path`` as a dict (UPPER_CASE keys)."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_path.resolve()}.")
    if config_path.suffix != ".json":
        raise ValueError(f"{config_path}: SPIDA config files are JSON (.env is no longer read).")
    with config_path.open("r") as f:
        return json.load(f)

@click.command(cls=RichCommand, help="Display the current configuration settings.")
@click.argument("config_path", type=click.Path(exists=True, dir_okay=False))
def display_config(config_path: str | Path):
    """Displays the current configuration settings from CONFIG_PATH."""
    config = _normalize_config_keys(load_config(config_path))
    for key, value in config.items():
        desc = CONFIG_DESCRIPTIONS.get(key, "")
        suffix = f" — {desc}" if desc else ""
        click.echo(click.style(key, fg='magenta', bold=True) + f": {click.style(value, fg='green')}{suffix}")


# Custom Click Group to preload config settings for all subcommands in group
class ConfigDefaultGroup(RichGroup):
    
    def make_context(self, info_name, args, parent=None, **extra):
        # quick pre-parse for --config / -c
        config_path = None
        for i, a in enumerate(args):
            if a.startswith('--config='):
                config_path = a.split('=', 1)[1]; break
            if a in ('--config', '-c') and i + 1 < len(args):
                config_path = args[i + 1]; break

        resolved = resolve_config(config_path=config_path, defaults=DEFAULT_CONFIG)

        # create the context using the parent implementation
        ctx = super().make_context(info_name, args, parent=parent, **extra)

        # Ensure ctx.obj exists
        try:
            ctx.ensure_object(dict)
        except Exception:
            if ctx.obj is None:
                ctx.obj = {}
        ctx.obj["config"] = load(resolved=resolved)

        # Robust extraction of extra args: find the invoked subcommand in the original
        # argv (`args`) and take the tokens after it as the raw extra args. This avoids
        # relying on click internals like ctx.args which may be empty depending on parsing.
        try:
            subcmd_index = None
            subcmd_name = None
            for i, token in enumerate(args):
                # skip option tokens
                if token.startswith("-"):
                    continue
                # first non-option token that matches a registered command name
                if token in self.commands:
                    subcmd_index = i
                    subcmd_name = token
                    break

            if subcmd_index is not None and subcmd_name is not None:
                # get the Command object
                cmd_obj = self.commands.get(subcmd_name)
                cmd_cs = getattr(cmd_obj, 'context_settings', {}) or {}

                if cmd_cs.get('allow_extra_args') or cmd_cs.get('ignore_unknown_options'):
                    raw_extra = args[subcmd_index + 1 :]
                    if raw_extra:
                        try:
                            parsed = parse_click_kwargs(raw_extra)
                        except Exception:
                            parsed = {'_raw_extra_args': raw_extra}
                        ctx.obj['extra_kwargs'] = parsed
        except Exception:
            # Don't raise from the group context creation if parsing extra args fails.
            pass

        return ctx
    
    # Function for allowing multiple command names (aliases) for a single function
    def command(self, *args, **kwargs):
        """Adds the ability to add `aliases` to commands.

        Usage: @group.command(name='real_name', aliases=['alias1', 'alias2'])
        """
        aliases = kwargs.pop("aliases", None)

        def decorator(f):
            # register the base command first
            cmd = super(ConfigDefaultGroup, self).command(*args, **kwargs)(f)

            if aliases:
                if not isinstance(aliases, (list, tuple)):
                    raise click.UsageError("`aliases` must be a list or tuple of strings.")

                for alias in aliases:
                    # Create a new Command object for the alias so we do not mutate
                    # the original command's name attribute when adding multiple names.
                    alias_cmd = cmd.__class__(
                        alias,
                        params=cmd.params,
                        callback=cmd.callback,
                        help=f"Alias for '{cmd.name}'.\n\n{(cmd.help or '')}",
                        hidden=True,
                    )
                    # add alias command under the alias name
                    self.add_command(alias_cmd, alias)

            return cmd

        return decorator

config_manager.add_command(setup_config)
config_manager.add_command(display_config)

if __name__ == "__main__":
    config_manager()