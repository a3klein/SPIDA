import os
import json
import click 
from rich_click import RichCommand, RichGroup, echo as rich_echo # type: ignore
from pathlib import Path
from spida.utilities.script_utils import parse_click_kwargs
from dotenv import load_dotenv  # type: ignore


RENAME_CONFIG_KEYS = {
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
    CLI options > config file > environment variables > defaults.

    The CLI layer is handled by Click automatically via default_map; this
    function resolves defaults, env, and config file values only.
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


@click.group(cls=RichGroup, help="Configuration management commands for SPIDA.")
@click.pass_context
def config_manager(ctx):
    ctx.ensure_object(dict)
    pass
    

try: 
    load_dotenv()
except Exception as e:
    click.echo(f"Error loading .env file: {e}")
    raise RuntimeError("Failed to load .env file. Please ensure it exists and is correctly formatted.")

@click.command(cls=RichCommand, help="Set up a configuration file with specified or default values.")
@click.option(
    "config_store_path",
    "--config_store_path",
    default=".env",
    prompt=True,
    type=click.Path(), 
    help="Path to the configuration store file (default: .env)"
)
@click.option(
    "rust_bin_path",
    "--rust_bin_path",
    default=lambda: _default_for("RUST_BIN_PATH"),
    show_default="RUST_BIN_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Path to the Rust binary (default: None, uses RUST_BIN_PATH in .env)"
)
@click.option(
    "vpt_bin_path",
    "--vpt_bin_path",
    default=lambda: _default_for("VPT_BIN_PATH"),
    show_default="VPT_BIN_PATH in .env",
    prompt=True,
    type=click.Path(), 
    help="Path to the VPT binary (default: None, uses VPT_BIN_PATH in .env)"
)
@click.option(
    "deconwolf_config_file",
    "--deconwolf_config_file",
    default=lambda: _default_for("DECONWOLF_CONFIG"),
    show_default="DECONWOLF_CONFIG in .env",
    prompt=True,
    type=click.Path(),
    help="Path to the Deconwolf configuration file (default: None, uses DECONWOLF_CONFIG in .env)"
)
@click.option(
    "zarr_storage_path",
    "--zarr_storage_path",
    default=lambda: _default_for("ZARR_STORAGE_PATH"),
    show_default="ZARR_STORAGE_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Path to the Zarr storage (default: None, uses ZARR_STORAGE_PATH in .env)"
)
@click.option(
    "root_path",
    "--root_path",
    default=lambda: _default_for("PROCESSED_ROOT_PATH"),
    show_default="PROCESSED_ROOT_PATH in .env",
    prompt=True, 
    type=click.Path(),
    help="Root path for processed data (default: None, uses PROCESSED_ROOT_PATH in .env)"
)
@click.option(
    "segmentation_output_dir",
    "--segmentation_output_dir",
    default=lambda: _default_for("SEGMENTATION_OUT_PATH"),
    show_default="SEGMENTATION_OUT_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Directory for segmentation outputs (default: None, uses SEGMENTATION_OUT_PATH in .env)"
)
@click.option(
    "anndata_store_dir",
    "--anndata_store_dir",
    default=lambda: _default_for("ANNDATA_STORE_PATH"),
    show_default="ANNDATA_STORE_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Directory for AnnData storage (default: None, uses ANNDATA_STORE_PATH in .env)"
)
@click.option(
    "annotation_store_dir",
    "--annotation_store_dir",
    default=lambda: _default_for("ANNOTATION_STORE_PATH"),
    show_default="ANNOTATION_STORE_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Directory for annotation storage (default: None, uses ANNOTATION_STORE_PATH in .env)"
)
@click.option(
    "image_store_path",
    "--image_store_path",
    default=lambda: _default_for("IMAGE_STORE_PATH"),
    show_default="IMAGE_STORE_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Path to the image storage (default: None, uses IMAGE_STORE_PATH in .env)"
)
@click.option(
    "cutoffs_path",
    "--cutoffs_path",
    default=lambda: _default_for("CUTOFFS_PATH"),
    show_default="CUTOFFS_PATH in .env",
    prompt=True,
    type=click.Path(),
    help="Path to the cutoffs file (default: None, uses CUTOFFS_PATH in .env)"
)
@click.option(
    "overwrite",
    "--overwrite",
    is_flag=True,
    default=False,
    prompt=True,
    help="Overwrite existing configuration file if it exists."
)
@click.option(
    "ext_type",
    "--ext_type",
    default="env",
    prompt=True,
    type=click.Choice(["env", "json"]),
    help="File extension type for the configuration file (default: env)."
)
def setup_config(
    config_store_path: str | Path = ".env",
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
    overwrite: bool = False,
    ext_type: str = "env",
):
    """
    Sets up the configuration file with provided or default values. 
    This is for setting up different projects without having to modify the .env file manually.
    """

    config = resolve_config(
        defaults=DEFAULT_CONFIG,
        config_path=None,
    )
    # CLI provided values should override resolved defaults
    config.update(
        {
            "RUST_BIN_PATH": rust_bin_path or config["RUST_BIN_PATH"],
            "VPT_BIN_PATH": vpt_bin_path or config["VPT_BIN_PATH"],
            "DECONWOLF_CONFIG": deconwolf_config_file or config["DECONWOLF_CONFIG"],
            "ZARR_STORAGE_PATH": zarr_storage_path or config["ZARR_STORAGE_PATH"],
            "PROCESSED_ROOT_PATH": root_path or config["PROCESSED_ROOT_PATH"],
            "SEGMENTATION_OUT_PATH": segmentation_output_dir or config["SEGMENTATION_OUT_PATH"],
            "ANNDATA_STORE_PATH": anndata_store_dir or config["ANNDATA_STORE_PATH"],
            "ANNOTATION_STORE_PATH": annotation_store_dir or config["ANNOTATION_STORE_PATH"],
            "IMAGE_STORE_PATH": image_store_path or config["IMAGE_STORE_PATH"],
            "CUTOFFS_PATH": cutoffs_path or config["CUTOFFS_PATH"],
        }
    )

    config_path = Path(config_store_path)
    if config_path.exists() and not overwrite:
        click.echo(f"Configuration file already exists at {config_path.resolve()}. Use --overwrite to replace it.")
        return
    if ext_type == "env":
        if not str(config_store_path).endswith(".env"):
            config_store_path = f"{config_store_path}.env"
        config_path = Path(config_store_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        with config_path.open("a") as f:
            for key, value in config.items():
                f.write(f"{key}={value}\n")
    elif ext_type == "json":
        if not str(config_store_path).endswith(".json"):
            config_store_path = f"{config_store_path}.json"
        config_path = Path(config_store_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        with config_path.open("a") as f:
            json.dump(config, f, indent=4)
    else: # TODO add yaml of .ini for peoples preferences
        click.echo("Unsupported file extension type. Use 'env' or 'json'.")
        return
    
    click.echo(f"Configuration file created at {config_path.resolve()}")

def load_config(
    config_path : str | Path  = ".env",
): 
    """
    Loads the configuration from the specified CONFIG_PATH file and returns it as a dictionary. 
    Parameters:
    - config_path (str | Path): The path to the configuration file (default: ".env").
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_path.resolve()}.")
    
    config = {}
    if str(config_path).endswith(".env"):
        with config_path.open("r") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    key, value = line.strip().split("=", 1)
                    config[key] = value
    elif str(config_path).endswith(".json"):
        with config_path.open("r") as f:
            config = json.load(f)
    else:
        raise ValueError("Unsupported configuration file format. Use '.env' or '.json'.")
    return config

def load_config_into_env(config): 
    for key, value in config.items():
        os.environ[key] = value

@click.command(cls=RichCommand, help="Display the current configuration settings.")
@click.argument(
    "config_path",
    type=click.Path(exists=True),
    default=".env"
)
def display_config(
    config_path : str | Path = ".env",
):
    """Displays the current configuration settings from CONFIG_PATH."""
    config = load_config(config_path)
    for key, value in config.items():
        click.echo(click.style(key, fg='magenta', bold=True) + f": {click.style(value, fg='green')}")


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
        defaults = {}
        for key, value in resolved.items():
            defaults[RENAME_CONFIG_KEYS.get(key, key.lower())] = value

        # build a per-subcommand default_map so defaults apply to every command
        mapping = {name: defaults.copy() for name in self.commands}
        extra['default_map'] = mapping

        # create the context using the parent implementation
        ctx = super().make_context(info_name, args, parent=parent, **extra)

        # Ensure ctx.obj exists
        try:
            ctx.ensure_object(dict)
        except Exception:
            if ctx.obj is None:
                ctx.obj = {}

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