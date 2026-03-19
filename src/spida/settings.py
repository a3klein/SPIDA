import os
import sys
import logging
from collections.abc import Iterable, Mapping
from typing import Any

# try to use psutil for parent process inspection; fall back to /proc if unavailable
try:
    import psutil  # type: ignore
except Exception:
    psutil = None

# try to import rich's RichHandler if available
try:
    from rich.logging import RichHandler  # type: ignore
    _RICH_AVAILABLE = True
except Exception:
    RichHandler = None
    _RICH_AVAILABLE = False


def _env_flag_true(name: str, default: bool = False) -> bool:
    """Parse boolean-like environment variables.

    Truthy values: 1, true, yes, on
    Falsy values: 0, false, no, off
    """
    value = os.getenv(name)
    if value is None:
        return default

    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    return default


def _get_parent_cmdline() -> str:
    """Return the parent's command line as a single string, best-effort.

    Uses psutil if available, otherwise reads from /proc/<ppid>/cmdline (Linux).
    Returns empty string on failure.
    """
    try:
        ppid = os.getppid()
    except Exception:
        return ""

    # psutil provides a robust API
    if psutil:
        try:
            p = psutil.Process(ppid)
            cmd = p.cmdline()
            return " ".join(cmd) if cmd else ""
        except Exception:
            pass

    # Fallback: read /proc
    try:
        with open(f"/proc/{ppid}/cmdline", "rb") as f:
            data = f.read()
            if not data:
                return ""
            # cmdline entries are null-separated
            parts = [p.decode("utf-8", errors="replace") for p in data.split(b"\x00") if p]
            return " ".join(parts)
    except Exception:
        return ""


def detect_runtime_env() -> str:
    """Detect runtime environment.

    Returns one of:
      - 'slurm-batch'       : running inside a Slurm batch job (sbatch, no TTY)
      - 'slurm-interactive' : running inside Slurm but interactive (srun/salloc or has a TTY)
      - 'interactive'       : not running under Slurm (regular CLI / REPL)

    Detection strategy:
      - If SLURM_* environment variables are absent -> 'interactive'
      - If SLURM_* present and stdout/stderr is a TTY -> 'slurm-interactive'
      - Inspect parent process command line for 'srun', 'salloc', or 'sbatch' when necessary.
      - Default to 'slurm-batch' when SLURM is present but no interactive hints are found.
    """
    slurm_keys = (
        "SLURM_JOB_ID",
        "SLURM_JOB_NODELIST",
        "SLURM_NTASKS",
        "SLURM_PROCID",
        "SLURM_ARRAY_TASK_ID",
    )

    if not any(k in os.environ for k in slurm_keys):
        return "interactive"

    # We're under Slurm (some SLURM_* present). Check for interactivity.
    try:
        if sys.stdout.isatty() or sys.stderr.isatty():
            return "slurm-interactive"
    except Exception:
        pass

    parent_cmd = _get_parent_cmdline().lower()
    if any(x in parent_cmd for x in ("srun", "salloc")):
        return "slurm-interactive"
    if "sbatch" in parent_cmd:
        return "slurm-batch"

    # Fallback default when SLURM env is present
    return "slurm-batch"


def configure_logging_for_runtime(level: int = logging.INFO, logger: logging.Logger | None = None, log_dir: str | None = None) -> str:
    """Configure logging handlers based on detected runtime and return the chosen environment string.

    Behavior:
      - interactive: use RichHandler (if available) or StreamHandler(sys.stdout) with colors
      - slurm-interactive: prefer RichHandler with markup disabled to avoid noisy escape sequences; fallback to plain StreamHandler
      - slurm-batch: write to a file named 'slurm_<JOBID>.log' in log_dir (or cwd) using FileHandler and also attach a plain StreamHandler to stderr

    Returns:
      The detected environment string (same values as detect_runtime_env).
    """
    env = detect_runtime_env()
    root = logger if logger is not None else logging.getLogger()

    # Remove existing handlers to avoid duplicate logs when called multiple times
    for h in list(root.handlers):
        try:
            root.removeHandler(h)
        except Exception:
            pass

    fmt = "[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s - %(message)s"
    datefmt = "%Y-%m-%dT%H:%M:%S%z"
    formatter = logging.Formatter(fmt, datefmt)

    if env == "interactive":
        if _RICH_AVAILABLE and RichHandler is not None:
            rh = RichHandler(rich_tracebacks=True, show_time=False, markup=True)
            rh.setLevel(level)
            root.addHandler(rh)
        else:
            sh = logging.StreamHandler(sys.stdout)
            sh.setLevel(level)
            sh.setFormatter(formatter)
            root.addHandler(sh)

    elif env == "slurm-interactive":
        # clusters often capture stdout/stderr to files; disable rich markup to avoid escape sequences
        if _RICH_AVAILABLE and RichHandler is not None:
            try:
                rh = RichHandler(rich_tracebacks=False, show_time=False, markup=False)
                rh.setLevel(level)
                root.addHandler(rh)
            except Exception:
                sh = logging.StreamHandler(sys.stdout)
                sh.setLevel(level)
                sh.setFormatter(formatter)
                root.addHandler(sh)
        else:
            sh = logging.StreamHandler(sys.stdout)
            sh.setLevel(level)
            sh.setFormatter(formatter)
            root.addHandler(sh)

    else:  # slurm-batch
        jobid = os.environ.get("SLURM_JOB_ID", "unknown")
        # If a log_dir was provided attempt to also write a copy there, otherwise rely on Slurm-captured stdout
        if log_dir:
            try:
                os.makedirs(log_dir, exist_ok=True)
                log_path = os.path.join(log_dir, f"slurm_{jobid}.log")
                fh = logging.FileHandler(log_path)
                fh.setLevel(level)
                fh.setFormatter(formatter)
                root.addHandler(fh)
            except Exception:
                pass

        # IMPORTANT: write to stdout so Slurm's job output file captures the logs
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(level)
        sh.setFormatter(formatter)
        root.addHandler(sh)

    root.setLevel(level)

    # Keep noisy dependency loggers from flooding output while still surfacing warnings/errors.
    logging.getLogger("spatialdata").setLevel(logging.WARNING)

    # Apply warning filters used across SPIDA entry points.
    config_warnings()

    return env




def config_warnings(
    extra_filters: Iterable[Mapping[str, Any]] | None = None,
    reset_existing_filters: bool = False,
) -> None:
    """Configure warning filters for known noisy dependencies.

    Parameters
    ----------
    extra_filters
        Optional iterable of warning filter mappings passed to
        ``warnings.filterwarnings``. This makes the behavior extensible without
        touching SPIDA source each time a new warning appears.
    reset_existing_filters
        If ``True``, call ``warnings.resetwarnings()`` before applying defaults.
    """
    import warnings

    # Env override for deep debugging:
    # SPIDA_STRICT_WARNINGS=1 disables SPIDA's warning suppression filters.
    if _env_flag_true("SPIDA_STRICT_WARNINGS", default=False):
        if reset_existing_filters:
            warnings.resetwarnings()
        return

    if reset_existing_filters:
        warnings.resetwarnings()

    # Optional warning classes from dependencies (safe fallback if unavailable).
    unstable_spec_warning: type[Warning] = Warning
    numba_perf_warning: type[Warning] = Warning
    try:
        from zarr.errors import UnstableSpecificationWarning  # type: ignore

        unstable_spec_warning = UnstableSpecificationWarning
    except Exception:
        pass

    try:
        from numba.core.errors import NumbaPerformanceWarning  # type: ignore

        numba_perf_warning = NumbaPerformanceWarning
    except Exception:
        pass

    default_filters: list[dict[str, Any]] = [
        # Existing broad filters.
        {"action": "ignore", "category": UserWarning, "module": r"^zarr(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^anndata(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^scanpy(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^matplotlib(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^xarray_schema(\.|$)"},
        {"action": "ignore", "category": FutureWarning, "module": r"^dask(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^ome_zarr(\.|$)"},
        {"action": "ignore", "category": SyntaxWarning, "module": r"^leidenalg(\.|$)"},
        {"action": "ignore", "category": UserWarning, "module": r"^xarray_schema(\.|$)"},
        {"action": "ignore", "category": FutureWarning, "module": r"^dask\.dataframe"},
        # {"action": "ignore", "category": RuntimeWarning, "module": r"^scipy\.sparse",
        #     "message": r".*divide by zero encountered in reciprocal.*",
        # },
        {"action": "ignore", "category": UserWarning, "module": r"^libpysal\.weights"},
        {"action": "ignore", "category": UserWarning, "module": r"^spatialdata_plot(\.|$)"},
        {"action": "ignore", "category": RuntimeWarning, "module": r"^scanpy\.tools\._rank_genes_groups$"},
        {"action": "ignore", "category": FutureWarning},
        {"action": "ignore", "category": unstable_spec_warning, "module": r"^zarr(\.|$)"},
        {"action": "ignore", "category": numba_perf_warning, "module": r"^scanpy(\.|$)"},
    ]

    for fw in default_filters:
        warnings.filterwarnings(**fw)

    if extra_filters:
        for fw in extra_filters:
            warnings.filterwarnings(**dict(fw))



# OLD
# def setup_logging(stdout=False, quiet=False, **kwargs):
#     stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
#     stream_handler.setLevel(logging.INFO if not quiet else logging.WARNING)
#     stream_handler.setFormatter(
#         logging.Formatter(
#             "[%(levelname)s | %(name)s | %(module)s | L%(lineno)d] %(asctime)s - %(message)s",
#             datefmt="%Y-%m-%dT%H:%M:%S%z",
#         )
#     )

#     logger.addHandler(stream_handler)
#     logger.setLevel(logging.INFO if not quiet else logging.WARNING)
