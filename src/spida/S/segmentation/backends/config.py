"""Shared across segmentation method backends. Contains run configurations

Every segmentation method owns a ``{Method}Config`` dataclass describing the
parameters its backend accepts. The classes deliberately do not share
fields — cellpose and proseg have almost nothing in common — but they share
this base so the mechanics are identical across methods:

* :meth:`BackendConfig.from_kwargs` builds a config from loose CLI kwargs,
  applying deprecated aliases and **rejecting unknown keys**. A typo must fail
  loudly rather than be ignored.
* :meth:`BackendConfig.validate` holds consistency rules.
* :meth:`BackendConfig.describe` renders the parameter table shown by
  ``segment-region <method> --list-params``, so the CLI stays small while the
  full parameter surface is discoverable.
* :meth:`BackendConfig.to_meta` produces a JSON-safe dict for run provenance.

To add a method: subclass :class:`BackendConfig`, set ``METHOD``, declare the
fields, and fill in ``HELP`` / ``CHOICES`` / ``ALIASES``. Nothing else is
required — the CLI and the validation machinery pick it up automatically.
"""

from __future__ import annotations

import difflib
import logging
import typing
import warnings
from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Any, ClassVar

logger = logging.getLogger(__package__)


class ConfigError(ValueError):
    """Raised when a backend run configuration is invalid.

    Deliberately a ``ValueError`` subclass so existing ``except ValueError``
    handlers keep working, but distinguishable when callers want to be precise.
    """


def _type_options(annotation: Any) -> tuple[list[type], list[Any]]:
    """Split a type annotation into concrete types and ``Literal`` values.

    ``int | Literal["auto"] | None`` -> ``([int, NoneType], ["auto"])``. Used by
    the coercion step to decide what a raw CLI string is allowed to become.
    """
    concrete: list[type] = []
    literals: list[Any] = []

    def _walk(ann: Any) -> None:
        origin = typing.get_origin(ann)
        if origin is typing.Literal:
            literals.extend(typing.get_args(ann))
        elif origin is not None and origin.__name__ in ("UnionType", "Union"):
            for arg in typing.get_args(ann):
                _walk(arg)
        elif isinstance(ann, type):
            concrete.append(ann)

    # `X | Y` unions land here as types.UnionType, which get_origin reports oddly
    # across versions; normalise by inspecting get_args when present.
    args = typing.get_args(annotation)
    if args and typing.get_origin(annotation) is not typing.Literal:
        for arg in args:
            _walk(arg)
    else:
        _walk(annotation)
    return concrete, literals


@dataclass(frozen=True)
class BackendConfig:
    """Base class for per-method run configurations. Not used directly."""

    #: Method name as registered in the backend registry (e.g. "cellpose").
    METHOD: ClassVar[str] = "base"
    #: Deprecated parameter name -> current name. Applied with a warning.
    ALIASES: ClassVar[dict[str, str]] = {}
    #: Field name -> one-line help string, shown by ``describe()``.
    HELP: ClassVar[dict[str, str]] = {}
    #: Field name -> tuple of permitted values, enforced in ``validate()``.
    CHOICES: ClassVar[dict[str, tuple]] = {}

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def field_names(cls) -> tuple[str, ...]:
        return tuple(f.name for f in fields(cls))

    @classmethod
    def _hints(cls) -> dict[str, Any]:
        # resolved lazily: `from __future__ import annotations` makes raw
        # dataclass field types strings, which we cannot introspect directly.
        return typing.get_type_hints(cls)

    @classmethod
    def _coerce(cls, name: str, value: Any, annotation: Any) -> Any:
        """Best-effort coercion of a raw CLI value to the declared field type.

        ``parse_click_kwargs`` already ``literal_eval``s most values, so this
        mainly handles the leftovers (e.g. ``--scale=4.0`` for an int field, or
        a quoted ``"3d_true"``). Anything that cannot be coerced is returned
        unchanged and left for :meth:`validate` to reject with a better message.
        """
        concrete, literals = _type_options(annotation)
        if value in literals:
            return value
        if value is None:
            return None
        if bool in concrete and isinstance(value, bool):
            return value
        if int in concrete and not isinstance(value, bool):
            if isinstance(value, int):
                return value
            if isinstance(value, float) and float(value).is_integer():
                return int(value)
            if isinstance(value, str):
                try:
                    return int(value)
                except ValueError:
                    pass
        if float in concrete and isinstance(value, (int, float)) and not isinstance(value, bool):
            return float(value)
        if float in concrete and isinstance(value, str):
            try:
                return float(value)
            except ValueError:
                pass
        return value

    @classmethod
    def _check_type(cls, name: str, value: Any, annotation: Any) -> None:
        """Reject a value whose type cannot match the declared field type.

        Needed because ``parse_click_kwargs`` requires ``--key value`` form: a
        bare boolean flag such as ``--apply_clahe`` swallows the *next* token as
        its value, so without this check ``apply_clahe`` could end up holding the
        string ``'--foo'`` and evaluate as truthy.

        Fields listed in ``CHOICES`` are skipped so ``validate`` can report them
        with its more specific message.
        """
        if name in cls.CHOICES:
            return
        concrete, literals = _type_options(annotation)
        if not concrete and not literals:
            return                                   # unannotated / Any
        if literals and value in literals:
            return

        allowed = tuple(t for t in concrete if t is not type(None))
        if value is None:
            if type(None) in concrete:
                return
            raise ConfigError(f"{cls.__name__}: {name} may not be None.")
        if not allowed:
            raise ConfigError(
                f"{cls.__name__}: {name}={value!r} is not one of "
                f"{', '.join(repr(literal) for literal in literals)}."
            )

        # bool is a subclass of int, so handle it explicitly in both directions:
        # a bool must not satisfy an int/float field, nor an int a bool field.
        if bool in allowed:
            if isinstance(value, bool):
                return
        elif isinstance(value, bool):
            raise ConfigError(
                f"{cls.__name__}: {name}={value!r} is a bool but "
                f"{' | '.join(t.__name__ for t in allowed)} was expected."
            )
        if float in allowed and isinstance(value, int):
            return                                   # ints widen to float
        if isinstance(value, allowed):
            return

        expected = " | ".join(t.__name__ for t in allowed)
        if literals:
            expected += " | " + " | ".join(repr(literal) for literal in literals)
        hint = ""
        if bool in allowed:
            hint = (f" Booleans must be passed as --{name}=true / --{name}=false; "
                    f"a bare --{name} consumes the next argument as its value.")
        raise ConfigError(
            f"{cls.__name__}: {name}={value!r} has type "
            f"{type(value).__name__}; expected {expected}.{hint}"
        )

    @classmethod
    def from_kwargs(cls, **kwargs: Any) -> "BackendConfig":
        """Build a config from loose kwargs, rejecting anything unrecognised.

        Deprecated names in ``ALIASES`` are translated with a
        ``DeprecationWarning``. Unknown names raise :class:`ConfigError` with
        the closest valid names suggested.
        """
        valid = set(cls.field_names())
        resolved: dict[str, Any] = {}

        for key, value in kwargs.items():
            target = key
            if key in cls.ALIASES:
                target = cls.ALIASES[key]
                warnings.warn(
                    f"{cls.__name__}: '{key}' is deprecated; use '{target}'.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                logger.warning("%s: parameter '%s' is deprecated; using '%s'",
                               cls.__name__, key, target)
            if target not in valid:
                close = difflib.get_close_matches(target, sorted(valid), n=3)
                hint = f" Did you mean: {', '.join(close)}?" if close else ""
                raise ConfigError(
                    f"{cls.__name__}: unknown parameter '{key}'.{hint}\n"
                    f"Run `segment-region {cls.METHOD} --list-params` "
                    f"to see every accepted parameter."
                )
            if target in resolved:
                raise ConfigError(
                    f"{cls.__name__}: '{target}' supplied twice "
                    f"(once via the deprecated alias '{key}')."
                )
            resolved[target] = value

        hints = cls._hints()
        coerced = {k: cls._coerce(k, v, hints.get(k, Any)) for k, v in resolved.items()}
        for key, value in coerced.items():
            cls._check_type(key, value, hints.get(key, Any))

        cfg = cls(**coerced)
        cfg.validate()
        return cfg

    # ------------------------------------------------------------------
    # validation
    # ------------------------------------------------------------------
    def validate(self) -> None:
        """Check ``CHOICES`` membership. Subclasses extend with their own rules."""
        for name, allowed in self.CHOICES.items():
            value = getattr(self, name)
            if value not in allowed:
                raise ConfigError(
                    f"{type(self).__name__}: {name}={value!r} is not one of "
                    f"{', '.join(repr(a) for a in allowed)}."
                )

    # ------------------------------------------------------------------
    # presentation / provenance
    # ------------------------------------------------------------------
    @classmethod
    def describe(cls) -> str:
        """Render the human-readable parameter table for ``--list-params``."""
        hints = cls._hints()
        lines = [
            f"{cls.__name__} — parameters accepted by `segment-region {cls.METHOD}`",
            "",
            "Pass any of these as `--name=value` after the positional arguments.",
            "Unknown names are rejected rather than ignored.",
            "",
            "Booleans need an explicit value: `--apply_clahe=true`. A bare",
            "`--apply_clahe` consumes the NEXT argument as its value.",
            "",
        ]
        width = max((len(n) for n in cls.field_names()), default=0)
        for f in fields(cls):
            concrete, literals = _type_options(hints.get(f.name, Any))
            if literals:
                type_str = " | ".join(repr(literal) for literal in literals)
                if concrete:
                    names = [t.__name__ for t in concrete if t is not type(None)]
                    if names:
                        type_str = " | ".join(names + [type_str])
            else:
                names = [t.__name__ for t in concrete if t is not type(None)]
                type_str = " | ".join(names) if names else "any"
            lines.append(f"  {f.name.ljust(width)}  {type_str}  (default: {f.default!r})")
            if f.name in cls.HELP:
                lines.append(f"  {' ' * width}      {cls.HELP[f.name]}")
            if f.name in cls.CHOICES:
                allowed = ", ".join(repr(a) for a in cls.CHOICES[f.name])
                lines.append(f"  {' ' * width}      choices: {allowed}")
            lines.append("")
        if cls.ALIASES:
            lines.append("Deprecated aliases (accepted, with a warning):")
            for old, new in sorted(cls.ALIASES.items()):
                lines.append(f"  {old} -> {new}")
            lines.append("")
        return "\n".join(lines)

    def to_meta(self) -> dict[str, Any]:
        """JSON-safe dict of the resolved configuration, for run provenance."""
        out: dict[str, Any] = {}
        for key, value in asdict(self).items():
            out[key] = str(value) if isinstance(value, Path) else value
        return out
