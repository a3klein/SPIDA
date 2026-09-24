"""The contract every per-tile transform implements.

A transform is a module exposing NAME, USES_GPU, DEFAULTS, validate() and run().
The orchestrator supplies a TileContext; everything else is the transform's own
parameters, declared in its run() signature.
"""
from __future__ import annotations

import inspect
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, Sequence


@dataclass(frozen=True)
class TileContext:
    """What the orchestrator knows about the tile being transformed.

    `gpu_lock` serialises access to the GPU and is only supplied to transforms that
    declare USES_GPU; CPU transforms get None and so run in parallel.
    """

    channel: str
    colors: Sequence
    z_step_um: float
    gpu: bool = False
    gpu_lock: Any = None


class Transform(Protocol):
    NAME: str
    USES_GPU: bool
    DEFAULTS: dict

    def validate(self, params: dict) -> dict: ...

    def run(self, src: Path, dst: Path, ctx: TileContext, **params) -> Path: ...


def check_params(run_fn, params: dict, name: str) -> dict:
    """Reject parameters the transform does not accept.

    Derived from run()'s signature rather than a hand-kept list, so the two cannot
    drift. This is why run() takes no **kwargs: a typo in a parameter name must be an
    error, not something silently dropped.
    """
    accepted = {
        p for p in inspect.signature(run_fn).parameters
        if p not in ("src", "dst", "ctx")
    }
    unknown = set(params) - accepted
    if unknown:
        raise ValueError(
            f"{name} does not accept {sorted(unknown)}; accepts {sorted(accepted)}"
        )
    return params
