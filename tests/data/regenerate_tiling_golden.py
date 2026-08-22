#!/usr/bin/env python
"""Regenerate ``tiling_golden.pkl`` — the golden master for ``utilities/tiling.py``.

The golden master pins ``tile_grid_with_overlap`` across its full parameter
surface so refactors of the shared tiling geometry stay behaviour-preserving for
its three consumers (deconvolution tiling, deconvolution reconstruction, and
segmentation polygon merging).

Run this ONLY when a change to the tiling geometry is deliberate. It refuses to
overwrite unless ``--accept`` is passed, prints exactly which configs move, and
by default asserts that no config gains a negative origin — the failure mode the
2026-08-04 clamp was added to eliminate.

    pixi run -e preprocessing python tests/data/regenerate_tiling_golden.py
    pixi run -e preprocessing python tests/data/regenerate_tiling_golden.py --accept
"""

from __future__ import annotations

import argparse
import ast
import itertools
import pickle
import sys
from pathlib import Path

GOLDEN = Path(__file__).with_name("tiling_golden.pkl")

#: Fallback parameter surface, used only when there is no golden file to replay.
#: Deliberately includes degenerate cases (tile larger than the extent, overlap
#: >= tile, non-square tiles, tiny arrays) — those are where the geometry has
#: actually gone wrong historically.
SHAPES = [(100, 100), (1200, 800), (1000, 1000), (37, 53), (512, 512),
          (7, 512, 512), (3, 100, 100)]
TILES = [64, 128, 300, (100, 50), (10, 20)]
OVERLAPS = [0, 16, 64, (10, 20)]


def capture(shape, tile_size, overlap):
    """Reduce one config to a comparable value, or a marker if it raises."""
    from spida.utilities.tiling import tile_grid_with_overlap

    try:
        info_list = tile_grid_with_overlap(shape, tile_size, overlap)
    except Exception as exc:                       # noqa: BLE001
        return f"RAISED {type(exc).__name__}"
    return [
        (info["tile_id"], tuple(info["position"]), tuple(info["end_position"]),
         tuple(info["actual_size"]))
        for info in info_list
    ]


def build(keys: list | None = None) -> list:
    """Replay ``keys`` if given, else the fallback surface.

    Replaying the RECORDED keys is what makes a re-baseline auditable: the
    config set stays fixed, so every difference is a genuine behaviour change
    rather than an artefact of the parameter grid drifting.
    """
    if keys is None:
        keys = [(shape, repr(tile), repr(overlap))
                for shape, tile, overlap in itertools.product(SHAPES, TILES, OVERLAPS)]
    out = []
    for key in keys:
        shape, tile_repr, overlap_repr = key
        out.append((key, capture(shape, ast.literal_eval(tile_repr),
                                 ast.literal_eval(overlap_repr))))
    return out


def has_negative_origin(value) -> bool:
    return isinstance(value, list) and any(min(pos) < 0 for _, pos, _, _ in value)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accept", action="store_true",
                    help="write the new golden master (default: report only)")
    ap.add_argument("--allow-negative-origin", action="store_true",
                    help="permit configs with a negative tile origin (they are bugs)")
    args = ap.parse_args()

    old = pickle.loads(GOLDEN.read_bytes()) if GOLDEN.is_file() else []
    new = build([key for key, _ in old] or None)
    old_map = dict(old)

    changed, gained_negative = [], []
    for key, value in new:
        if key in old_map and old_map[key] != value:
            changed.append((key, old_map[key], value))
        if has_negative_origin(value):
            gained_negative.append(key)

    print(f"configs             : {len(new)} ({len(old)} recorded)")
    print(f"changed             : {len(changed)}")
    print(f"negative origins    : "
          f"{sum(1 for _, v in old if has_negative_origin(v))} before -> "
          f"{len(gained_negative)} now")
    for key, before, after in changed:
        def _fmt(v):
            return v if isinstance(v, str) else f"{len(v)} tiles"
        print(f"  {key[0]!s:>14} tile={key[1]:<10} overlap={key[2]:<9} "
              f"{_fmt(before)} -> {_fmt(after)}")

    if gained_negative and not args.allow_negative_origin:
        print(f"\nREFUSING: {len(gained_negative)} config(s) still produce a "
              f"negative tile origin, which silently shifts merged polygons and "
              f"breaks tile reconstruction. Fix the geometry, or pass "
              f"--allow-negative-origin if this is genuinely intended.",
              file=sys.stderr)
        return 1

    if not args.accept:
        print("\nreport only; pass --accept to overwrite the golden master")
        return 0

    GOLDEN.write_bytes(pickle.dumps(new))
    print(f"\nwrote {GOLDEN} ({len(new)} configs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
