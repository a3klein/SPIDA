"""
FILENAME: dark_sectioning.py
DESCRIPTION: 
Python port of "Dark sectioning" (Cao, Xi et al., Nature Methods 2025), 
official repo (github.com/Cao-ruijie/Dark-sectioning, Apache). 
Library module — no side effects on import.
The core dehazing math is Stephen Tierney's 2014 dark-channel dehazing.

NUMERICAL PROVENANCE. This is a verbatim copy of the bench port, which is validated
against the upstream MATLAB run under Octave to max|diff|/range = 2.98e-08,
r = 1.000000000000, across four optical regimes.

If used, dark sectioning must run before deconvolution.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import imageio.v3 as iio
import numpy as np
from scipy.ndimage import minimum_filter, uniform_filter
from scipy.special import j1

from spida.S.image_prep.transforms.base import TileContext, check_params

logger = logging.getLogger(__package__)

# `background` presets from dark_sectioning.m. 0 is the default single pass;
# 1 is the "severe" three-pass mode, which re-feeds its own output without
# re-normalising to 0-255 between passes (so `thres` means less each pass).
# NOTE on the shipped values: upstream annotates `deg` as "3-10" (6 is in range),
# `dep` as "0.7-2" (it ships 3 -- ABOVE the range) and `hl` as "3-8" (it ships 1 --
# BELOW the range), in both presets. `dep` scales the atmosphere, i.e. how much is
# subtracted; `hl` divides the processed low band in `Lo/hl + Hi`. Both are
# therefore magnitude-of-correction knobs shipped outside their own documented
# bounds. Overridable via dark_sectioning_stack(dep=..., hl=...).
BACKGROUND_PRESETS: dict[int, dict[str, tuple[float, ...]]] = {
    0: {"deg": (6.0,), "dep": (3.0,), "hl": (1.0,)},
    1: {"deg": (6.0, 3.0, 1.2), "dep": (3.0, 2.0, 2.0), "hl": (1.0, 1.0, 1.0)},
}

PAD_SIZE = 40
OMEGA = 0.95
GUIDED_RADIUS = 15
GUIDED_EPS = 0.001
DIVIDE = 0.5
THRES = 50.0

# `thres` splits "background" from "information" for the purpose of estimating
# min_atm -- the LOW end of the atmosphere. Upstream fixes it at 50 on the
# 0-255-normalised image, and that is not transferable: the normalisation is by
# stack min/max, so a FOV whose max is set by a few bright nuclei has almost all
# of its tissue pushed below 50. Measured on DAPI z=3 -- fraction of Lo under 50:
#   THM1 (sparse)  96.19%   Lo median  4.4
#   VePo (dense)   30.78%   Lo median 88.7
# The same number therefore selects "nearly the whole image" in one FOV and "under
# a third" in the other, so min_atm is estimated from tissue in the first case.
#
# THRES_MODES:
#   "abs"   the published behaviour, thres taken literally
#   "pct"   thres is a PERCENTILE of Lo. Scale-free, and it makes the estimator
#           see the same RANK range in every FOV, which is what transfers.
#   "otsu"  data-driven bimodal split. Note it does not stabilise coverage either
#           (96.7% vs 52.6% here) because the two FOVs genuinely differ in tissue
#           fraction -- it is included to test that, not because it is obviously right.
THRES_MODES = ("abs", "pct", "otsu")
# The mask must not be degenerate. An empty mask gives min_atm = 0 and an
# all-inclusive one gives min_atm = max_atm; both break the atmosphere estimate.
MASK_COVER_MIN, MASK_COVER_MAX = 0.01, 0.99


# ---------------------------------------------------------------------------
# Frequency-domain helpers
# ---------------------------------------------------------------------------

def lp_gauss(h: int, w: int, sigma: float) -> np.ndarray:
    """MATLAB `lpgauss(H, W, SIGMA)`: Gaussian lowpass in fft (unshifted) order."""
    kcx = float(sigma)
    kcy = (h / w) * float(sigma)
    x = np.arange(-(w // 2), (w - 1) // 2 + 1)
    y = np.arange(-(h // 2), (h - 1) // 2 + 1)
    xx, yy = np.meshgrid(x, y)
    return np.fft.ifftshift(np.exp(-(xx**2 / kcx**2 + yy**2 / kcy**2)))


def separate_hi_lo(img, emwavelength, na, pixel_size, factor, deg, divide=DIVIDE):
    """MATLAB `separateHiLo`. Returns (Hi, Lo, lp, EL)."""
    ny, nx = img.shape
    res = 0.5 * emwavelength / na / factor
    k_m = nx / (res / pixel_size)          # MATLAB Ny = columns; see fidelity note 2
    kc = np.floor(k_m * 0.2)
    sigma_lp = kc * 2 / 2.355

    lp = lp_gauss(ny, nx, sigma_lp * 2 * divide)
    hp = 1.0 - lp
    elp = lp_gauss(ny, nx, sigma_lp / deg)

    f = np.fft.fft2(img)
    hi = np.real(np.fft.ifft2(f * hp))
    lo = np.real(np.fft.ifft2(f * lp))
    el = np.real(np.fft.ifft2(f * elp))
    return hi, lo, lp, el


def psf_generator(emwavelength: float, pixel_size: float, na: float, w: int,
                  factor: float) -> np.ndarray:
    """MATLAB `PSF_Generator`: an Airy intensity PSF on a w-by-w grid."""
    x = np.arange(w)
    xx, yy = np.meshgrid(x, x)
    scale = 2 * np.pi * na / emwavelength * pixel_size * factor
    r = np.sqrt(np.minimum(xx, np.abs(xx - w)) ** 2
                + np.minimum(yy, np.abs(yy - w)) ** 2)
    eps = np.finfo(float).eps
    psf = np.abs(2 * j1(scale * r + eps) / (scale * r + eps)) ** 2
    return np.fft.fftshift(psf / psf.sum())


def confirm_block(emwavelength, na, pixel_size, factor, n, lp) -> int:
    """MATLAB `confirm_block`: dark-channel window = lowpassed-PSF 1% radius."""
    psf = psf_generator(emwavelength, pixel_size, na, n, factor)
    psf_lo = np.abs(np.fft.ifft2(np.fft.fftshift(np.fft.fft2(psf))
                                 * np.fft.fftshift(lp)))
    psf_lo = psf_lo / psf_lo.max()
    center = n // 2 - 1                     # MATLAB floor(Nx/2), 1-based
    count_x = n - 1
    for count_x in range(center, n):
        if psf_lo[count_x, center] < 0.01:
            break
    return count_x - center


# ---------------------------------------------------------------------------
# Dark-channel dehazing (Tierney 2014, as vendored by the paper)
# ---------------------------------------------------------------------------

def get_dark_channel(image: np.ndarray, win_size: int) -> np.ndarray:
    """Local minimum over a win_size box.

    MATLAB pads with +Inf by floor(win_size/2) and reads the window starting at
    each output pixel. In original coordinates that window is
    [i - floor(w/2), i + ceil(w/2) - 1] — exactly scipy's centred window for both
    odd and even `win_size` — so a centred `minimum_filter` with an infinite
    constant border is the same operation without the explicit pad.

    NaN is mapped to +Inf first. This is not a guard bolted on: MATLAB's `min`
    IGNORES NaN, numpy's propagates it, and the caller divides by an atmosphere
    that can legitimately reach 0 (giving 0/0). Substituting +Inf reproduces
    MATLAB exactly — a NaN can never win a min there either — whereas letting it
    propagate would smear NaN over a whole window and then across the image via
    the guided filter. The one case that still differs is an all-NaN window,
    where MATLAB returns NaN and this returns +Inf; both are degenerate.
    """
    return minimum_filter(np.nan_to_num(image, nan=np.inf), size=win_size,
                          mode="constant", cval=np.inf)


def get_atmosphere(image: np.ndarray, dark_channel: np.ndarray) -> float:
    """Mean image value over the top 1% of dark-channel pixels.

    Column-major ravel + stable descending sort to match MATLAB tie ordering
    (see fidelity note 3) — the dark channel is heavily tied, so this is not
    cosmetic.
    """
    n_search = int(np.floor(image.size * 0.01))
    if n_search < 1:
        n_search = 1
    dark_vec = dark_channel.ravel(order="F")
    img_vec = image.ravel(order="F")
    idx = np.argsort(-dark_vec, kind="stable")[:n_search]
    return float(img_vec[idx].mean())


def guided_filter(guide, target, radius=GUIDED_RADIUS, eps=GUIDED_EPS):
    """MATLAB `guided_filter` on `window_sum_filter` — truncated-window means."""
    win = 2 * radius + 1
    denom = uniform_filter(np.ones_like(guide, dtype=np.float64), win,
                           mode="constant", cval=0.0)

    def box_mean(a: np.ndarray) -> np.ndarray:
        return uniform_filter(a, win, mode="constant", cval=0.0) / denom

    mean_g = box_mean(guide)
    mean_t = box_mean(target)
    var_g = box_mean(guide * guide) - mean_g * mean_g
    cov_gt = box_mean(guide * target) - mean_g * mean_t

    a = cov_gt / (var_g + eps)
    b = mean_t - a * mean_g
    return box_mean(a) * guide + box_mean(b)


def resolve_thres(image: np.ndarray, thres: float, mode: str) -> float:
    """Turn (thres, mode) into an absolute level on this slice's Lo band.

    Clamped so the resulting background mask covers between MASK_COVER_MIN and
    MASK_COVER_MAX of the slice — outside that the min_atm estimate is degenerate
    rather than merely different.
    """
    if mode == "abs":
        level = float(thres)
    elif mode == "pct":
        level = float(np.percentile(image, thres))
    elif mode == "otsu":
        from skimage.filters import threshold_otsu
        level = float(threshold_otsu(image))
    else:
        raise ValueError(f"thres_mode must be one of {THRES_MODES}, got {mode!r}")

    cover = float((image < level).mean())
    if cover < MASK_COVER_MIN:
        level = float(np.percentile(image, MASK_COVER_MIN * 100))
    elif cover > MASK_COVER_MAX:
        level = float(np.percentile(image, MASK_COVER_MAX * 100))
    return level


def dehaze_fast2(image, omega, win_size, el, dep, thres, thres_mode="abs",
                 atmosphere_bounds=None):
    """MATLAB `dehaze_fast2`: spatially varying atmosphere from the very-low band.

    `atmosphere_bounds` overrides the derived `(min_atm, max_atm)`. `get_atmosphere`
    is global over the array it is given, so deriving it per tile subtracts a
    different background level in each one. None reproduces upstream behaviour.
    """
    level = resolve_thres(image, thres, thres_mode)
    mask = (image < level).astype(np.float64)

    if atmosphere_bounds is None:
        min_atm = get_atmosphere(image * mask, get_dark_channel(image * mask, win_size))
        max_atm = get_atmosphere(image, get_dark_channel(image, win_size))
    else:
        min_atm, max_atm = (float(v) for v in atmosphere_bounds)

    el = el - el.min()
    atmosphere = dep * (el / el.max() * (max_atm - min_atm) + min_atm)

    # No epsilon on the divisor: MATLAB divides straight through, and adding one
    # would turn 0/0 (NaN, ignored by MATLAB's min) into 0, which WOULD win the
    # min and change trans_est. get_dark_channel handles NaN/Inf the way MATLAB
    # does instead. On real data min_atm > 0 so this never triggers; it only
    # arises when the whole `image < thres` mask comes up empty.
    with np.errstate(divide="ignore", invalid="ignore"):
        trans_est = 1 - omega * get_dark_channel(image / atmosphere, win_size)
    transmission = guided_filter(image, trans_est)
    return (image - atmosphere) / np.maximum(transmission, 0.1) + atmosphere


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def _pad_square(img: np.ndarray) -> np.ndarray:
    """MATLAB pads the SHORT axis with zeros up to square before filtering."""
    ny, nx = img.shape
    if ny == nx:
        return img
    n = max(ny, nx)
    out = np.zeros((n, n), dtype=img.dtype)
    out[:ny, :nx] = img
    return out


def _mirror_pad(img: np.ndarray, pad_size: int = PAD_SIZE):
    ny, nx = img.shape
    py = int(np.floor(ny / pad_size)) + 1
    px = int(np.floor(nx / pad_size)) + 1
    return np.pad(img, ((py, py), (px, px)), mode="symmetric"), py, px


def _one_pass(image_sq, emwavelength, na, pixel_size, factor,
              deg, dep, hl, thres, divide, pad_size, thres_mode="abs",
              atmosphere_bounds=None):
    """One dark-sectioning pass over an already-square slice. Returns (result, block)."""
    padded, py, px = _mirror_pad(image_sq, pad_size)
    hi, lo, lp, el = separate_hi_lo(padded, emwavelength, na, pixel_size,
                                    factor, deg, divide)
    block = confirm_block(emwavelength, na, pixel_size, factor, padded.shape[0], lp)
    lo_proc = dehaze_fast2(lo, OMEGA, block, el, dep, thres, thres_mode,
                           atmosphere_bounds=atmosphere_bounds)
    result = lo_proc / hl + hi
    ny, nx = image_sq.shape
    return result[py:py + ny, px:px + nx], block


def dark_sectioning_stack(stack, emwavelength, na, pixel_size, factor,
                          background=0, thres=THRES, divide=DIVIDE,
                          pad_size=PAD_SIZE, workers=1, verbose=False,
                          dep=None, hl=None, deg=None, thres_mode="abs",
                          input_range=None, atmosphere_bounds=None,
                          output_max=None):
    """Dark sectioning on a (Z, Y, X) stack — or a single (Y, X) slice.

    Parameters mirror `dark_sectioning.m`: `emwavelength` in nm, `pixel_size` in
    nm, `na` unitless, `factor` the resolution factor (1 for widefield; the
    paper's own demo uses 2 for its super-resolved input).

    `background=0` is the published default single pass; `background=1` is the
    three-pass "severe background" preset.

    `workers` fans the slices of one pass across processes. Slices are
    independent within a pass and the only cross-slice coupling is the final
    whole-stack rescale, so this changes nothing numerically. Budget ~2.5 GB per
    worker at 4000x4000.

    Returns (result, meta). `result` is float32 scaled to 0-65535 by the max over
    the WHOLE stack (fidelity note 6) and may contain negatives, which is why it
    is not written as uint16. Shape matches the input.
    """
    single = stack.ndim == 2
    vol = stack[None] if single else stack
    if vol.ndim != 3:
        raise ValueError(f"expected (Z, Y, X) or (Y, X), got {stack.shape}")

    presets = BACKGROUND_PRESETS[background]
    # Per-pass overrides. Upstream ships dep=3 against its own annotated range of
    # 0.7-2 and hl=1 against 3-8 (see BACKGROUND_PRESETS), so these two are worth
    # being able to move without editing the presets. A scalar override is applied
    # to every pass.
    if dep is not None:
        presets = {**presets, "dep": (float(dep),) * len(presets["dep"])}
    if hl is not None:
        presets = {**presets, "hl": (float(hl),) * len(presets["hl"])}
    if deg is not None:
        presets = {**presets, "deg": (float(deg),) * len(presets["deg"])}
    ny0, nx0 = vol.shape[1:]

    # MATLAB normalises the whole stack to 0-255 once, before any passes.
    # `input_range` overrides the per-array (min, max): this is an affine remap, so
    # per-tile bounds give per-tile gain and a brightness step at every seam.
    vol = vol.astype(np.float64)
    if input_range is None:
        lo, hi = float(vol.min()), float(vol.max())
    else:
        lo, hi = (float(v) for v in input_range)
    work = np.stack([_pad_square(s) for s in 255.0 * (vol - lo) / (hi - lo + 1e-12)])

    blocks: list[int] = []
    for pass_i, (deg, dep, hl) in enumerate(
            zip(presets["deg"], presets["dep"], presets["hl"])):
        fn = partial(_one_pass, emwavelength=emwavelength, na=na,
                     pixel_size=pixel_size, factor=factor, deg=deg, dep=dep,
                     hl=hl, thres=thres, divide=divide, pad_size=pad_size,
                     thres_mode=thres_mode,
                     atmosphere_bounds=atmosphere_bounds)
        if workers > 1 and work.shape[0] > 1:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                results = list(ex.map(fn, list(work)))
        else:
            results = [fn(s) for s in work]

        out = np.empty_like(work)
        for jz, (slice_out, block) in enumerate(results):
            out[jz] = slice_out
            blocks.append(block)
            if verbose:
                print(f"  [dark] pass {pass_i} z{jz} block_size={block}", flush=True)
        work = out
        del results

    work = work[:, :ny0, :nx0]              # undo the square pad
    # `output_max`: the other half of the affine remap, same reasoning as input_range.
    _omax = float(output_max) if output_max is not None else float(work.max())
    result = 65535.0 * work / (_omax + 1e-12)
    meta = {
        "dep": presets["dep"][0],
        "hl": presets["hl"][0],
        "deg": presets["deg"][0],
        "emwavelength_nm": emwavelength,
        "na": na,
        "pixel_size_nm": pixel_size,
        "factor": factor,
        "background": background,
        "thres": thres,
        "thres_mode": thres_mode,
        "divide": divide,
        "pad_size": pad_size,
        "block_size": blocks[0] if blocks else None,
        "n_passes": len(presets["deg"]),
        "input_min": lo,
        "input_max": hi,
    }
    return (result[0] if single else result).astype(np.float32), meta


def derived_parameters(shape, emwavelength, na, pixel_size, factor,
                       deg=6.0, divide=DIVIDE, pad_size=PAD_SIZE) -> dict:
    """Report the frequency cut-offs and dark-channel window for a given geometry.

    Cheap (one PSF, no image), so it is worth printing before committing a run —
    `block_size` is the parameter that actually sets how much background is
    removed, and it is derived, not chosen.
    """
    ny, nx = shape[-2:]
    n = max(ny, nx)
    py = int(np.floor(n / pad_size)) + 1
    npad = n + 2 * py
    res = 0.5 * emwavelength / na / factor
    k_m = npad / (res / pixel_size)
    kc = np.floor(k_m * 0.2)
    sigma_lp = kc * 2 / 2.355
    lp = lp_gauss(npad, npad, sigma_lp * 2 * divide)
    return {
        "padded_size": npad,
        "resolution_nm": res,
        "k_m": k_m,
        "kc": float(kc),
        "sigma_lp": sigma_lp,
        "sigma_lp_used": sigma_lp * 2 * divide,
        "sigma_elp_used": sigma_lp / deg,
        "block_size": confirm_block(emwavelength, na, pixel_size, factor, npad, lp),
    }


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

NAME = "dark_sectioning"
USES_GPU = False

# `dep=1` rather than the upstream 3, which sits above its own annotated range
# (`% 0.7-2`) and hollows out nuclei. The rest are as shipped.
# See dev_notebooks/cellpose_3d_test/dark_staining/README.md
DS_DEFAULTS = {"dep": 1.0, "factor": 1.0, "hl": 1.0,
               "thres": 50.0, "thres_mode": "abs", "background": 0}
DEFAULTS = dict(DS_DEFAULTS)


def validate(params: dict) -> dict:
    return check_params(run, params, NAME)


def run(
    src: Path,
    dst: Path,
    ctx: TileContext,
    dep: float = 1.0,
    factor: float = 1.0,
    hl: float = 1.0,
    thres: float = 50.0,
    thres_mode: str = "abs",
    background: int = 0,
    emwavelength: float | None = None,
    na: float | None = None,
    pixel_size: float | None = None,
    workers: int = 1,
    input_range=None,
    atmosphere_bounds=None,
    output_max=None,
):
    """Dark-section one tile to `dst`, returning that path.

    Optics default to the `[Optics]` config section; `emwavelength` falls back to the
    channel. `input_range` / `atmosphere_bounds` / `output_max` take whole-mosaic
    values; None derives them per tile.
    """
    from spida.S.image_prep.config import get_optics

    img_path, decon_path = Path(src), Path(dst)
    channel = ctx.channel

    # `colors` carries either a wavelength (405/488, as in the PSF filenames) or a
    # channel name.
    if emwavelength is None and ctx.colors is not None:
        c = ctx.colors[0] if isinstance(ctx.colors, (list, np.ndarray)) else ctx.colors
        try:
            emwavelength = float(c)
        except (TypeError, ValueError):
            channel = str(c)

    optics = get_optics(channel)
    if emwavelength is None:
        emwavelength = optics.get("emission_nm")
    if emwavelength is None:
        raise ValueError(
            "dark_sectioning needs an emission wavelength: pass emwavelength=, or "
            "channel= matching an emission_nm_<channel> key in the [Optics] config"
        )
    na = float(na if na is not None else optics["na"])
    pixel_size = float(pixel_size if pixel_size is not None
                       else optics["pixel_size_nm"])

    params = {"dep": dep, "factor": factor, "hl": hl,
              "thres": thres, "thres_mode": thres_mode, "background": background}
    # Whole-mosaic overrides; None derives them per tile. See dark_sectioning_stack.
    globals_ = {"input_range": input_range,
                "atmosphere_bounds": atmosphere_bounds,
                "output_max": output_max}

    stack = iio.imread(img_path)
    if stack.ndim == 4:
        raise ValueError(
            f"dark_sectioning expects (Z, Y, X) or (Y, X), got {stack.shape}; "
            "the tile loop runs one channel at a time"
        )

    logger.info(
        "dark_sectioning %s: shape=%s emwavelength=%s na=%s pixel_size=%s %s",
        img_path.name, stack.shape, emwavelength, na, pixel_size, params,
    )
    result, meta = dark_sectioning_stack(
        stack, emwavelength=emwavelength, na=na, pixel_size=pixel_size,
        workers=workers, **params, **globals_,
    )

    # `result` may be negative. Shift rather than clip: clipping makes an exactly-zero
    # population, and zero is absorbing under Richardson-Lucy's multiplicative update.
    lo = float(np.nanmin(result))
    if lo < 0:
        result = result - lo
        logger.info("dark_sectioning %s: shifted floor by %.4g (min was negative)",
                    img_path.name, -lo)
    mx = float(np.nanmax(result))
    out = (65535.0 * result / mx) if mx > 0 else result
    iio.imwrite(decon_path, out.astype(np.uint16))
    logger.info("dark_sectioning %s -> %s (block_size=%s)",
                img_path.name, Path(decon_path).name, meta.get("block_size"))
    return decon_path
