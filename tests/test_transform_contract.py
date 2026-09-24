"""Every registered transform must satisfy the Transform protocol."""
import inspect

import pytest

from spida.S.image_prep.transforms import TRANSFORMS, TileContext, resolve_transform
from spida.S.image_prep.transforms.base import check_params


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_transform_declares_the_contract(name):
    tf = TRANSFORMS[name]
    assert tf.NAME == name
    assert isinstance(tf.USES_GPU, bool)
    assert isinstance(tf.DEFAULTS, dict)
    assert callable(tf.validate) and callable(tf.run)


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_run_takes_the_standard_positional_args(name):
    params = list(inspect.signature(TRANSFORMS[name].run).parameters)
    assert params[:3] == ["src", "dst", "ctx"]


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_run_takes_no_var_kwargs(name):
    """**kwargs would silently swallow a mistyped parameter, which is what the
    contract exists to prevent."""
    kinds = [p.kind for p in inspect.signature(TRANSFORMS[name].run).parameters.values()]
    assert inspect.Parameter.VAR_KEYWORD not in kinds


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_defaults_are_accepted_by_run(name):
    tf = TRANSFORMS[name]
    assert tf.validate(dict(tf.DEFAULTS)) is not None


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_validate_rejects_unknown_params(name):
    with pytest.raises(ValueError, match="does not accept"):
        TRANSFORMS[name].validate({"definitely_not_a_param": 1})


def test_resolve_transform_names_the_alternatives():
    with pytest.raises(ValueError, match="deconwolf"):
        resolve_transform("nope")


def test_tile_context_defaults_to_no_gpu():
    ctx = TileContext(channel="DAPI", colors=[405], z_step_um=1.5)
    assert ctx.gpu is False and ctx.gpu_lock is None


@pytest.mark.parametrize("name", sorted(TRANSFORMS))
def test_pool_dispatch_is_picklable(name):
    """The orchestrator hands fresh_run to mp.Pool, so every argument must pickle.
    Passing the transform module rather than its name failed here at runtime."""
    import multiprocessing as mp
    import pickle
    from functools import partial

    from spida.S.image_prep.orchestrator import fresh_run

    with mp.Manager() as manager:
        ctx = TileContext(channel="DAPI", colors=[405], z_step_um=1.5,
                          gpu=True, gpu_lock=manager.Lock())
        fn = partial(fresh_run, _transform_name=name, _ctx=ctx,
                     **TRANSFORMS[name].DEFAULTS)
        assert pickle.loads(pickle.dumps(fn)).keywords["_transform_name"] == name


def test_check_params_ignores_the_context_args():
    def run(src, dst, ctx, alpha=1):
        ...
    assert check_params(run, {"alpha": 2}, "x") == {"alpha": 2}
    with pytest.raises(ValueError):
        check_params(run, {"ctx": 2}, "x")
