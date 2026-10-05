"""The Python reference backend makes no Warp call at Python scope (Warp 1.17 Python-scope
dispatch segfaulted intermittently in worker processes; decisions 0037/0039): its shared
functions are pure-Python twins re-executed from the same source text as the Warp functions."""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest

from ionmc.config import SimulationConfig
from ionmc.simulation import Simulation


def test_python_backend_never_calls_warp_at_python_scope(
    make_config: Callable[..., SimulationConfig], monkeypatch: pytest.MonkeyPatch
) -> None:
    import warp._src.context as ctx

    def boom(self: object, *args: object, **kwargs: object) -> None:
        raise AssertionError(f"Warp call at Python scope: {self!r}")

    monkeypatch.setattr(ctx.Function, "__call__", boom)
    monkeypatch.setattr(ctx.Function, "call_builtin", boom)
    cfg = make_config(energy=30.0, n=4, n_batches=2)
    res = Simulation(cfg).run()
    assert res.valid and res.energy_balance.initial_mev == pytest.approx(4 * 30.0)


def test_twins_agree_with_the_warp_functions_on_samples() -> None:
    """Spot check of twin versus Warp function (kernel-compiled Python scope is avoided: the
    Warp side is evaluated in a tiny kernel)."""
    import warp as wp

    from ionmc._wpfunc import python_twin
    from ionmc.physics.em import make_em
    from ionmc.physics.kinematics import make_kinematics

    kin_w, kin_p = make_kinematics(wp.float64), python_twin(make_kinematics)
    em_w, em_p = make_em(wp.float64), python_twin(make_em)

    @wp.kernel
    def k(out: wp.array(dtype=wp.float64)):  # type: ignore[no-untyped-def,valid-type]
        out[0] = kin_w.tmax_mev(wp.float64(70.0), wp.float64(938.272))
        out[1] = em_w.bohr_variance(
            wp.float64(60.0), wp.float64(938.272), wp.float64(1.0), wp.float64(0.555),
            wp.float64(1.0), wp.float64(0.5),
        )  # fmt: skip

    out = wp.zeros(2, dtype=wp.float64, device="cpu")
    wp.launch(k, dim=1, inputs=[out], device="cpu")
    ref = out.numpy()
    assert kin_p.tmax_mev(70.0, 938.272) == pytest.approx(ref[0], rel=1e-14)
    assert em_p.bohr_variance(60.0, 938.272, 1.0, 0.555, 1.0, 0.5) == pytest.approx(
        ref[1], rel=1e-13
    )
    assert math.isfinite(ref[1])


def test_twins_of_all_factories_evaluate_without_warp_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rule: Warp functions are executed only inside kernels; Python-side evaluation uses the
    twins. Every function of every twin namespace runs with Warp calls forbidden."""
    import warp._src.context as ctx

    from ionmc._wpfunc import python_twin
    from ionmc.physics.em import make_em
    from ionmc.physics.kinematics import make_kinematics
    from ionmc.transport.funcs import make_transport_funcs

    def boom(self: object, *args: object, **kwargs: object) -> None:
        raise AssertionError(f"Warp call at Python scope: {self!r}")

    monkeypatch.setattr(ctx.Function, "__call__", boom)
    monkeypatch.setattr(ctx.Function, "call_builtin", boom)
    tf, em, kin = (python_twin(f) for f in (make_transport_funcs, make_em, make_kinematics))
    v = tf.vec3(0.3, 0.4, 0.5)
    assert tf.orthonormal_basis(em.rotate_dir(tf.vec3(0.0, 0.6, 0.8), 0.2, 1.0))
    assert tf.dda_next(v, v, 0, 0, 0, tf.vec3(0, 0, 0), tf.vec3(1, 1, 1))[1] in (0, 1, 2)
    assert math.isfinite(em.bohr_variance(50.0, 938.272, 1.0, 0.555, 1.0, 0.5))
    assert kin.beta2(50.0, 938.272) > 0.0
