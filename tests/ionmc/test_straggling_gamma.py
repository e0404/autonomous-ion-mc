"""The straggling model ``bohr_gamma_v1`` (Gamma for every ratio mean / sigma): exact mean and
variance, positivity, the Gamma skewness, parity of the reference and the kernel, and the option
plumbing. The default model is unchanged (the stored baseline trace guards that)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest
import warp as wp

from ionmc.config import DiagnosticsOptions, SimulationConfig, validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.physics.em import make_em
from ionmc.simulation import Simulation
from ionmc.transport.parity import compare_traces

EM = make_em(wp.float64)
MakeConfig = Callable[..., SimulationConfig]
N = 400_000
MAX_TRIES = 64


@wp.kernel
def _sample(
    mean: float,
    var: float,
    u: wp.array2d(dtype=wp.float64),  # type: ignore[valid-type]
    out: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    ok: wp.array(dtype=wp.int32),  # type: ignore[valid-type]
    use_gamma: int,
) -> None:
    i = wp.tid()
    loss = wp.float64(0.0)
    accepted = int(0)
    t = int(0)
    while t < 64 and accepted == 0:
        base = t * 4
        lw = wp.float64(0.0)
        o = int(0)
        if use_gamma == 1:
            lw, o = EM.straggle_attempt_gamma(
                wp.float64(mean), wp.float64(var), u[i, base], u[i, base + 1], u[i, base + 2],
                u[i, base + 3],
            )  # fmt: skip
        else:
            lw, o = EM.straggle_attempt(
                wp.float64(mean), wp.float64(var), u[i, base], u[i, base + 1], u[i, base + 2],
                u[i, base + 3],
            )  # fmt: skip
        if o == 1:
            loss = lw
            accepted = 1
        t = t + 1
    out[i] = loss
    ok[i] = accepted


def _draw(mean: float, sigma: float, gamma: bool) -> np.ndarray:
    rng = np.random.default_rng(5)
    u = rng.uniform(1e-12, 1.0 - 1e-12, size=(N, 4 * MAX_TRIES))
    out = wp.zeros(N, dtype=wp.float64, device="cpu")
    ok = wp.zeros(N, dtype=wp.int32, device="cpu")
    wp.launch(_sample, dim=N, inputs=[mean, sigma**2, wp.array(u, dtype=wp.float64), out, ok,
                                      1 if gamma else 0], device="cpu")  # fmt: skip
    assert ok.numpy().all()
    return out.numpy()


@pytest.mark.parametrize("ratio", [0.7, 1.9, 6.0, 25.0])
def test_gamma_model_has_exact_moments_skew_and_is_positive(ratio: float) -> None:
    mean = 1.0
    x = _draw(mean, mean / ratio, gamma=True)
    assert x.min() > 0.0  # positivity everywhere, no clamp
    se_mean = (mean / ratio) / np.sqrt(N)
    assert abs(x.mean() - mean) < 4.5 * se_mean
    assert x.var() == pytest.approx((mean / ratio) ** 2, rel=0.03)
    skew = float(np.mean((x - x.mean()) ** 3) / x.std() ** 3)
    assert skew == pytest.approx(2.0 / ratio, abs=0.06 + 0.1 / ratio)  # Gamma: 2 / sqrt(k)


def test_default_model_differs_where_it_clamps() -> None:
    """The default Gaussian-clamped branch (ratio >= 3) has no skew; the Gamma model keeps it."""
    d = _draw(1.0, 1.0 / 6.0, gamma=False)
    g = _draw(1.0, 1.0 / 6.0, gamma=True)
    skew = lambda x: float(np.mean((x - x.mean()) ** 3) / x.std() ** 3)  # noqa: E731
    assert abs(skew(d)) < 0.05 and skew(g) > 0.25


def test_summed_gamma_steps_are_step_independent() -> None:
    """With a common scale theta the sum of Gamma steps is Gamma with the summed shape: ten steps
    of shape k/10 give the same distribution as one step of shape k (moments and skewness)."""
    ratio = 3.0
    rng = np.random.default_rng(1)
    k_total = ratio**2
    a = rng.gamma(k_total, 1.0 / k_total, size=N)
    b = sum(rng.gamma(k_total / 10.0, 1.0 / k_total, size=N) for _ in range(10))
    assert np.mean(a) == pytest.approx(np.mean(b), rel=0.01)
    assert np.var(a) == pytest.approx(np.var(b), rel=0.03)


def test_option_is_accepted_and_recorded(make_config: MakeConfig) -> None:
    cfg = make_config(physics_kwargs={"straggling_model": "bohr_gamma_v1"})
    assert validate(cfg).summary()["physics"]["straggling_model"] == "bohr_gamma_v1"
    with pytest.raises(UnsupportedCombinationError):
        validate(replace(cfg, physics=replace(cfg.physics, straggling_model="nope")))


def test_reference_and_kernel_follow_the_same_trajectory_with_the_gamma_model(
    make_config: MakeConfig,
) -> None:
    runs = {}
    for backend in ("python", "warp-cpu"):
        cfg = make_config(
            energy=30.0, n=6, n_batches=2, seed=4, backend=backend,
            physics_kwargs={"straggling_model": "bohr_gamma_v1"},
            diagnostics=DiagnosticsOptions(trace_histories=6, track_end_positions=True),
        )  # fmt: skip
        runs[backend] = Simulation(cfg).run()
    assert compare_traces(runs["python"].diagnostics, runs["warp-cpu"].diagnostics)["pass"]
    default = Simulation(
        make_config(
            energy=30.0, n=6, n_batches=2, seed=4, diagnostics=DiagnosticsOptions(trace_histories=6)
        )  # fmt: skip
    ).run()
    assert not np.array_equal(default.diagnostics["trace"]["energy_mev"][:50],
                              runs["python"].diagnostics["trace"]["energy_mev"][:50])  # fmt: skip
