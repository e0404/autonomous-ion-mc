"""CUDA backend tests (``cuda`` marker: skipped without a GPU; ``IONMC_REQUIRE_CUDA=1`` makes a
missing device a failure). They run on the GPU host runner."""

from __future__ import annotations

from collections.abc import Callable

import pytest

from ionmc.config import SimulationConfig
from ionmc.errors import UnsupportedCombinationError
from ionmc.simulation import Simulation
from ionmc.transport.parity import compare_partition, format_partition_verdict

pytestmark = pytest.mark.cuda
MakeConfig = Callable[..., SimulationConfig]


def test_cuda_smoke_float32(make_config: MakeConfig) -> None:
    cfg = make_config(
        energy=60.0, n=2000, n_batches=4, seed=3, backend="warp-cuda", precision="float32"
    )
    res = Simulation(cfg).run()
    assert res.valid and res.device.startswith("cuda")
    assert res.energy_balance.relative_residual < 1e-5
    part = res.transport_report["partials"][0]
    assert part["compile_s"] >= 0.0 and "register_count" in part
    assert max(part["chunk_seconds"]) < 60.0  # TDR guard: record, never silently pass long chunks


def test_cuda_chunk_sizes_partition_invariance(make_config: MakeConfig) -> None:
    """T13 (HR): chunk sizes 2**10 and 2**18 give identical tallies and counters and a deposit
    grid within 1e-5."""

    def run(chunk: int):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=60.0, n=20000, n_batches=4, seed=9, backend="warp-cuda",
            precision="float32", run_kwargs={"chunk_histories": chunk},
        )  # fmt: skip
        return Simulation(cfg).run()

    v = compare_partition(run(2**10), run(2**18))
    assert v["tallies_identical"] and v["counters_identical"] and v["pass"], (
        format_partition_verdict(v)
    )


def test_cuda_cpu_workers_rejected(make_config: MakeConfig) -> None:
    with pytest.raises(UnsupportedCombinationError):
        Simulation(
            make_config(
                n=4, backend="warp-cuda", precision="float32", run_kwargs={"cpu_workers": 2}
            )
        )


def test_cuda_and_cpu_float32_agree_statistically(make_config: MakeConfig) -> None:
    def run(backend: str):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=60.0, n=20000, n_batches=10, seed=4, backend=backend, precision="float32"
        )
        return Simulation(cfg).run()

    a, b = run("warp-cpu"), run("warp-cuda")
    assert a.energy_balance.initial_mev == b.energy_balance.initial_mev
    assert a.energy_balance.step_deposit_mev == pytest.approx(
        b.energy_balance.step_deposit_mev, rel=1e-3
    )
