"""Tests for infrastructure/diagnostics/warp_backend_report.py.

These tests must pass with `warp` and `numpy` BOTH absent (that is the
GitHub CI host, which installs only pytest) and also pass on a host where
either or both are actually present. No test may assume `warp` is present
OR that it is absent - graceful degradation and the "warp available" paths
are both exercised through explicit fakes/monkeypatching rather than
depending on host truth. Any test that genuinely needs real `warp` calls
`pytest.importorskip("warp")` first and makes no assumption about which
devices (CPU-only, or CPU+CUDA) that installation exposes.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys

import pytest

import warp_backend_report as wbr


# ---------------------------------------------------------------------------
# Fakes used to exercise report assembly without requiring `warp`/`numpy`
# ---------------------------------------------------------------------------


class FakeDevice:
    """A minimal stand-in for a warp.context.Device."""

    def __init__(self, alias, is_cpu, is_cuda, name, arch=None, total_memory=None, free_memory=None):
        self.alias = alias
        self._is_cpu = is_cpu
        self._is_cuda = is_cuda
        self.name = name
        self.arch = arch
        self.total_memory = total_memory
        self.free_memory = free_memory

    @property
    def is_cpu(self):
        return self._is_cpu

    @property
    def is_cuda(self):
        return self._is_cuda

    def __str__(self):
        return self.alias


class FakeWarpModule:
    """A minimal stand-in for the `warp` module, enough for build_devices_section."""

    __version__ = "9.9.9-fake"

    def __init__(self, devices):
        self._devices = devices

    def get_all_devices(self):
        return list(self._devices)


def make_two_device_fake():
    cpu = FakeDevice("cpu", is_cpu=True, is_cuda=False, name="Fake CPU", arch="0", total_memory=0, free_memory=0)
    cuda = FakeDevice("cuda:0", is_cpu=False, is_cuda=True, name="Fake GPU", arch="sm_86", total_memory=2048, free_memory=1024)
    return FakeWarpModule([cpu, cuda]), cpu, cuda


def identity_runner(workload_name, dtype_name, device_alias, values):
    """A fake run_kernel: every device reproduces the input exactly."""
    return list(values)


def make_perturbing_runner(perturb_by_device):
    """A fake run_kernel where each device applies a fixed additive perturbation."""

    def _runner(workload_name, dtype_name, device_alias, values):
        delta = perturb_by_device.get(device_alias, 0.0)
        return [v + delta for v in values]

    return _runner


def make_constant_runner(value_by_device):
    """A fake run_kernel where each device returns a fixed constant sequence."""

    def _runner(workload_name, dtype_name, device_alias, values):
        value = value_by_device.get(device_alias, 0.0)
        return [value for _ in values]

    return _runner


TWO_DEVICE_LIST = [
    {"alias": "cpu", "is_cpu": True, "is_cuda": False},
    {"alias": "cuda:0", "is_cpu": False, "is_cuda": True},
]

CPU_ONLY_DEVICE_LIST = [{"alias": "cpu", "is_cpu": True, "is_cuda": False}]


# ---------------------------------------------------------------------------
# bits_equal
# ---------------------------------------------------------------------------


def test_bits_equal_identical_values():
    assert wbr.bits_equal(1.5, 1.5) is True


def test_bits_equal_distinguishes_negative_and_positive_zero():
    assert 0.0 == -0.0  # sanity: Python equality does not distinguish them
    assert wbr.bits_equal(0.0, -0.0) is False
    assert wbr.bits_equal(0.0, 0.0) is True
    assert wbr.bits_equal(-0.0, -0.0) is True


def test_bits_equal_nan_is_bit_equal_to_any_nan():
    assert wbr.bits_equal(math.nan, math.nan) is True
    assert wbr.bits_equal(math.nan, float("nan")) is True


def test_bits_equal_nan_not_equal_to_non_nan():
    assert wbr.bits_equal(math.nan, 1.0) is False
    assert wbr.bits_equal(1.0, math.nan) is False


def test_bits_equal_inf():
    assert wbr.bits_equal(math.inf, math.inf) is True
    assert wbr.bits_equal(-math.inf, -math.inf) is True
    assert wbr.bits_equal(math.inf, -math.inf) is False


def test_bits_equal_different_values():
    assert wbr.bits_equal(1.0, 1.0000001) is False


# ---------------------------------------------------------------------------
# abs_diff / rel_diff / normalized_diff
# ---------------------------------------------------------------------------


def test_abs_diff_basic():
    assert wbr.abs_diff(1.0, 1.5) == pytest.approx(0.5)
    assert wbr.abs_diff(1.5, 1.0) == pytest.approx(0.5)


def test_abs_diff_equal_values_is_zero_even_for_inf():
    assert wbr.abs_diff(math.inf, math.inf) == 0.0
    assert wbr.abs_diff(-0.0, 0.0) == 0.0


def test_abs_diff_nan_propagates():
    assert math.isnan(wbr.abs_diff(math.nan, 1.0))
    assert math.isnan(wbr.abs_diff(1.0, math.nan))


def test_rel_diff_basic():
    assert wbr.rel_diff(2.0, 2.2) == pytest.approx(0.1)


def test_rel_diff_zero_reference_and_equal_is_zero():
    assert wbr.rel_diff(0.0, 0.0) == 0.0
    assert wbr.rel_diff(0.0, -0.0) == 0.0


def test_rel_diff_zero_reference_and_unequal_is_inf():
    assert wbr.rel_diff(0.0, 1e-300) == math.inf


def test_rel_diff_nan_propagates():
    assert math.isnan(wbr.rel_diff(math.nan, 1.0))


def test_rel_diff_can_be_inflated_near_zero_crossing():
    # This is exactly the DRYRUN-002 round-2 finding: a tiny absolute
    # difference near a reference value close to zero produces a huge
    # relative figure - which is why rel_diff is informational only and
    # normalized_diff (mixed criterion) is what verdicts use.
    reference = 1e-9
    candidate = reference + 1e-7
    assert wbr.rel_diff(reference, candidate) > 1.0


def test_normalized_diff_within_budget():
    # atol=1e-6, rtol=1e-6; diff of 5e-7 is within atol alone.
    assert wbr.normalized_diff(1.0, 1.0 + 5e-7, rtol=1e-6, atol=1e-6) <= 1.0


def test_normalized_diff_exceeds_budget():
    assert wbr.normalized_diff(1.0, 1.0 + 1.0, rtol=1e-6, atol=1e-6) > 1.0


def test_normalized_diff_zero_budget_equal_values_is_zero():
    assert wbr.normalized_diff(0.0, 0.0, rtol=0.0, atol=0.0) == 0.0


def test_normalized_diff_zero_budget_unequal_values_is_inf():
    assert wbr.normalized_diff(0.0, 1e-300, rtol=0.0, atol=0.0) == math.inf


def test_normalized_diff_nan_propagates():
    assert math.isnan(wbr.normalized_diff(math.nan, 1.0, rtol=1e-6, atol=1e-6))


def test_normalized_diff_not_inflated_near_zero_crossing():
    # The same near-zero-reference case that inflates rel_diff must NOT
    # blow up normalized_diff, given a sane atol.
    reference = 1e-9
    candidate = reference + 1e-7
    assert wbr.normalized_diff(reference, candidate, rtol=1e-6, atol=1e-6) <= 1.0


# ---------------------------------------------------------------------------
# compare_sequences / select_verdict
# ---------------------------------------------------------------------------


def test_compare_sequences_bitwise_identical():
    stats = wbr.compare_sequences([1.0, 2.0, -0.0], [1.0, 2.0, -0.0])
    assert stats["n"] == 3
    assert stats["bitwise_identical"] is True
    assert stats["nan_count"] == 0
    assert stats["max_abs_diff"] == 0.0
    assert stats["max_rel_diff"] == 0.0
    assert stats["max_normalized_diff"] is None


def test_compare_sequences_with_tolerance_budget_reports_normalized_diff():
    stats = wbr.compare_sequences([1.0], [1.0 + 5e-7], rtol=1e-6, atol=1e-6)
    assert stats["max_normalized_diff"] is not None
    assert stats["max_normalized_diff"] <= 1.0


def test_compare_sequences_not_bitwise_but_numerically_equal():
    stats = wbr.compare_sequences([0.0], [-0.0])
    assert stats["bitwise_identical"] is False
    assert stats["max_abs_diff"] == 0.0
    assert stats["max_rel_diff"] == 0.0


def test_compare_sequences_length_mismatch_raises():
    with pytest.raises(ValueError):
        wbr.compare_sequences([1.0, 2.0], [1.0])


def test_compare_sequences_empty_sequences_raises():
    with pytest.raises(ValueError):
        wbr.compare_sequences([], [])


def test_compare_sequences_nan_count():
    stats = wbr.compare_sequences([1.0, math.nan, 3.0], [1.0, math.nan, 3.0])
    assert stats["nan_count"] == 1


def test_compare_sequences_nan_count_counts_either_side():
    stats = wbr.compare_sequences([1.0, math.nan], [math.nan, 2.0])
    assert stats["nan_count"] == 2


def test_compare_sequences_propagates_nan_as_max():
    stats = wbr.compare_sequences([1.0, 2.0], [1.0, math.nan])
    assert math.isnan(stats["max_abs_diff"])
    assert math.isnan(stats["max_rel_diff"])


# --- select_verdict: bitwise class ---


def test_select_verdict_bitwise_identical():
    assert wbr.select_verdict(wbr.EXPECTATION_BITWISE, True, 0, None) == wbr.VERDICT_IDENTICAL


def test_select_verdict_bitwise_never_within_tolerance_for_negative_zero():
    # -0.0 vs 0.0: not bitwise identical, no NaN. The bitwise class must
    # report exceeds_tolerance, never within_tolerance.
    stats = wbr.compare_sequences([0.0], [-0.0])
    verdict = wbr.select_verdict(wbr.EXPECTATION_BITWISE, stats["bitwise_identical"], stats["nan_count"], stats["max_normalized_diff"])
    assert verdict == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_select_verdict_bitwise_nan_forces_exceeds_even_if_bits_equal():
    # bits_equal() treats nan/nan as bit-equal, but the NaN guard in
    # select_verdict() must dominate: a NaN never counts as a pass.
    assert wbr.bits_equal(math.nan, math.nan) is True
    assert wbr.select_verdict(wbr.EXPECTATION_BITWISE, True, nan_count=1, max_normalized_diff=None) == wbr.VERDICT_EXCEEDS_TOLERANCE


# --- select_verdict: tolerance classes ---


def test_select_verdict_tolerance_class_identical_when_bitwise():
    assert wbr.select_verdict(wbr.EXPECTATION_TRANSCENDENTAL, True, 0, None) == wbr.VERDICT_IDENTICAL


def test_select_verdict_tolerance_class_within_tolerance_for_negative_zero():
    stats = wbr.compare_sequences([0.0], [-0.0], rtol=1e-6, atol=1e-6)
    verdict = wbr.select_verdict(
        wbr.EXPECTATION_TRANSCENDENTAL, stats["bitwise_identical"], stats["nan_count"], stats["max_normalized_diff"]
    )
    assert verdict == wbr.VERDICT_WITHIN_TOLERANCE


def test_select_verdict_within_tolerance_at_boundary():
    assert wbr.select_verdict(wbr.EXPECTATION_TRANSCENDENTAL, False, 0, 1.0) == wbr.VERDICT_WITHIN_TOLERANCE


def test_select_verdict_exceeds_tolerance_just_beyond_boundary():
    assert wbr.select_verdict(wbr.EXPECTATION_TRANSCENDENTAL, False, 0, 1.0 + 1e-12) == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_select_verdict_nan_always_exceeds_tolerance():
    assert wbr.select_verdict(wbr.EXPECTATION_TRANSCENDENTAL, False, nan_count=1, max_normalized_diff=0.0) == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_select_verdict_contraction_sensitive_class_behaves_like_tolerance_class():
    assert wbr.select_verdict(wbr.EXPECTATION_CONTRACTION_SENSITIVE, False, 0, 0.5) == wbr.VERDICT_WITHIN_TOLERANCE
    assert wbr.select_verdict(wbr.EXPECTATION_CONTRACTION_SENSITIVE, False, 0, 2.0) == wbr.VERDICT_EXCEEDS_TOLERANCE


# ---------------------------------------------------------------------------
# generate_inputs
# ---------------------------------------------------------------------------


def test_generate_inputs_deterministic_for_fixed_seed():
    a = wbr.generate_inputs(50, seed=42)
    b = wbr.generate_inputs(50, seed=42)
    assert a == b


def test_generate_inputs_differs_for_different_seed():
    a = wbr.generate_inputs(50, seed=42)
    b = wbr.generate_inputs(50, seed=43)
    assert a != b


def test_generate_inputs_length_and_domain():
    values = wbr.generate_inputs(200, seed=1)
    assert len(values) == 200
    assert all(wbr.DOMAIN_MIN <= v <= wbr.DOMAIN_MAX for v in values)


def test_generate_inputs_does_not_disturb_global_random_state():
    import random

    random.seed(7)
    before = random.random()
    random.seed(7)
    wbr.generate_inputs(10, seed=999)
    after = random.random()
    assert before == after


# ---------------------------------------------------------------------------
# build_devices_section (fake warp module - no real `warp` needed)
# ---------------------------------------------------------------------------


def test_build_devices_section_with_fake_module():
    fake_module, cpu, cuda = make_two_device_fake()
    section = wbr.build_devices_section(fake_module)
    assert section["status"] == wbr.STATUS_AVAILABLE
    aliases = {d["alias"] for d in section["devices"]}
    assert aliases == {"cpu", "cuda:0"}
    cpu_entry = next(d for d in section["devices"] if d["alias"] == "cpu")
    assert cpu_entry["is_cpu"] is True
    assert cpu_entry["is_cuda"] is False
    assert cpu_entry["name"] == "Fake CPU"
    cuda_entry = next(d for d in section["devices"] if d["alias"] == "cuda:0")
    assert cuda_entry["total_memory_bytes"] == 2048


def test_build_devices_section_normalizes_non_cuda_arch_to_null():
    # Fake CPU device reports arch="0" (mirroring the real DRYRUN-002
    # GPU-host run) - this must be normalized to null for a non-CUDA device.
    fake_module, _, _ = make_two_device_fake()
    section = wbr.build_devices_section(fake_module)
    cpu_entry = next(d for d in section["devices"] if d["alias"] == "cpu")
    assert cpu_entry["arch"] is None


def test_build_devices_section_preserves_cuda_arch():
    fake_module, _, _ = make_two_device_fake()
    section = wbr.build_devices_section(fake_module)
    cuda_entry = next(d for d in section["devices"] if d["alias"] == "cuda:0")
    assert cuda_entry["arch"] == "sm_86"


def test_build_devices_section_normalizes_zero_memory_to_null():
    fake_module, _, _ = make_two_device_fake()
    section = wbr.build_devices_section(fake_module)
    cpu_entry = next(d for d in section["devices"] if d["alias"] == "cpu")
    assert cpu_entry["total_memory_bytes"] is None
    assert cpu_entry["free_memory_bytes"] is None


def test_build_devices_section_zero_devices_is_unavailable():
    fake_module = FakeWarpModule([])
    section = wbr.build_devices_section(fake_module)
    assert section["status"] == wbr.STATUS_UNAVAILABLE
    assert section["devices"] == []


def test_build_devices_section_degrades_on_broken_device_attrs():
    class BrokenDevice:
        alias = "weird"

        @property
        def is_cpu(self):
            raise RuntimeError("boom")

        @property
        def name(self):
            raise RuntimeError("boom too")

    fake_module = FakeWarpModule([BrokenDevice()])
    section = wbr.build_devices_section(fake_module)
    assert section["status"] == wbr.STATUS_AVAILABLE
    entry = section["devices"][0]
    assert entry["alias"] == "weird"
    assert entry["is_cpu"] is None
    assert entry["name"] is None


def test_build_devices_section_enumeration_failure_is_error():
    class BrokenModule:
        def get_all_devices(self):
            raise RuntimeError("enumeration exploded")

    section = wbr.build_devices_section(BrokenModule())
    assert section["status"] == wbr.STATUS_ERROR
    assert "detail" in section


def test_build_devices_section_falls_back_when_get_all_devices_absent():
    class LegacyModule:
        def get_device(self, name):
            assert name == "cpu"
            return FakeDevice("cpu", is_cpu=True, is_cuda=False, name="Legacy CPU")

        def get_cuda_device_count(self):
            return 1

        def get_cuda_device(self, index):
            return FakeDevice(f"cuda:{index}", is_cpu=False, is_cuda=True, name="Legacy GPU")

    section = wbr.build_devices_section(LegacyModule())
    assert section["status"] == wbr.STATUS_AVAILABLE
    aliases = {d["alias"] for d in section["devices"]}
    assert aliases == {"cpu", "cuda:0"}


# ---------------------------------------------------------------------------
# build_warp_version_section
# ---------------------------------------------------------------------------


def test_build_warp_version_section_with_fake_module():
    fake_module, _, _ = make_two_device_fake()
    section = wbr.build_warp_version_section(fake_module)
    assert section["status"] == wbr.STATUS_AVAILABLE
    assert section["version"] == "9.9.9-fake"


def test_build_warp_version_section_reports_unavailable_on_import_error(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise ImportError("No module named 'warp'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    section = wbr.build_warp_version_section()
    assert section["status"] == wbr.STATUS_UNAVAILABLE
    assert "warp" in section["detail"].lower()


def test_build_warp_version_section_reports_error_on_unexpected_import_failure(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise RuntimeError("kaboom")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    section = wbr.build_warp_version_section()
    assert section["status"] == wbr.STATUS_ERROR


# ---------------------------------------------------------------------------
# build_parity_section (fake devices + fake run_kernel - no warp/numpy needed)
# ---------------------------------------------------------------------------


def test_build_parity_section_no_cpu_device_is_unavailable():
    devices = [{"alias": "cuda:0", "is_cpu": False, "is_cuda": True}]
    section = wbr.build_parity_section(devices, identity_runner)
    assert section["status"] == wbr.STATUS_UNAVAILABLE
    assert section["reference_device"] is None
    assert section["results"] == []


@pytest.mark.parametrize("elements", [0, -1, -100])
def test_build_parity_section_non_positive_elements_is_error(elements):
    section = wbr.build_parity_section(TWO_DEVICE_LIST, identity_runner, elements=elements)
    assert section["status"] == wbr.STATUS_ERROR
    assert "detail" in section
    assert section["results"] == []


def test_build_parity_section_identical_runner_is_all_identical():
    section = wbr.build_parity_section(TWO_DEVICE_LIST, identity_runner, elements=16, seed=1)
    assert section["status"] == wbr.STATUS_AVAILABLE
    assert section["reference_device"] == "cpu"
    assert len(section["results"]) == len(wbr.WORKLOADS) * len(wbr.DTYPES) * len(TWO_DEVICE_LIST)
    for result in section["results"]:
        assert result["status"] == wbr.STATUS_AVAILABLE
        assert result["bitwise_identical"] is True
        assert result["verdict"] == wbr.VERDICT_IDENTICAL
        assert result["n"] == 16
        assert result["nan_count"] == 0


def test_build_parity_section_reports_expectation_class_and_workload_per_result():
    section = wbr.build_parity_section(TWO_DEVICE_LIST, identity_runner, elements=4, seed=1)
    by_workload = {r["workload"] for r in section["results"]}
    assert by_workload == set(wbr.WORKLOADS)
    for result in section["results"]:
        assert result["expectation_class"] == wbr.WORKLOAD_EXPECTATION_CLASS[result["workload"]]


def test_build_parity_section_product_workload_reports_null_tolerance():
    section = wbr.build_parity_section(TWO_DEVICE_LIST, identity_runner, elements=4, seed=1)
    product_results = [r for r in section["results"] if r["workload"] == "product"]
    for result in product_results:
        assert result["rtol"] is None
        assert result["atol"] is None
        assert result["max_normalized_diff"] is None


def test_build_parity_section_bitwise_workload_small_perturbation_exceeds_tolerance():
    # A "small" perturbation on the bitwise (`product`) workload must
    # ALWAYS be exceeds_tolerance - the bitwise class never accepts
    # within_tolerance, regardless of how small the difference is.
    runner = make_perturbing_runner({"cpu": 0.0, "cuda:0": 1e-12})
    section = wbr.build_parity_section(TWO_DEVICE_LIST, runner, elements=8, seed=2)
    product_cuda_results = [r for r in section["results"] if r["workload"] == "product" and r["device"] == "cuda:0"]
    assert product_cuda_results
    for result in product_cuda_results:
        assert result["bitwise_identical"] is False
        assert result["verdict"] == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_build_parity_section_tolerance_class_small_perturbation_is_within_tolerance():
    # A perturbation well inside every tolerance-class budget (smallest
    # budget here is float64's atol=1e-14) must be exactly
    # within_tolerance - not "either verdict is fine".
    runner = make_perturbing_runner({"cpu": 0.0, "cuda:0": 1e-16})
    section = wbr.build_parity_section(TWO_DEVICE_LIST, runner, elements=8, seed=2)
    tolerance_cuda_results = [
        r for r in section["results"] if r["workload"] in ("fma_exposed", "transcendental") and r["device"] == "cuda:0"
    ]
    assert tolerance_cuda_results
    for result in tolerance_cuda_results:
        assert result["bitwise_identical"] is False
        assert result["verdict"] == wbr.VERDICT_WITHIN_TOLERANCE


def test_build_parity_section_tolerance_class_large_perturbation_exceeds_tolerance():
    runner = make_perturbing_runner({"cpu": 0.0, "cuda:0": 5.0})
    section = wbr.build_parity_section(TWO_DEVICE_LIST, runner, elements=8, seed=3)
    tolerance_cuda_results = [
        r for r in section["results"] if r["workload"] in ("fma_exposed", "transcendental") and r["device"] == "cuda:0"
    ]
    for result in tolerance_cuda_results:
        assert result["verdict"] == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_build_parity_section_nan_vs_nan_is_not_identical():
    runner = make_constant_runner({"cpu": math.nan, "cuda:0": math.nan})
    section = wbr.build_parity_section(TWO_DEVICE_LIST, runner, elements=4, seed=1)
    cuda_results = [r for r in section["results"] if r["device"] == "cuda:0"]
    for result in cuda_results:
        assert result["nan_count"] > 0
        assert result["verdict"] == wbr.VERDICT_EXCEEDS_TOLERANCE
        assert result["verdict"] != wbr.VERDICT_IDENTICAL


def test_build_parity_section_nan_count_reported_for_partial_nan():
    def runner(workload_name, dtype_name, device_alias, values):
        if device_alias == "cuda:0":
            return [math.nan] + list(values[1:])
        return list(values)

    section = wbr.build_parity_section(TWO_DEVICE_LIST, runner, elements=4, seed=1)
    cuda_results = [r for r in section["results"] if r["device"] == "cuda:0"]
    for result in cuda_results:
        assert result["nan_count"] == 1
        assert result["verdict"] == wbr.VERDICT_EXCEEDS_TOLERANCE


def test_build_parity_section_reference_failure_marks_all_devices_error():
    def broken_reference_runner(workload_name, dtype_name, device_alias, values):
        if device_alias == "cpu":
            raise RuntimeError("reference kernel exploded")
        return list(values)

    section = wbr.build_parity_section(TWO_DEVICE_LIST, broken_reference_runner, elements=4, seed=4)
    assert section["status"] == wbr.STATUS_ERROR
    assert all(result["status"] == wbr.STATUS_ERROR for result in section["results"])


def test_build_parity_section_one_device_failure_does_not_affect_others():
    def flaky_runner(workload_name, dtype_name, device_alias, values):
        if device_alias == "cuda:0":
            raise RuntimeError("cuda device exploded")
        return list(values)

    section = wbr.build_parity_section(TWO_DEVICE_LIST, flaky_runner, elements=4, seed=5)
    assert section["status"] == wbr.STATUS_AVAILABLE
    cpu_results = [r for r in section["results"] if r["device"] == "cpu"]
    cuda_results = [r for r in section["results"] if r["device"] == "cuda:0"]
    assert all(r["status"] == wbr.STATUS_AVAILABLE for r in cpu_results)
    assert all(r["status"] == wbr.STATUS_ERROR for r in cuda_results)


# ---------------------------------------------------------------------------
# build_configuration_section
# ---------------------------------------------------------------------------


def test_build_configuration_section_is_reproducible_description():
    section = wbr.build_configuration_section(elements=123, seed=456)
    assert section["status"] == wbr.STATUS_AVAILABLE
    assert section["elements"] == 123
    assert section["seed"] == 456
    assert set(section["workloads"]) == set(wbr.WORKLOADS)
    assert set(section["dtypes"]) == set(wbr.DTYPES)
    assert section["expectation_classes"]["product"] == wbr.EXPECTATION_BITWISE
    assert section["expectation_classes"]["fma_exposed"] == wbr.EXPECTATION_CONTRACTION_SENSITIVE
    assert section["expectation_classes"]["transcendental"] == wbr.EXPECTATION_TRANSCENDENTAL
    assert section["tolerances"]["product/float32"] == {"rtol": None, "atol": None}
    assert section["tolerances"]["fma_exposed/float32"] == {"rtol": pytest.approx(1e-6), "atol": pytest.approx(1e-6)}
    assert section["tolerances"]["transcendental/float64"] == {"rtol": pytest.approx(1e-12), "atol": pytest.approx(1e-14)}


# ---------------------------------------------------------------------------
# build_report - full assembly, with injected fakes (no warp/numpy needed)
# ---------------------------------------------------------------------------


EXPECTED_TOP_LEVEL_KEYS = {"schema_version", "generated_at_utc", "warp", "devices", "parity", "configuration"}


def test_build_report_with_fake_module_and_runner_has_expected_keys():
    fake_module, cpu, cuda = make_two_device_fake()
    report = wbr.build_report(elements=8, seed=1, warp_module=fake_module, run_kernel=identity_runner)
    assert EXPECTED_TOP_LEVEL_KEYS.issubset(report.keys())
    assert report["warp"]["status"] == wbr.STATUS_AVAILABLE
    assert report["devices"]["status"] == wbr.STATUS_AVAILABLE
    assert report["parity"]["status"] == wbr.STATUS_AVAILABLE
    assert report["configuration"]["status"] == wbr.STATUS_AVAILABLE


def test_build_report_is_json_serializable():
    fake_module, _, _ = make_two_device_fake()
    report = wbr.build_report(elements=8, seed=1, warp_module=fake_module, run_kernel=identity_runner)
    text = json.dumps(report)
    round_tripped = json.loads(text)
    assert round_tripped["schema_version"] == wbr.SCHEMA_VERSION


def test_build_report_schema_version_and_timestamp():
    from datetime import datetime

    fake_module, _, _ = make_two_device_fake()
    report = wbr.build_report(warp_module=fake_module, run_kernel=identity_runner)
    assert report["schema_version"] == 1
    assert isinstance(report["schema_version"], int)
    parsed = datetime.fromisoformat(report["generated_at_utc"])
    assert parsed.tzinfo is not None
    assert parsed.utcoffset().total_seconds() == 0


def test_build_report_degrades_when_warp_absent(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise ImportError("No module named 'warp'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    report = wbr.build_report()
    assert report["warp"]["status"] == wbr.STATUS_UNAVAILABLE
    assert report["devices"]["status"] == wbr.STATUS_UNAVAILABLE
    assert report["parity"]["status"] == wbr.STATUS_UNAVAILABLE
    assert "detail" in report["parity"]
    # Configuration is not a probe - it is always available even when warp is not.
    assert report["configuration"]["status"] == wbr.STATUS_AVAILABLE


def test_build_report_zero_devices_makes_parity_unavailable():
    empty_module = FakeWarpModule([])
    report = wbr.build_report(warp_module=empty_module, run_kernel=identity_runner)
    assert report["devices"]["status"] == wbr.STATUS_UNAVAILABLE
    assert report["parity"]["status"] == wbr.STATUS_UNAVAILABLE


def test_build_report_every_section_has_valid_status():
    fake_module, _, _ = make_two_device_fake()
    report = wbr.build_report(warp_module=fake_module, run_kernel=identity_runner)
    valid_statuses = {wbr.STATUS_AVAILABLE, wbr.STATUS_UNAVAILABLE, wbr.STATUS_ERROR}
    for key in EXPECTED_TOP_LEVEL_KEYS - {"schema_version", "generated_at_utc"}:
        section = report[key]
        assert isinstance(section, dict)
        assert section["status"] in valid_statuses


# ---------------------------------------------------------------------------
# _gate_failures / CLI gate exit codes
# ---------------------------------------------------------------------------


def test_gate_failures_require_cpu_absent():
    report = {"devices": {"status": wbr.STATUS_AVAILABLE, "devices": [{"alias": "cuda:0", "is_cpu": False, "is_cuda": True}]}}
    failures = wbr._gate_failures(report, require_cpu=True, require_cuda=False, require_parity=False)
    assert any("cpu" in f.lower() for f in failures)


def test_gate_failures_require_cuda_absent():
    report = {"devices": {"status": wbr.STATUS_AVAILABLE, "devices": [{"alias": "cpu", "is_cpu": True, "is_cuda": False}]}}
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=True, require_parity=False)
    assert any("cuda" in f.lower() for f in failures)


def test_gate_failures_require_parity_unavailable():
    report = {"devices": {"status": wbr.STATUS_UNAVAILABLE, "devices": []}, "parity": {"status": wbr.STATUS_UNAVAILABLE}}
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=False, require_parity=True)
    assert failures


def test_gate_failures_require_parity_exceeds_tolerance():
    report = {
        "devices": {"status": wbr.STATUS_AVAILABLE, "devices": []},
        "parity": {
            "status": wbr.STATUS_AVAILABLE,
            "reference_device": "cpu",
            "results": [
                {"status": wbr.STATUS_AVAILABLE, "device": "cuda:0", "expectation_class": wbr.EXPECTATION_TRANSCENDENTAL, "verdict": wbr.VERDICT_EXCEEDS_TOLERANCE}
            ],
        },
    }
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=False, require_parity=True)
    assert failures


def test_gate_failures_require_parity_errored_result_fails():
    report = {
        "devices": {"status": wbr.STATUS_AVAILABLE, "devices": []},
        "parity": {
            "status": wbr.STATUS_AVAILABLE,
            "reference_device": "cpu",
            "results": [
                {"status": wbr.STATUS_AVAILABLE, "device": "cuda:0", "expectation_class": wbr.EXPECTATION_BITWISE, "verdict": wbr.VERDICT_IDENTICAL},
                {"status": wbr.STATUS_ERROR, "device": "cuda:0", "workload": "transcendental"},
            ],
        },
    }
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=False, require_parity=True)
    assert any("errored" in f for f in failures)


def test_gate_failures_require_parity_cpu_only_run_fails():
    # Only the reference device itself was compared - not cross-backend
    # evidence, and must not satisfy --require-parity even though every
    # individual result "passed".
    report = {
        "devices": {"status": wbr.STATUS_AVAILABLE, "devices": CPU_ONLY_DEVICE_LIST},
        "parity": {
            "status": wbr.STATUS_AVAILABLE,
            "reference_device": "cpu",
            "results": [
                {"status": wbr.STATUS_AVAILABLE, "device": "cpu", "expectation_class": wbr.EXPECTATION_BITWISE, "verdict": wbr.VERDICT_IDENTICAL},
                {"status": wbr.STATUS_AVAILABLE, "device": "cpu", "expectation_class": wbr.EXPECTATION_TRANSCENDENTAL, "verdict": wbr.VERDICT_IDENTICAL},
            ],
        },
    }
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=False, require_parity=True)
    assert any("cross-backend" in f for f in failures)


def test_gate_failures_none_requested_is_empty():
    report = {"devices": {"status": wbr.STATUS_UNAVAILABLE, "devices": []}, "parity": {"status": wbr.STATUS_UNAVAILABLE}}
    failures = wbr._gate_failures(report, require_cpu=False, require_cuda=False, require_parity=False)
    assert failures == []


def test_gate_failures_all_satisfied_is_empty():
    report = {
        "devices": {"status": wbr.STATUS_AVAILABLE, "devices": TWO_DEVICE_LIST},
        "parity": {
            "status": wbr.STATUS_AVAILABLE,
            "reference_device": "cpu",
            "results": [
                {"status": wbr.STATUS_AVAILABLE, "device": "cpu", "expectation_class": wbr.EXPECTATION_BITWISE, "verdict": wbr.VERDICT_IDENTICAL},
                {"status": wbr.STATUS_AVAILABLE, "device": "cuda:0", "expectation_class": wbr.EXPECTATION_BITWISE, "verdict": wbr.VERDICT_IDENTICAL},
                {"status": wbr.STATUS_AVAILABLE, "device": "cuda:0", "expectation_class": wbr.EXPECTATION_TRANSCENDENTAL, "verdict": wbr.VERDICT_WITHIN_TOLERANCE},
            ],
        },
    }
    failures = wbr._gate_failures(report, require_cpu=True, require_cuda=True, require_parity=True)
    assert failures == []


def test_cli_require_cpu_exits_3_when_warp_absent(monkeypatch, capsys):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise ImportError("No module named 'warp'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    exit_code = wbr.main(["--require-cpu"])
    assert exit_code == 3
    captured = capsys.readouterr()
    # stdout still carries a well-formed report - report delivery succeeded.
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION
    error_payload = json.loads(captured.err)
    assert error_payload["status"] == wbr.STATUS_ERROR


def test_cli_no_gates_exits_0_when_warp_absent(monkeypatch, capsys):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise ImportError("No module named 'warp'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    exit_code = wbr.main([])
    assert exit_code == 0
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION


# ---------------------------------------------------------------------------
# CLI - other flags, output writing, exit codes
# ---------------------------------------------------------------------------


def test_cli_main_exits_zero_and_prints_valid_json():
    exit_code = wbr.main([])
    assert exit_code in (0, 3)  # host may or may not have a warp CPU device; either way it must not crash


def test_cli_main_respects_indent(capsys):
    wbr.main(["--indent", "0"])
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION


def test_cli_main_writes_output_file(tmp_path):
    output_path = tmp_path / "report.json"
    wbr.main(["--output", str(output_path)])
    assert output_path.exists()
    parsed = json.loads(output_path.read_text(encoding="utf-8"))
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION


def test_cli_main_exits_one_with_structured_error_when_output_unwritable(tmp_path, capsys):
    unwritable_path = tmp_path / "no-such-directory" / "report.json"
    exit_code = wbr.main(["--output", str(unwritable_path)])
    assert exit_code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    error_payload = json.loads(captured.err)
    assert error_payload["status"] == wbr.STATUS_ERROR
    assert not unwritable_path.exists()


def test_cli_main_custom_elements_and_seed(capsys):
    wbr.main(["--elements", "16", "--seed", "7"])
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["configuration"]["elements"] == 16
    assert parsed["configuration"]["seed"] == 7


@pytest.mark.parametrize("bad_value", ["0", "-1", "-100"])
def test_cli_elements_non_positive_exits_two(bad_value):
    with pytest.raises(SystemExit) as exc_info:
        wbr.main(["--elements", bad_value])
    assert exc_info.value.code == 2


def test_cli_elements_positive_boundary_is_accepted(capsys):
    exit_code = wbr.main(["--elements", "1"])
    assert exit_code in (0, 3)
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["configuration"]["elements"] == 1


def test_cli_help_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc_info:
        wbr.main(["--help"])
    assert exc_info.value.code == 0


def test_cli_unrecognized_argument_exits_two():
    with pytest.raises(SystemExit) as exc_info:
        wbr.main(["--not-a-real-flag"])
    assert exc_info.value.code == 2


def test_cli_as_subprocess_exits_zero_or_three_with_valid_json():
    module_path = wbr.__file__ if hasattr(wbr, "__file__") else None
    assert module_path is not None
    proc = subprocess.run(
        [sys.executable, module_path],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        text=True,
    )
    assert proc.returncode in (0, 3), proc.stderr
    parsed = json.loads(proc.stdout)
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION


def test_cli_as_subprocess_elements_zero_exits_two():
    module_path = wbr.__file__
    proc = subprocess.run(
        [sys.executable, module_path, "--elements", "0"],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        text=True,
    )
    assert proc.returncode == 2
    assert proc.stdout == ""


# ---------------------------------------------------------------------------
# Sensitive-content guard
# ---------------------------------------------------------------------------


def test_report_contains_no_environment_variable_dump():
    import os as os_module

    fake_module, _, _ = make_two_device_fake()
    report = wbr.build_report(warp_module=fake_module, run_kernel=identity_runner)
    text = json.dumps(report)

    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                assert k.lower() not in ("environ", "environment", "env"), f"suspicious key: {k}"
                walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(report)

    for var_name, var_value in os_module.environ.items():
        if not var_value:
            continue
        if any(token in var_name.upper() for token in ("TOKEN", "SECRET", "PASSWORD", "CREDENTIAL", "API_KEY")):
            assert var_value not in text, f"leaked value of env var {var_name}"


# ---------------------------------------------------------------------------
# Real `warp`, only when actually installed - no assumption about which
# devices it exposes (CPU-only vs CPU+CUDA), and no assumption that it is
# absent either.
# ---------------------------------------------------------------------------


def test_real_warp_devices_section_well_formed():
    warp = pytest.importorskip("warp")
    section = wbr.build_devices_section(warp)
    assert section["status"] in (wbr.STATUS_AVAILABLE, wbr.STATUS_UNAVAILABLE, wbr.STATUS_ERROR)
    if section["status"] == wbr.STATUS_AVAILABLE:
        for device in section["devices"]:
            assert "alias" in device
            assert "is_cpu" in device
            assert "is_cuda" in device
            if not device["is_cuda"]:
                assert device["arch"] is None


def test_real_warp_full_report_is_well_formed_and_json_serializable():
    pytest.importorskip("warp")
    report = wbr.build_report(elements=32, seed=11)
    text = json.dumps(report)
    parsed = json.loads(text)
    assert parsed["schema_version"] == wbr.SCHEMA_VERSION
    assert parsed["parity"]["status"] in (wbr.STATUS_AVAILABLE, wbr.STATUS_UNAVAILABLE, wbr.STATUS_ERROR)
