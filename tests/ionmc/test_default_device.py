"""The test session pins Warp's default device to the CPU (host independence)."""

import warp as wp


def test_default_device_is_cpu_in_the_suite() -> None:
    assert wp.get_device().is_cpu
    arr = wp.array([1.0, 2.0], dtype=wp.float64)  # no device given: lands on the default
    assert arr.device.is_cpu
