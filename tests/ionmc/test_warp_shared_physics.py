"""The shared step physics compiled by Warp must match its Python execution and
advance random state across calls (regression for the by-value state defect)."""

import numpy as np
import pytest
import warp as wp
from warp_test_kernels import f_dm_kernel, sequential_losses, two_losses

from ionmc.transport.shared import physics


@pytest.mark.parametrize("dtype", [wp.float32, wp.float64])
def test_consecutive_loss_samples_differ_and_have_correct_moments(dtype):
    out = wp.zeros((4000, 2), dtype=dtype, device="cpu")
    wp.launch(two_losses, dim=4000, inputs=[out, dtype(0.7), dtype(0.01)], device="cpu")
    x = out.numpy().astype(np.float64)
    assert np.all(x[:, 0] != x[:, 1])
    seq = wp.zeros(50000, dtype=dtype, device="cpu")
    wp.launch(
        sequential_losses,
        dim=1,
        inputs=[seq, dtype(0.7), dtype(0.01), 50000],
        device="cpu",
    )
    s = seq.numpy().astype(np.float64)
    assert s.mean() == pytest.approx(0.7, rel=0.005)
    assert s.var() == pytest.approx(0.01, rel=0.05)
    assert len(np.unique(s[:1000])) > 990


@pytest.mark.parametrize("dtype", [wp.float32, wp.float64])
def test_shared_scattering_factor_matches_python(dtype):
    py = physics("python")
    out = wp.zeros(1, dtype=dtype, device="cpu")
    wp.launch(
        f_dm_kernel, dim=1, inputs=[out, dtype(200.0), dtype(280.0)], device="cpu"
    )
    assert float(out.numpy()[0]) == pytest.approx(py.f_dm(200.0, 280.0), rel=1e-6)
