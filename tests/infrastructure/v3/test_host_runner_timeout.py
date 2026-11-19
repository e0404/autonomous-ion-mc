"""V3-INFRA-008: the controlled host runner's timeout ceiling.

Raised from 7200 s to 14400 s under the operator's authorization
IR-20261118-201459-7CFC95 (option A) so that one frozen V3-013 ionmc sample
(1e6 He histories, warp-cpu float64, measured at up to ~8600 s) completes in a
single controlled run. The bound is still enforced fail-closed before any
run directory is created or any process is launched.
"""

from __future__ import annotations

import pytest

from infrastructure.host_runner import host_runner


def test_timeout_ceiling_is_four_hours() -> None:
    assert host_runner.MAX_TIMEOUT_SECONDS == 14400


@pytest.mark.parametrize("timeout", [0, -1, 14401, 2 * 14400])
def test_timeouts_outside_the_bound_are_rejected_before_any_side_effect(
    monkeypatch: pytest.MonkeyPatch, timeout: int
) -> None:
    monkeypatch.setattr(
        host_runner, "create_run_dir", lambda *a, **k: pytest.fail("created run")
    )
    monkeypatch.setattr(
        host_runner.subprocess, "run", lambda *a, **k: pytest.fail("launched process")
    )
    with pytest.raises(ValueError, match="between 1 and 14400"):
        host_runner.run_validation("V3-TEST", ["true"], timeout)
