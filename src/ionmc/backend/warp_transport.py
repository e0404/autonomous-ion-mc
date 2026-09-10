"""Warp kernel for continuous-slowing-down proton depth-dose transport.

One history per thread with an in-kernel step loop and atomic scoring into a
1-D depth grid (decision ``0009``). The kernel contains no physics: the step
length and the midpoint energy loss come from the shared-source functions in
:mod:`ionmc.physics.transport`, which also run as the float64 reference driver.
Runs in float32 on the CPU and CUDA devices.
"""

# mypy: disable-error-code="valid-type"

from __future__ import annotations

from typing import Any

import numpy as np
import warp as wp

from ionmc.data.stopping_tables import StoppingTable
from ionmc.physics import transport


@wp.kernel
def csda_depth_dose_kernel(
    energy0: wp.array(dtype=float),
    z0: wp.array(dtype=float),
    weight: wp.array(dtype=float),
    density: float,
    max_fraction: float,
    max_step_mm: float,
    bin_width_mm: float,
    geom_depth_mm: float,
    energy_cut_mev: float,
    max_steps: int,
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    n: int,
    n_steps: int,
    n_bins: int,
    edep: wp.array(dtype=float),
    truncated: wp.array(dtype=int),
):
    i = wp.tid()
    e = energy0[i]
    z = z0[i]
    w = weight[i]
    # dynamic (mutable) loop variables: int(...) keeps Warp from treating them
    # as constants inside the while loop (noqa keeps ruff from stripping int()).
    step = int(0)  # noqa: UP018, RUF046
    alive = int(1)  # noqa: UP018, RUF046
    while alive == 1 and step < max_steps:
        cur_bin = int(wp.floor(z / bin_width_mm))
        boundary = float(cur_bin + 1) * bin_width_mm
        dl_e = transport.energy_loss_step_length(
            e, max_fraction, max_step_mm, density, table_e, table_s, table_d, n, n_steps
        )
        dl = wp.min(dl_e, boundary - z + 1.0e-6)
        dl = wp.min(dl, geom_depth_mm - z)
        de = transport.midpoint_energy_loss(
            e, dl, density, table_e, table_s, table_d, n, n_steps
        )
        if cur_bin >= 0 and cur_bin < n_bins:
            wp.atomic_add(edep, cur_bin, w * de)
        z = z + dl
        e = e - de
        step = step + 1
        if e <= energy_cut_mev:
            dep_bin = int(wp.floor(z / bin_width_mm))
            if dep_bin >= 0 and dep_bin < n_bins:
                wp.atomic_add(edep, dep_bin, w * e)
            e = 0.0
            alive = 0
        if z >= geom_depth_mm:
            alive = 0
    if step >= max_steps and alive == 1:
        wp.atomic_add(truncated, 0, 1)


class DepthDoseKernel:
    """Launches the CSDA depth-dose kernel for one table on one device."""

    def __init__(self, table: StoppingTable, device: str = "cpu") -> None:
        wp.init()
        self.device = device
        self.table = table
        self.n = int(table.size)
        self.n_steps = int(table.bisection_steps)
        self.tables: dict[str, Any] = {
            name: wp.array(
                np.asarray(getattr(table, attr), dtype=np.float32),
                dtype=float,
                device=device,
            )
            for name, attr in (
                ("e", "energy_mev"),
                ("s", "stopping_mev_cm2_per_g"),
                ("d", "slope"),
            )
        }

    def run(
        self,
        energy0: np.ndarray,
        z0: np.ndarray,
        weight: np.ndarray,
        density: float,
        max_fraction: float,
        max_step_mm: float,
        bin_width_mm: float,
        geom_depth_mm: float,
        energy_cut_mev: float,
        max_steps: int,
        n_bins: int,
    ) -> tuple[np.ndarray, int]:
        """Return (energy deposited per bin [MeV], number of truncated histories)."""
        d = self.device
        n_hist = int(energy0.shape[0])
        e0: Any = wp.array(
            np.ascontiguousarray(energy0, dtype=np.float32), dtype=float, device=d
        )
        zz: Any = wp.array(
            np.ascontiguousarray(z0, dtype=np.float32), dtype=float, device=d
        )
        ww: Any = wp.array(
            np.ascontiguousarray(weight, dtype=np.float32), dtype=float, device=d
        )
        edep = wp.zeros(n_bins, dtype=float, device=d)
        truncated = wp.zeros(1, dtype=int, device=d)
        t = self.tables
        wp.launch(
            csda_depth_dose_kernel,
            dim=n_hist,
            inputs=[
                e0,
                zz,
                ww,
                float(density),
                float(max_fraction),
                float(max_step_mm),
                float(bin_width_mm),
                float(geom_depth_mm),
                float(energy_cut_mev),
                int(max_steps),
                t["e"],
                t["s"],
                t["d"],
                self.n,
                self.n_steps,
                int(n_bins),
                edep,
                truncated,
            ],
            device=d,
        )
        wp.synchronize_device(d)
        return edep.numpy().astype(np.float64), int(truncated.numpy()[0])
