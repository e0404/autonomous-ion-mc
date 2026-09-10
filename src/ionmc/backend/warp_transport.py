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
    rng_state0: wp.array(dtype=wp.uint32),
    density: float,
    max_fraction: float,
    max_step_mm: float,
    bin_width_mm: float,
    geom_depth_mm: float,
    energy_cut_mev: float,
    max_steps: int,
    rest_energy_mev: float,
    charge: float,
    za_ratio: float,
    straggling: int,
    straggling_floor_mev: float,
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    n: int,
    n_steps: int,
    n_bins: int,
    edep: wp.array(dtype=wp.float64),
    truncated: wp.array(dtype=int),
    final_z: wp.array(dtype=float),
    final_status: wp.array(dtype=int),
):
    i = wp.tid()
    e = energy0[i]
    z = z0[i]
    w = weight[i]
    rng = rng_state0[i]
    # dynamic (mutable) loop variables: int(...) keeps Warp from treating them
    # as constants inside the while loop (noqa keeps ruff from stripping int()).
    step = int(0)  # noqa: UP018, RUF046
    alive = int(1)  # noqa: UP018, RUF046
    while alive == 1 and step < max_steps:
        # physics step: fractional energy loss, not bin-limited (decision 0010),
        # so the straggling and clamp are unbiased; deposition is split across
        # the bins the step spans (mirrors the reference driver).
        dl_e = transport.energy_loss_step_length(
            e, max_fraction, max_step_mm, density, table_e, table_s, table_d, n, n_steps
        )
        dl = wp.min(dl_e, geom_depth_mm - z)
        de = transport.midpoint_energy_loss(
            e, dl, density, table_e, table_s, table_d, n, n_steps
        )
        if straggling == 1 and e > straggling_floor_mev:
            sigma = transport.bohr_straggling_sigma(
                e, rest_energy_mev, dl, density, za_ratio, charge
            )
            variate = wp.randn(rng)
            de = transport.straggled_energy_loss(de, sigma, variate, e)
        deposit = w * de
        z1 = z + dl
        pos = z
        b = int(wp.floor(z / bin_width_mm))
        inv_dl = 1.0 / wp.max(dl, 1.0e-12)
        while pos < z1 - 1.0e-12:
            bin_end = float(b + 1) * bin_width_mm
            seg_end = wp.min(bin_end, z1)
            if b >= 0 and b < n_bins:
                wp.atomic_add(edep, b, wp.float64(deposit * (seg_end - pos) * inv_dl))
            pos = seg_end
            b = b + 1
        z = z + dl
        e = e - de
        step = step + 1
        if e <= energy_cut_mev:
            dep_bin = int(wp.floor(z / bin_width_mm))
            if dep_bin >= 0 and dep_bin < n_bins:
                wp.atomic_add(edep, dep_bin, wp.float64(w * e))
            e = 0.0
            alive = 0
        if z >= geom_depth_mm:
            alive = 0
    # per-history outcome: status 1 stopped, 2 escaped, 3 truncated
    final_z[i] = z
    if step >= max_steps and alive == 1:
        wp.atomic_add(truncated, 0, 1)
        final_status[i] = 3
    elif z >= geom_depth_mm:
        final_status[i] = 2
    else:
        final_status[i] = 1


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
        rng_state: np.ndarray,
        density: float,
        max_fraction: float,
        max_step_mm: float,
        bin_width_mm: float,
        geom_depth_mm: float,
        energy_cut_mev: float,
        max_steps: int,
        n_bins: int,
        rest_energy_mev: float,
        charge: float,
        za_ratio: float,
        straggling: bool,
        straggling_floor_mev: float,
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray]:
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
        rng: Any = wp.array(
            np.ascontiguousarray(rng_state, dtype=np.uint32),
            dtype=wp.uint32,
            device=d,
        )
        edep = wp.zeros(n_bins, dtype=wp.float64, device=d)
        truncated = wp.zeros(1, dtype=int, device=d)
        final_z = wp.zeros(n_hist, dtype=float, device=d)
        final_status = wp.zeros(n_hist, dtype=int, device=d)
        t = self.tables
        wp.launch(
            csda_depth_dose_kernel,
            dim=n_hist,
            inputs=[
                e0,
                zz,
                ww,
                rng,
                float(density),
                float(max_fraction),
                float(max_step_mm),
                float(bin_width_mm),
                float(geom_depth_mm),
                float(energy_cut_mev),
                int(max_steps),
                float(rest_energy_mev),
                float(charge),
                float(za_ratio),
                (1 if straggling else 0),
                float(straggling_floor_mev),
                t["e"],
                t["s"],
                t["d"],
                self.n,
                self.n_steps,
                int(n_bins),
                edep,
                truncated,
                final_z,
                final_status,
            ],
            device=d,
        )
        wp.synchronize_device(d)
        return (
            edep.numpy().astype(np.float64),  # already float64 on device
            int(truncated.numpy()[0]),
            final_z.numpy().astype(np.float64),
            final_status.numpy().astype(np.int32),
        )
