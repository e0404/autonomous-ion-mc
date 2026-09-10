"""Warp kernels evaluating the shared-source stopping-power physics.

This module is the only place in the stopping-power stack that imports Warp
unconditionally. The kernels contain no physics of their own: they unpack a
parameter struct, read the per-thread energy and call the *same*
``ionmc.physics.stopping`` functions that the Python and numpy paths execute.
Kernels run in float32 (decision ``0005``).
"""

# mypy: disable-error-code="valid-type"
# (Warp kernel signatures use ``wp.array(dtype=...)`` call expressions as
# annotations, which mypy cannot interpret; Warp requires this form.)

from __future__ import annotations

from typing import Any

import numpy as np
import warp as wp

from ionmc.data.stopping_tables import StoppingTable
from ionmc.physics import stopping, tabulated
from ionmc.stopping_power import StoppingPowerParameters


@wp.struct
class StoppingParams:
    """Scalar parameters of :func:`ionmc.physics.stopping.mass_stopping_power`."""

    rest_energy: float
    charge: float
    za_ratio: float
    i_mev: float
    cbar: float
    x0: float
    x1: float
    a: float
    mexp: float
    delta0: float
    shell_m2: float
    shell_m3: float
    n_elem: int
    n_table: int
    use_density: float
    use_shell: float
    use_barkas: float
    use_bloch: float


@wp.kernel
def mass_stopping_power_kernel(
    energies: wp.array(dtype=float),
    p: StoppingParams,
    elem_z: wp.array(dtype=float),
    elem_b: wp.array(dtype=float),
    elem_f: wp.array(dtype=float),
    table_w: wp.array(dtype=float),
    table_f: wp.array(dtype=float),
    out: wp.array(dtype=float),
):
    i = wp.tid()
    out[i] = stopping.mass_stopping_power(
        energies[i],
        p.rest_energy,
        p.charge,
        p.za_ratio,
        p.i_mev,
        p.cbar,
        p.x0,
        p.x1,
        p.a,
        p.mexp,
        p.delta0,
        p.shell_m2,
        p.shell_m3,
        elem_z,
        elem_b,
        elem_f,
        p.n_elem,
        table_w,
        table_f,
        p.n_table,
        p.use_density,
        p.use_shell,
        p.use_barkas,
        p.use_bloch,
    )


@wp.kernel
def csda_range_kernel(
    energies: wp.array(dtype=float),
    energy_floor: float,
    n_steps: int,
    p: StoppingParams,
    elem_z: wp.array(dtype=float),
    elem_b: wp.array(dtype=float),
    elem_f: wp.array(dtype=float),
    table_w: wp.array(dtype=float),
    table_f: wp.array(dtype=float),
    out: wp.array(dtype=float),
):
    i = wp.tid()
    out[i] = stopping.csda_range_increment(
        energy_floor,
        energies[i],
        n_steps,
        p.rest_energy,
        p.charge,
        p.za_ratio,
        p.i_mev,
        p.cbar,
        p.x0,
        p.x1,
        p.a,
        p.mexp,
        p.delta0,
        p.shell_m2,
        p.shell_m3,
        elem_z,
        elem_b,
        elem_f,
        p.n_elem,
        table_w,
        table_f,
        p.n_table,
        p.use_density,
        p.use_shell,
        p.use_barkas,
        p.use_bloch,
    )


class StoppingPowerKernels:
    """Launches the stopping-power kernels for one parameter set on one device."""

    def __init__(
        self, parameters: StoppingPowerParameters, device: str = "cpu"
    ) -> None:
        wp.init()
        self.device = device
        self.parameters = parameters
        p = StoppingParams()
        p.rest_energy = parameters.rest_energy_mev
        p.charge = parameters.charge
        p.za_ratio = parameters.za_ratio
        p.i_mev = parameters.i_mev
        p.cbar = parameters.cbar
        p.x0 = parameters.x0
        p.x1 = parameters.x1
        p.a = parameters.a
        p.mexp = parameters.mexp
        p.delta0 = parameters.delta0
        p.shell_m2 = parameters.shell_m2
        p.shell_m3 = parameters.shell_m3
        p.n_elem = parameters.n_elem
        p.n_table = parameters.n_table
        p.use_density = parameters.use_density
        p.use_shell = parameters.use_shell
        p.use_barkas = parameters.use_barkas
        p.use_bloch = parameters.use_bloch
        self.struct = p
        self.tables: dict[str, Any] = {
            name: wp.array(
                np.asarray(getattr(parameters, name), dtype=np.float32),
                dtype=float,
                device=device,
            )
            for name in ("elem_z", "elem_b", "elem_f", "table_w", "table_f")
        }

    def _launch(
        self, kernel: wp.Kernel, energies: np.ndarray, extra: list
    ) -> np.ndarray:
        e32 = np.ascontiguousarray(np.asarray(energies, dtype=np.float32))
        e_arr: Any = wp.array(e32, dtype=float, device=self.device)
        out = wp.zeros(e32.shape[0], dtype=float, device=self.device)
        t = self.tables
        inputs = [
            e_arr,
            *extra,
            self.struct,
            t["elem_z"],
            t["elem_b"],
            t["elem_f"],
            t["table_w"],
            t["table_f"],
            out,
        ]
        wp.launch(kernel, dim=e32.shape[0], inputs=inputs, device=self.device)
        wp.synchronize_device(self.device)
        return out.numpy().astype(np.float64)

    def mass_stopping_power(self, energies: np.ndarray) -> np.ndarray:
        """Mass stopping power [MeV cm^2/g] evaluated in float32 on the device."""
        return self._launch(mass_stopping_power_kernel, energies, [])

    def csda_range(
        self, energies: np.ndarray, energy_floor: float, n_steps: int
    ) -> np.ndarray:
        """CSDA range increment [g/cm^2] evaluated in float32 on the device."""
        return self._launch(
            csda_range_kernel, energies, [float(energy_floor), int(n_steps)]
        )


@wp.kernel
def tabulated_stopping_power_kernel(
    energies: wp.array(dtype=float),
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    n: int,
    n_steps: int,
    out: wp.array(dtype=float),
):
    i = wp.tid()
    out[i] = tabulated.tabulated_mass_stopping_power(
        energies[i], table_e, table_s, table_d, n, n_steps
    )


@wp.kernel
def tabulated_csda_range_kernel(
    energies: wp.array(dtype=float),
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    table_range: wp.array(dtype=float),
    n: int,
    n_steps: int,
    out: wp.array(dtype=float),
):
    i = wp.tid()
    out[i] = tabulated.tabulated_csda_range(
        energies[i], table_e, table_s, table_d, table_range, n, n_steps
    )


class TabulatedKernels:
    """Launches the table-lookup kernels for one prepared table on one device."""

    def __init__(self, table: StoppingTable, device: str = "cpu") -> None:
        wp.init()
        self.device = device
        self.table = table
        self.n = int(table.size)
        self.n_steps = int(table.bisection_steps)
        self.arrays: dict[str, Any] = {
            name: wp.array(
                np.asarray(getattr(table, attr), dtype=np.float32),
                dtype=float,
                device=device,
            )
            for name, attr in (
                ("e", "energy_mev"),
                ("s", "stopping_mev_cm2_per_g"),
                ("d", "slope"),
                ("range", "csda_range_g_per_cm2"),
            )
        }

    def _launch(
        self, kernel: wp.Kernel, energies: np.ndarray, tables: list
    ) -> np.ndarray:
        e32 = np.ascontiguousarray(np.asarray(energies, dtype=np.float32))
        e_arr: Any = wp.array(e32, dtype=float, device=self.device)
        out = wp.zeros(e32.shape[0], dtype=float, device=self.device)
        inputs = [e_arr, *tables, self.n, self.n_steps, out]
        wp.launch(kernel, dim=e32.shape[0], inputs=inputs, device=self.device)
        wp.synchronize_device(self.device)
        return out.numpy().astype(np.float64)

    def mass_stopping_power(self, energies: np.ndarray) -> np.ndarray:
        """Mass stopping power [MeV cm^2/g] from the table, float32 on the device."""
        a = self.arrays
        return self._launch(
            tabulated_stopping_power_kernel, energies, [a["e"], a["s"], a["d"]]
        )

    def csda_range(self, energies: np.ndarray) -> np.ndarray:
        """CSDA range [g/cm^2] from the table's floor energy, float32 on the device."""
        a = self.arrays
        return self._launch(
            tabulated_csda_range_kernel, energies, [a["e"], a["s"], a["d"], a["range"]]
        )
