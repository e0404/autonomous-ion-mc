"""Device-side packing of the nuclear table for the Warp kernels (decision 0041, V3-005B).

:func:`pack_nuclear` flattens the runtime rows of a loaded ``NuclearTable`` and of the per-material
rows (``NuclearTable.material_rows``) into the flat arrays that ``make_nuclear(real).sample_event``
and ``choose_target`` index (layouts in the docstring of ``sample_event``): table grid, per-target
multiplicity means, bin edges, pre-compound fractions and recoil energies, the per-target constant
rows (masses, separation energies, ``z_t``, ``a_t``), the residual-mass table, and per material the
total ``Sigma_mass``, the window and end-of-range majorants and the cumulative target ``Sigma``.
:class:`NuclearDevice` uploads them as ``wp.array`` of the kernel precision (float32 or float64;
integer index rows as int32) on an explicit device, as ``TransportTables.to_warp`` does for the
transport tables, and carries a ``sha256`` identity of the uploaded bytes
(:meth:`NuclearDevice.device_sha256` re-hashes the bytes read back from the device). The same
packed host arrays (float64) drive the pure-Python twin of the sampler, so a parity test compares
one data set through two execution paths.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import numpy as np
import warp as wp
from numpy.typing import NDArray

from ionmc.data.ame import AmeEntry
from ionmc.materials import Material
from ionmc.nuclear.build import bounds_from_array
from ionmc.nuclear.events import build_event_model, event_constants
from ionmc.nuclear.tables import DEFAULT_F_E, NuclearTable
from ionmc.physics.nuclear import DA_MAX, DZ_MAX, N_BINS, N_SPECIES, TCONST_STRIDE

REAL_FIELDS = (
    "grid", "lam", "edges", "rpre", "recoil", "tconst", "m_res",
    "sigma", "sigma_win", "sigma_end", "cum_sigma",
)  # fmt: skip
"""Arrays stored in the kernel precision."""
INT_FIELDS = ("mat_ntargets", "mat_target")
"""Integer index arrays (int32): targets per material and table-target index per material row."""


@dataclass(frozen=True)
class NuclearHost:
    """Packed float64/int64 arrays of a table and its materials (flat layouts of ``sample_event``).

    ``n_grid`` nodes, ``n_targets`` table targets, ``kmax`` the largest number of targets of a
    material, ``n_materials``; ``bounds`` the table's energy and path bounds (floats)."""

    table_id: str
    n_grid: int
    n_targets: int
    kmax: int
    n_materials: int
    arrays: dict[str, NDArray[Any]]
    bounds: dict[str, Any]


def pack_nuclear(
    table: NuclearTable,
    materials: tuple[Material, ...],
    ame: dict[tuple[int, int], AmeEntry] | None = None,
    f_e: float = DEFAULT_F_E,
    rows: Any = None,
) -> NuclearHost:
    """Pack ``table`` and the per-material rows of ``materials`` (``ame``: AME2020 entries for
    the event constants; default: the process-cached verified file of the data cache; ``rows``: the
    ``MaterialNuclear`` rows of ``materials``, default ``table.material_rows(m, f_e)``)."""
    if ame is None:
        ame = cached_ame2020()
    a = table.arrays
    grid = np.asarray(a["grid_e_mev"], dtype=np.float64)
    n = grid.size
    targets = table.info["targets"]
    n_t = len(targets)
    tconst = np.zeros((n_t, TCONST_STRIDE))
    m_res = np.full((n_t, DZ_MAX + 1, DA_MAX + 1), np.inf)
    for t, info in enumerate(targets):
        model = build_event_model(ame, int(info["z"]), int(info["a"]))
        tconst[t] = event_constants(model)
        m_res[t] = model.m_res_mev
    rows = [table.material_rows(m, f_e) for m in materials] if rows is None else list(rows)
    n_m = max(len(rows), 1)
    kmax = max([len(r.target_index) for r in rows] + [1])
    sigma, win, end = (np.zeros((n_m, n)) for _ in range(3))
    cum = np.zeros((n_m, kmax, n))
    n_k = np.zeros(n_m, dtype=np.int64)
    tgt = np.full((n_m, kmax), -1, dtype=np.int64)
    for i, r in enumerate(rows):
        sigma[i], win[i], end[i] = r.sigma_mass_cm2_g, r.sigma_hat_window, r.sigma_hat_end
        k = len(r.target_index)
        n_k[i] = k
        tgt[i, :k] = r.target_index
        cum[i, :k] = r.cum_sigma_mass_cm2_g
    arrays: dict[str, NDArray[Any]] = {
        "grid": grid,
        "lam": np.ascontiguousarray(a["lam"], dtype=np.float64).ravel(),
        "edges": np.ascontiguousarray(a["edges_mev"], dtype=np.float64).ravel(),
        "rpre": np.ascontiguousarray(a["r_pre"], dtype=np.float64).ravel(),
        "recoil": np.ascontiguousarray(a["recoil_t_cm_mev"], dtype=np.float64).ravel(),
        "tconst": tconst.ravel(),
        "m_res": m_res.ravel(),
        "sigma": sigma.ravel(),
        "sigma_win": win.ravel(),
        "sigma_end": end.ravel(),
        "cum_sigma": cum.ravel(),
        "mat_ntargets": n_k,
        "mat_target": tgt.ravel(),
    }
    expected = {
        "lam": n_t * N_SPECIES * n, "edges": n_t * N_SPECIES * n * (N_BINS + 1),
        "rpre": n_t * N_SPECIES * n * N_BINS, "recoil": n_t * n,
    }  # fmt: skip
    for name, size in expected.items():
        if arrays[name].size != size:
            raise ValueError(
                f"nuclear table array {name!r} has size {arrays[name].size}, not {size}"
            )
    return NuclearHost(
        table.table_id,
        n,
        n_t,
        kmax,
        len(rows),
        arrays,
        dict(bounds_from_array(a["bounds"])),
    )


_AME_CACHE: dict[str, dict[tuple[int, int], AmeEntry]] = {}
_DEVICE_CACHE: dict[tuple[Any, ...], NuclearDevice] = {}
_MAX_CACHED_DEVICES = 4


def cached_ame2020() -> dict[tuple[int, int], AmeEntry]:
    """AME2020 entries of the verified cache file, verified and parsed ONCE per process.

    Key: the resolved cache directory (``cache.verify`` is content-addressed). The returned
    mapping is shared: callers must not modify it."""
    from ionmc.data import cache
    from ionmc.data.ame import load_ame2020

    cdir = cache.resolve_cache_dir(None)
    key = str(cdir)
    if key not in _AME_CACHE:
        path = cache.verify("ame2020-mass", cdir)
        _AME_CACHE.clear()  # one directory at a time: the cache never grows
        _AME_CACHE[key] = load_ame2020(path.read_text(encoding="ascii"))
    return _AME_CACHE[key]


def _np_dtype(real: Any) -> Any:
    if real is wp.float32:
        return np.float32
    if real is wp.float64:
        return np.float64
    raise ValueError("real must be wp.float32 or wp.float64")


def _hash_arrays(table_id: str, named: dict[str, NDArray[Any]]) -> str:
    h = hashlib.sha256()
    h.update(table_id.encode())
    for name in (*REAL_FIELDS, *INT_FIELDS):
        arr = np.ascontiguousarray(named[name])
        h.update(f"{name}|{arr.dtype.str}|{arr.shape}".encode())
        h.update(arr.tobytes())
    return h.hexdigest()


def host_sha256(host: NuclearHost, precision: str) -> str:
    """``NuclearDevice.sha256`` of ``host`` in ``precision`` ("float32" or "float64") computed
    without a device (the same bytes: arrays cast to the kernel precision, integer rows int32)."""
    np_real = {"float32": np.float32, "float64": np.float64}[precision]
    cast: dict[str, NDArray[Any]] = {
        k: np.ascontiguousarray(host.arrays[k], dtype=np_real) for k in REAL_FIELDS
    }
    cast.update({k: np.ascontiguousarray(host.arrays[k], dtype=np.int32) for k in INT_FIELDS})
    return _hash_arrays(host.table_id, cast)


class NuclearDevice:
    """The packed nuclear arrays on one device in one precision.

    Attributes are ``wp.array`` objects named as :data:`REAL_FIELDS` / :data:`INT_FIELDS` plus
    ``n_grid``, ``n_targets``, ``kmax``, ``n_materials``, ``bounds``, ``real``, ``device``,
    ``table_id`` and ``sha256`` (identity of the bytes uploaded: the table id and, per array in
    the fixed order of the two field tuples, name, dtype, shape and bytes after the cast to the
    kernel precision)."""

    def __init__(self, host: NuclearHost, real: Any, device: str) -> None:
        np_real = _np_dtype(real)
        cast: dict[str, NDArray[Any]] = {
            k: np.ascontiguousarray(host.arrays[k], dtype=np_real) for k in REAL_FIELDS
        }
        cast.update({k: np.ascontiguousarray(host.arrays[k], dtype=np.int32) for k in INT_FIELDS})
        self.real = real
        self.device = device
        self.table_id = host.table_id
        self.n_grid, self.n_targets = host.n_grid, host.n_targets
        self.kmax, self.n_materials = host.kmax, host.n_materials
        self.bounds = dict(host.bounds)
        self.sha256 = _hash_arrays(host.table_id, cast)
        for k in REAL_FIELDS:
            setattr(self, k, wp.array(cast[k], dtype=real, device=device))
        for k in INT_FIELDS:
            setattr(self, k, wp.array(cast[k], dtype=wp.int32, device=device))

    @classmethod
    def from_table(
        cls,
        table: NuclearTable,
        materials: tuple[Material, ...],
        *,
        real: Any,
        device: str,
        ame: dict[tuple[int, int], AmeEntry] | None = None,
        f_e: float = DEFAULT_F_E,
        rows: Any = None,
    ) -> NuclearDevice:
        """Pack and upload ``table`` with the rows of ``materials`` (see :func:`pack_nuclear`)."""
        return cls(pack_nuclear(table, materials, ame, f_e, rows), real, device)

    def readback(self) -> dict[str, NDArray[Any]]:
        """The uploaded arrays copied back to the host (their device dtype)."""
        return {k: getattr(self, k).numpy() for k in (*REAL_FIELDS, *INT_FIELDS)}

    def device_sha256(self) -> str:
        """``sha256`` recomputed from the bytes read back from the device (equals ``sha256``)."""
        return _hash_arrays(self.table_id, self.readback())


def _rows_digest(rows: Any) -> str:
    """sha256 of the per-material rows that are packed (so that a device never outlives a change
    of the
    rows it was built from)."""
    h = hashlib.sha256()
    for r in rows:
        for a in (r.target_index, r.sigma_mass_cm2_g, r.sigma_hat_window, r.sigma_hat_end,
                  r.cum_sigma_mass_cm2_g):  # fmt: skip
            arr = np.ascontiguousarray(a)
            h.update(f"{arr.dtype.str}|{arr.shape}".encode())
            h.update(arr.tobytes())
    return h.hexdigest()


def nuclear_device_key(
    table: NuclearTable,
    materials: tuple[Material, ...],
    f_e: float,
    real: Any,
    device: str,
    rows: Any = None,
) -> tuple[Any, ...]:
    """Cache key of a packed device.

    The table identity (``table_id`` and the npz sha256 of its sidecar), the material names in
    order, ``f_e``, the real dtype name, the device string and the digest of the packed material
    rows (default ``table.material_rows``)."""
    if rows is None:
        rows = [table.material_rows(m, f_e) for m in materials]
    return (
        table.table_id,
        table.info.get("npz_sha256"),
        tuple(m.name for m in materials),
        float(f_e),
        np.dtype(_np_dtype(real)).name,
        str(device),
        _rows_digest(rows),
    )


def cached_nuclear_device(
    table: NuclearTable,
    materials: tuple[Material, ...],
    *,
    real: Any,
    device: str,
    f_e: float = DEFAULT_F_E,
    rows: Any = None,
) -> NuclearDevice:
    """The packed :class:`NuclearDevice` of ``table`` / ``materials`` for ``real`` on ``device``.

    Packed and uploaded ONCE per process and key (``nuclear_device_key``) and reused by every later
    call (``rows``: the effective ``MaterialNuclear`` rows, which the key digests and the pack
    uses). The device arrays are read-only for the kernels (the nuclear sampler only reads
    ``grid``, ``lam``, ``edges``, ``rpre``, ``recoil``, ``tconst``, ``m_res``, ``sigma*``,
    ``cum_sigma`` and the material rows); the per-run scratch arrays (stack, event buffers, traces,
    counters) belong to ``NucData`` and are allocated per call by the driver. ``sha256`` is the
    unchanged identity of the uploaded bytes. A different key adds an entry; at most
    ``_MAX_CACHED_DEVICES`` are kept. Rebuilding the device for every ``run_range`` call (18000
    times per V7 shard) caused intermittent native crashes from allocation churn."""
    if rows is None:
        rows = [table.material_rows(m, f_e) for m in materials]
    key = nuclear_device_key(table, materials, f_e, real, device, rows)
    dev = _DEVICE_CACHE.get(key)
    if dev is None:
        dev = NuclearDevice.from_table(
            table, materials, real=real, device=device, f_e=f_e, rows=rows
        )
        while len(_DEVICE_CACHE) >= _MAX_CACHED_DEVICES:
            _DEVICE_CACHE.pop(next(iter(_DEVICE_CACHE)))
        _DEVICE_CACHE[key] = dev
    return dev


def clear_nuclear_device_cache() -> None:
    """Drop the per-process caches of the packed device and of the parsed AME2020 table."""
    _DEVICE_CACHE.clear()
    _AME_CACHE.clear()


__all__ = [
    "NuclearDevice",
    "NuclearHost",
    "cached_ame2020",
    "cached_nuclear_device",
    "clear_nuclear_device_cache",
    "host_sha256",
    "nuclear_device_key",
    "pack_nuclear",
]
