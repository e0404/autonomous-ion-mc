"""Device-side packing of the elastic table for the Warp kernels (decision 0041 slice C, V3-005C).

Mirrors :mod:`ionmc.transport.nuclear_device`: :func:`pack_elastic` flattens the table and the
per-material rows (``ElasticTable.material_rows``) into flat arrays (grid; per-target ``sigma`` and
mu_CM edges ``[target, node, n_q + 1]``; target masses ``[m_p, M_t...]``; per material the total
``Sigma_mass``, the window and end-of-range majorants and the cumulative target ``Sigma``).
:class:`ElasticDevice` uploads them as ``wp.array`` of the kernel precision (float32 or float64;
index rows int32) with a ``sha256`` identity of the uploaded bytes (:meth:`device_sha256` re-hashes
the bytes read back). :func:`cached_elastic_device` packs and uploads ONCE per process and key
(table
id and npz sha256, material names, ``f_e``, dtype, device and the digest of the packed rows),
keeps at
most four devices and never rebuilds per ``run_range`` call (allocation churn). No kernel uses these
arrays yet (the elastic channel lands in C4/C5).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import numpy as np
import warp as wp
from numpy.typing import NDArray

from ionmc.materials import Material
from ionmc.nuclear.elastic_tables import ElasticTable
from ionmc.nuclear.tables import DEFAULT_F_E
from ionmc.transport.nuclear_device import _np_dtype, _rows_digest

REAL_FIELDS = (
    "grid",
    "edges",
    "sigma_target",
    "target_mass",
    "e_min_shape",
    "e_min_pp",
    "sigma",
    "sigma_win",
    "sigma_end",
    "cum_sigma",
)
INT_FIELDS = ("mat_ntargets", "mat_target")
_MAX_CACHED_DEVICES = 4
_DEVICE_CACHE: dict[tuple[Any, ...], ElasticDevice] = {}


@dataclass(frozen=True)
class ElasticHost:
    """Packed float64/int64 arrays of a table and its materials."""

    table_id: str
    n_grid: int
    n_targets: int
    n_edges: int
    kmax: int
    n_materials: int
    arrays: dict[str, NDArray[Any]]


def pack_elastic(
    table: ElasticTable,
    materials: tuple[Material, ...],
    f_e: float = DEFAULT_F_E,
    rows: Any = None,
) -> ElasticHost:
    a = table.arrays
    grid = np.asarray(a["grid_e_mev"], dtype=np.float64)
    n = grid.size
    n_t, n_e = len(table.target_names), a["edges_mu"].shape[2]
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
        n_k[i], tgt[i, :k], cum[i, :k] = k, r.target_index, r.cum_sigma_mass_cm2_g
    arrays: dict[str, NDArray[Any]] = {
        "grid": grid,
        "edges": np.ascontiguousarray(a["edges_mu"], dtype=np.float64).ravel(),
        "sigma_target": np.ascontiguousarray(a["sigma_barn"], dtype=np.float64).ravel(),
        "target_mass": np.asarray(a["target_mass_mev"], dtype=np.float64),
        # per-target lower domain bound of the p + A shape (H-1 entry 0) and E_min,pp (size 1), so
        # that the transport can count elastic_below_domain and pp_below_domain
        "e_min_shape": np.concatenate(([0.0], np.asarray(a["target_e_min_mev"][1:], np.float64))),
        "e_min_pp": np.asarray(a["target_e_min_mev"][:1], dtype=np.float64),
        "sigma": sigma.ravel(), "sigma_win": win.ravel(), "sigma_end": end.ravel(),
        "cum_sigma": cum.ravel(), "mat_ntargets": n_k, "mat_target": tgt.ravel(),
    }  # fmt: skip
    if arrays["edges"].size != n_t * n * n_e or arrays["sigma_target"].size != n_t * n:
        raise ValueError("elastic table array sizes are inconsistent")
    return ElasticHost(table.table_id, n, n_t, n_e, kmax, len(rows), arrays)


def _hash_arrays(table_id: str, named: dict[str, NDArray[Any]]) -> str:
    h = hashlib.sha256()
    h.update(table_id.encode())
    for name in (*REAL_FIELDS, *INT_FIELDS):
        arr = np.ascontiguousarray(named[name])
        h.update(f"{name}|{arr.dtype.str}|{arr.shape}".encode())
        h.update(arr.tobytes())
    return h.hexdigest()


def _cast(host: ElasticHost, np_real: Any) -> dict[str, NDArray[Any]]:
    cast: dict[str, NDArray[Any]] = {
        k: np.ascontiguousarray(host.arrays[k], dtype=np_real) for k in REAL_FIELDS
    }
    cast.update({k: np.ascontiguousarray(host.arrays[k], dtype=np.int32) for k in INT_FIELDS})
    return cast


def host_sha256(host: ElasticHost, precision: str) -> str:
    """``ElasticDevice.sha256`` of ``host`` in ``precision`` ("float32" / "float64") without a
    device."""
    return _hash_arrays(
        host.table_id, _cast(host, {"float32": np.float32, "float64": np.float64}[precision])
    )


class ElasticDevice:
    """The packed elastic arrays on one device in one precision (attributes named as
    :data:`REAL_FIELDS` / :data:`INT_FIELDS` plus ``n_grid``, ``n_targets``, ``n_edges``, ``kmax``,
    ``n_materials``, ``real``, ``device``, ``table_id`` and ``sha256``)."""

    def __init__(self, host: ElasticHost, real: Any, device: str) -> None:
        cast = _cast(host, _np_dtype(real))
        self.real, self.device, self.table_id = real, device, host.table_id
        self.n_grid, self.n_targets, self.n_edges = host.n_grid, host.n_targets, host.n_edges
        self.kmax, self.n_materials = host.kmax, host.n_materials
        self.sha256 = _hash_arrays(host.table_id, cast)
        for k in REAL_FIELDS:
            setattr(self, k, wp.array(cast[k], dtype=real, device=device))
        for k in INT_FIELDS:
            setattr(self, k, wp.array(cast[k], dtype=wp.int32, device=device))

    def readback(self) -> dict[str, NDArray[Any]]:
        return {k: getattr(self, k).numpy() for k in (*REAL_FIELDS, *INT_FIELDS)}

    def device_sha256(self) -> str:
        return _hash_arrays(self.table_id, self.readback())


def elastic_device_key(
    table: ElasticTable,
    materials: tuple[Material, ...],
    f_e: float,
    real: Any,
    device: str,
    rows: Any = None,
) -> tuple[Any, ...]:
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


def cached_elastic_device(
    table: ElasticTable,
    materials: tuple[Material, ...],
    *,
    real: Any,
    device: str,
    f_e: float = DEFAULT_F_E,
    rows: Any = None,
) -> ElasticDevice:
    """The packed device of ``table`` / ``materials``, packed and uploaded once per process and
    key."""
    if rows is None:
        rows = [table.material_rows(m, f_e) for m in materials]
    key = elastic_device_key(table, materials, f_e, real, device, rows)
    dev = _DEVICE_CACHE.get(key)
    if dev is None:
        dev = ElasticDevice(pack_elastic(table, materials, f_e, rows), real, device)
        while len(_DEVICE_CACHE) >= _MAX_CACHED_DEVICES:
            _DEVICE_CACHE.pop(next(iter(_DEVICE_CACHE)))
        _DEVICE_CACHE[key] = dev
    return dev


def clear_elastic_device_cache() -> None:
    _DEVICE_CACHE.clear()


__all__ = [
    "ElasticDevice", "ElasticHost", "cached_elastic_device", "clear_elastic_device_cache",
    "elastic_device_key", "host_sha256", "pack_elastic",
]  # fmt: skip
