"""V3-005C C4: the hadronic elastic channel in the Warp nuclear kernel (warp-cpu; CUDA variants are
marked ``cuda`` and skipped without a device).

Rows: P5-ext (event parity on recorded inputs: the ``@wp.func`` twins of
:mod:`ionmc.physics.elastic` called from a Warp kernel replay ``meta["elastic_trace"]`` of the
python reference; float64 bit-equal, float32 within the maxima reported), V8-LV (trajectory parity
python vs warp-cpu float64 with the channel on: identical events, species, genealogy, counters;
continuous columns <= 1e-10), counters parity, ``elastic=False`` / Sigma_el = 0 bit identity of
the Warp path, P7 (chi2 of the Warp samplers), the forced-counter C1-ext cases and the device
cache / fail-closed checks. Data-backed (``IONMC_CACHE_DIR``).

float32 note: in the float32 variant the target choice and the CM-cosine sampling run in float32
on the float32 table; the two-body kinematics and the outgoing states are evaluated in float64
from the (float32) cosine, because energies and the per-event ledger are float64 in every variant
(as the non-elastic ledger of the nuclear kernel). The float32 event parity therefore measures the
float32 sampling error of mu alone, reported separately from the outgoing-energy error of the
recorded cosine.

P5-ext coverage by twin (all on warp-cpu against the python twin / reference):

* ``channel_is_elastic``: ``test_channel_choice_parity_synthetic_*`` (elastic / non-elastic, the
  boundary ``u Sigma_tot == Sigma_el`` exactly and one ulp either side in the precision under test,
  Sigma_tot = 0, Sigma_el = 0, Sigma_el = Sigma_tot, u = 0, u -> 1; float64 and float32 bit-equal
  to the python twin evaluated in the same precision) and
  ``test_channel_choice_parity_recorded_candidates`` (candidate uniforms u3 and the sigmas
  re-derived from the recorded ``elastic_trace``: every recorded elastic event is chosen elastic
  by the Warp twin and every other block of the same Philox streams gives the same decision as the
  twin; float64 and float32);
* ``faster_is_recoil``: ``test_faster_is_recoil_parity_*`` (synthetic ties / ordering and the
  recorded (T_primary, T_other) of the trace, both orders);
* ``sample_mu_edges``: P7 ``test_p7_warp_sampler_chi2`` plus the replay of the recorded events;
* ``elastic_sample``, ``elastic_two_body`` (inside ``elastic_outgoing``), ``frame_dir`` and
  ``elastic_outgoing``: the recorded-input replay ``test_event_parity_recorded_inputs_*`` (the
  float64 replay is bit-equal, including the p-p swap rule and both directions);
* run level: ``test_trajectory_parity_python_vs_warp_cpu_float64_elastic_on`` (identical event
  records, counters and tallies, hence the same number of elastic events as the python reference;
  ``test_elastic_event_counts_equal_python`` states the equality of the event counters).

Recovery of the candidate inputs from ``elastic_trace`` (the columns are not extended): the trace
carries the history, the genealogy id, ``u_target`` (slot 2 of the candidate block ``nc``) and the
event-block uniforms (``u_mu``, ``u_phi`` = slots 0, 1 of block ``nc + 1``); the candidate block
index itself is not recorded, so the test regenerates the Philox stream of (history, genealogy id)
and finds the block ``nc`` with ``u[nc][2] == u_target`` and ``u[nc + 1][:2] == (u_mu, u_phi)``
(exact float equality). ``u3 = u[nc][3]``; ``Sigma_el(T_in)`` and ``Sigma_nonel(T_in)`` come from
the material rows through the same lookup as the reference (single-material configuration)."""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np
import pytest
import warp as wp

import ionmc.transport.kernels_nuclear as kn
from ionmc._wpfunc import python_twin
from ionmc.config import DiagnosticsOptions
from ionmc.nuclear.elastic_tables import ElasticTable, chi2_equiprobable
from ionmc.physics.elastic import make_elastic
from ionmc.physics.nuclear import make_nuclear
from ionmc.physics.projectiles import PROTON
from ionmc.rng.philox import PURPOSE_NUCLEAR, draw_block, key_from_seed, u01_py
from ionmc.simulation import Simulation
from ionmc.transport import elastic_device as ed
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.reference import run_reference_range
from ionmc.transport.tally import ELASTIC_COUNTER_NAMES
from ionmc.transport.warp_driver import run_warp_range
from tests.ionmc.test_elastic_transport import _cfg, _pin, _scale_elastic
from tests.ionmc.test_nuclear_transport_loop import _table_id
from tests.ionmc.test_nuclear_warp import _compare, _sums

M_P = PROTON.mass_mev
NU = python_twin(make_nuclear)
N_TRAJ = 32
TOL = 1e-10


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


@pytest.fixture(scope="module")
def table() -> ElasticTable:
    return ElasticTable.load(None, _pin())


def _warp_cfg(tid: str, energy: float, n: int, precision: str, *, trace: int = 0, **phys: Any):
    allow = bool(phys.pop("allow_invalid", False))
    cfg = _cfg(tid, energy, n, trace=trace, **phys)
    return replace(
        cfg,
        run=replace(cfg.run, backend="warp-cpu", precision=precision, allow_invalid_result=allow),
        diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=trace),
    )


@pytest.fixture(scope="module")
def pair(tid: str) -> tuple[Any, Any, Any]:
    """Python reference and warp-cpu float64 runs of the same dense configuration (Sigma_el x 30),
    with traces: ``(effective_python, python_part, warp_part)``."""
    n = N_TRAJ
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        cfg = _cfg(tid, 150.0, n, trace=n)
        cfg = replace(
            cfg, diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=n)
        )
        ep = Simulation(cfg).effective
        py = run_reference_range(ep, 0, n)
        ew = Simulation(_warp_cfg(tid, 150.0, n, "float64", trace=n)).effective
        wr = run_warp_range(ew, 0, n, "cpu")
    return ep, py, wr


# ---- V8-LV trajectory parity, counters, balance -------------------------------------------------
def test_trajectory_parity_python_vs_warp_cpu_float64_elastic_on(
    pair: tuple[Any, Any, Any],
) -> None:
    """Identical events (step trace, end codes, voxels), species and genealogy of the non-elastic
    secondaries, tallies including the elastic columns and counters; continuous columns <= 1e-10;
    the elastic domain counters are part of the counter parity."""
    _, py, wr = pair
    m = _compare(py, wr)
    sp, cp = _sums(py)
    sw, cw = _sums(wr)
    n_el = len(ELASTIC_COUNTER_NAMES)
    assert cp == cw and len(cp) >= n_el
    assert cp[-1] > 0 and cp[-2] > 0  # both domain counters fired (diagnostics only)
    assert sp[-1] > 0 and sp[-2] > 0 and sp[-3] > 0  # recoil energy, p-p and p + A events
    print("C4 trajectory parity maxima", m, "elastic tallies", sp[-3:], "counters", cp)
    assert max(m.values()) <= TOL, m
    assert all(a == b for a, b in zip(py.edep[0].ravel(), wr.edep[0].ravel(), strict=True))


# ---- P5-ext: the Warp twins on the recorded inputs ----------------------------------------------
def _event_kernel(real: Any) -> Any:
    name = "float32" if real is wp.float32 else "float64"
    EL = make_elastic(real)
    ELD = make_elastic(wp.float64)
    F = make_transport_funcs(real)
    v3 = F.vec3
    two_pi = wp.constant(wp.float64(2.0 * math.pi))

    def ev(
        rows: wp.array2d(dtype=wp.float64),
        grid: wp.array(dtype=real),
        sigma: wp.array(dtype=real),
        cum: wp.array(dtype=real),
        edges: wp.array(dtype=real),
        mass: wp.array(dtype=wp.float64),
        mat_target: wp.array(dtype=wp.int32),
        mat_nt: wp.array(dtype=wp.int32),
        kmax: int,
        n_grid: int,
        n_e: int,
        m1: wp.float64,
        out: wp.array2d(dtype=wp.float64),
    ):
        i = wp.tid()
        t1 = rows[i, 4]
        d = v3(real(rows[i, 8]), real(rows[i, 9]), real(rows[i, 10]))
        e1, e2 = F.orthonormal_basis(d)
        tgt, mu = EL.elastic_sample(
            real(rows[i, 5]),
            real(rows[i, 6]),
            real(t1),
            0,
            int(mat_nt[0]),
            kmax,
            n_grid,
            n_e,
            grid,
            sigma,
            cum,
            edges,
            mat_target,
        )
        pp = 0
        if tgt == 0:
            pp = 1
        m2 = wp.float64(mass[tgt])
        phi = two_pi * rows[i, 7]
        ta, ax, ay, az, tb, bx, by, bz = ELD.elastic_outgoing(
            t1,
            m1,
            m2,
            wp.float64(mu),
            phi,
            pp,
            wp.float64(e1[0]),
            wp.float64(e1[1]),
            wp.float64(e1[2]),
            wp.float64(e2[0]),
            wp.float64(e2[1]),
            wp.float64(e2[2]),
            wp.float64(d[0]),
            wp.float64(d[1]),
            wp.float64(d[2]),
        )
        out[i, 0] = wp.float64(tgt)
        out[i, 1] = wp.float64(mu)
        out[i, 2] = phi
        out[i, 3] = ta
        out[i, 4] = ax
        out[i, 5] = ay
        out[i, 6] = az
        out[i, 7] = tb
        out[i, 8] = bx
        out[i, 9] = by
        out[i, 10] = bz
        # outgoing states from the RECORDED cosine (isolates the kinematics from the sampling)
        ta2, ax2, ay2, az2, tb2, bx2, by2, bz2 = ELD.elastic_outgoing(
            t1,
            m1,
            m2,
            rows[i, 11],
            phi,
            pp,
            wp.float64(e1[0]),
            wp.float64(e1[1]),
            wp.float64(e1[2]),
            wp.float64(e2[0]),
            wp.float64(e2[1]),
            wp.float64(e2[2]),
            wp.float64(d[0]),
            wp.float64(d[1]),
            wp.float64(d[2]),
        )
        out[i, 11] = ta2
        out[i, 12] = ax2
        out[i, 13] = ay2
        out[i, 14] = az2
        out[i, 15] = tb2
        out[i, 16] = bx2
        out[i, 17] = by2
        out[i, 18] = bz2

    ev.__name__ = f"el_event_parity_{name}"
    ev.__qualname__ = ev.__name__
    return wp.kernel(enable_backward=False, module="unique")(ev)


def _replay(eff: Any, trace: np.ndarray, real: Any) -> np.ndarray:
    el = eff.nuclear.elastic
    host = ed.pack_elastic(el.table, eff.geometry.materials, rows=el.rows)
    dev = ed.ElasticDevice(host, real, "cpu")
    mass = dev.el_mass  # the cached float64 masses the transport kernel consumes
    rows = wp.array(np.ascontiguousarray(trace), dtype=wp.float64, device="cpu")
    out = wp.zeros((trace.shape[0], 19), dtype=wp.float64, device="cpu")
    wp.launch(
        _event_kernel(real), dim=trace.shape[0],
        inputs=[rows, dev.grid, dev.sigma, dev.cum_sigma, dev.edges, mass,
                dev.mat_target, dev.mat_ntargets, dev.kmax, dev.n_grid, dev.n_edges,
                wp.float64(M_P), out],
        device="cpu",
    )  # fmt: skip
    return out.numpy()


def test_event_parity_recorded_inputs_float64_bit_equal(pair: tuple[Any, Any, Any]) -> None:
    """P5-ext: the Warp twins of ``elastic_sample`` / ``elastic_outgoing`` (float64) reproduce the
    recorded python events (p-p and p + A) bit for bit: target, mu_CM, phi, outgoing energies and
    directions of both bodies."""
    ep, py, _ = pair
    tr = py.meta["elastic_trace"]
    out = _replay(ep, tr, wp.float64)
    pp = tr[:, 3] == 0
    assert pp.sum() > 0 and (~pp).sum() > 0
    np.testing.assert_array_equal(out[:, 0], tr[:, 3])
    np.testing.assert_array_equal(out[:, 1:3], tr[:, 11:13])
    np.testing.assert_array_equal(out[:, 3:11], tr[:, 13:21])
    np.testing.assert_array_equal(out[:, 11:19], tr[:, 13:21])
    print(
        "P5-ext float64: events", len(tr), "pp", int(pp.sum()), "pA", int((~pp).sum()), "bit-equal"
    )


def test_event_parity_recorded_inputs_float32_reported(pair: tuple[Any, Any, Any]) -> None:
    """P5-ext, float32 (see the module docstring): discrete outputs (target) equal; the maxima are
    reported separately for the sampled mu_CM and for the outgoing energies / directions, the
    latter from the recorded cosine (pure float64 kinematics: equal to the float64 replay) and
    from the float32-sampled cosine."""
    ep, py, _ = pair
    tr = py.meta["elastic_trace"]
    out = _replay(ep, tr, wp.float32)
    np.testing.assert_array_equal(out[:, 0], tr[:, 3])
    d_mu = float(np.max(np.abs(out[:, 1] - tr[:, 11])))
    e_rec = float(np.max(np.abs(out[:, 11] - tr[:, 13]) / tr[:, 4]))
    e_rec_o = float(np.max(np.abs(out[:, 15] - tr[:, 17]) / tr[:, 4]))
    e_smp = float(np.max(np.abs(out[:, 3] - tr[:, 13]) / tr[:, 4]))
    e_smp_o = float(np.max(np.abs(out[:, 7] - tr[:, 17]) / tr[:, 4]))
    print(f"P5-ext float32: max|dmu|={d_mu:.3e}; energies/T_in: recorded-mu primary {e_rec:.3e} "
          f"other {e_rec_o:.3e}; sampled-mu primary {e_smp:.3e} other {e_smp_o:.3e}")  # fmt: skip
    assert d_mu <= 2e-6 and e_rec <= 1e-12 and e_rec_o <= 1e-12
    assert e_smp <= 1e-5 and e_smp_o <= 1e-5
    # the float64 kinematics of the float32 variant conserve energy to double rounding
    assert np.max(np.abs(tr[:, 4] - out[:, 3] - out[:, 7])) <= 1e-9


# ---- P5-ext: channel choice and the faster-proton rule -----------------------------------------
EL_PY = python_twin(make_elastic)
_REAL = {"float64": (wp.float64, np.float64), "float32": (wp.float32, np.float32)}


def _channel_kernel(real: Any) -> Any:
    EL = make_elastic(real)

    def ch(
        u: wp.array(dtype=real),
        s_el: wp.array(dtype=real),
        s_tot: wp.array(dtype=real),
        t1: wp.array(dtype=real),
        t2: wp.array(dtype=real),
        out: wp.array(dtype=wp.int32),
        rec: wp.array(dtype=wp.int32),
    ):
        i = wp.tid()
        out[i] = wp.int32(EL.channel_is_elastic(real(u[i]), s_el[i], s_tot[i]))
        rec[i] = wp.int32(EL.faster_is_recoil(t1[i], t2[i]))

    ch.__name__ = f"el_channel_parity_{'f32' if real is wp.float32 else 'f64'}"
    ch.__qualname__ = ch.__name__
    return wp.kernel(enable_backward=False, module="unique")(ch)


def _warp_channel(precision: str, u, s_el, s_tot, t1=None, t2=None):
    real, npr = _REAL[precision]
    u, s_el, s_tot = (np.asarray(a, dtype=npr) for a in (u, s_el, s_tot))
    t1 = np.zeros_like(u) if t1 is None else np.asarray(t1, dtype=npr)
    t2 = np.zeros_like(u) if t2 is None else np.asarray(t2, dtype=npr)
    out = wp.zeros(u.size, dtype=wp.int32, device="cpu")
    rec = wp.zeros(u.size, dtype=wp.int32, device="cpu")
    wp.launch(
        _channel_kernel(real), dim=u.size, device="cpu",
        inputs=[wp.array(a, dtype=real, device="cpu") for a in (u, s_el, s_tot, t1, t2)]
        + [out, rec],
    )  # fmt: skip
    return out.numpy(), rec.numpy()


def _twin_channel(precision: str, u, s_el, s_tot) -> np.ndarray:
    npr = _REAL[precision][1]
    return np.array(
        [
            EL_PY.channel_is_elastic(npr(a), npr(b), npr(c))
            for a, b, c in zip(u, s_el, s_tot, strict=True)
        ],
        dtype=np.int32,
    )


def _synthetic_channel_cases(npr: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(u, Sigma_el, Sigma_tot)`` in the precision ``npr``: ordinary elastic / non-elastic
    decisions, the boundary ``u Sigma_tot == Sigma_el`` (product rounded in the precision) and
    its two neighbours, Sigma_tot = 0, Sigma_el = 0, Sigma_el = Sigma_tot, u = 0, u -> 1."""
    rng = np.random.default_rng(5)
    rows: list[tuple[float, float, float]] = []
    for _ in range(64):  # ordinary
        tot = float(rng.uniform(0.01, 5.0))
        rows.append((float(rng.uniform()), float(rng.uniform()) * tot, tot))
    for _ in range(64):  # boundary: Sigma_el = fl(u Sigma_tot) and the neighbours either side
        u, tot = npr(rng.uniform(0.05, 0.95)), npr(rng.uniform(0.1, 4.0))
        p = npr(u * tot)
        for sel in (p, np.nextafter(p, npr(0)), np.nextafter(p, npr(np.inf))):
            rows.append((float(u), float(sel), float(tot)))
    one_minus = float(np.nextafter(npr(1), npr(0)))
    rows += [
        (0.5, 0.5, 0.0), (0.0, 0.0, 0.0), (0.3, 0.0, 0.0),  # Sigma_tot = 0 (never elastic)
        (0.4, 0.0, 2.0), (0.0, 0.0, 2.0),  # Sigma_el = 0
        (0.99, 2.0, 2.0), (one_minus, 2.0, 2.0), (0.0, 2.0, 2.0),  # Sigma_el = Sigma_tot
        (0.0, 1e-30, 2.0), (0.0, 0.5, 2.0),  # u = 0: elastic iff Sigma_el > 0
        (one_minus, 1.0, 1.0), (one_minus, np.nextafter(npr(1), npr(2)), 1.0),  # u -> 1
        (one_minus, np.nextafter(npr(1), npr(0)), 1.0),
        (0.5, 1.0, 2.0), (0.5, np.nextafter(npr(1), npr(2)), 2.0),  # exactly equal, tie side
    ]  # fmt: skip
    a = np.array(rows, dtype=npr)
    return a[:, 0], a[:, 1], a[:, 2]


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_channel_choice_parity_synthetic(precision: str) -> None:
    """P5-ext ``channel_is_elastic``: Warp (warp-cpu) bit-equal to the python twin in the same
    precision on synthetic candidates incl. the exact boundary, Sigma_tot = 0, Sigma_el = 0,
    Sigma_el = Sigma_tot, u = 0, u -> 1; both decisions occur and the boundary is strict."""
    u, se, st = _synthetic_channel_cases(_REAL[precision][1])
    got, _ = _warp_channel(precision, u, se, st)
    want = _twin_channel(precision, u, se, st)
    np.testing.assert_array_equal(got, want)
    assert 0 < int(got.sum()) < got.size
    npr = _REAL[precision][1]
    exact = (npr(0.5), npr(1.0), npr(2.0))  # u Sigma_tot == Sigma_el exactly: not elastic
    g, _ = _warp_channel(precision, *[[v] for v in exact])
    assert int(g[0]) == 0
    below = [[npr(0.5)], [np.nextafter(npr(1), npr(2))], [npr(2.0)]]
    assert int(_warp_channel(precision, *below)[0][0]) == 1
    assert int(_warp_channel(precision, [0.7], [0.0], [0.0])[0][0]) == 0  # Sigma_tot = 0


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_faster_is_recoil_parity_synthetic(precision: str) -> None:
    """P5-ext ``faster_is_recoil``: strict ``T2 > T1`` (ties keep the first particle); bit-equal to
    the python twin, in the precision under test (a float32 pair that rounds to a tie is a tie)."""
    npr = _REAL[precision][1]
    rng = np.random.default_rng(9)
    t1 = rng.uniform(0.1, 150.0, 200).astype(npr)
    t2 = rng.uniform(0.1, 150.0, 200).astype(npr)
    t2[:20] = t1[:20]  # exact ties
    t2[20:40] = np.nextafter(t1[20:40], npr(1e9))  # one ulp above
    t2[40:60] = np.nextafter(t1[40:60], npr(0))  # one ulp below
    _, got = _warp_channel(precision, t1, t1, t1, t1, t2)
    want = np.array(
        [EL_PY.faster_is_recoil(npr(a), npr(b)) for a, b in zip(t1, t2, strict=True)],
        dtype=np.int32,
    )
    np.testing.assert_array_equal(got, want)
    assert np.all(got[:20] == 0) and np.all(got[20:40] == 1) and np.all(got[40:60] == 0)


def _row_at(rows: Any, e: float) -> float:
    """``Sigma_mass(E)`` of a material's rows as the reference lookup (step cell by the shared
    ``grid_locate``, lin-lin in E, clamped)."""
    g = rows.grid_e_mev
    k = int(NU.grid_locate(e, g, g.size))
    t = min(max((e - g[k]) / (g[k + 1] - g[k]), 0.0), 1.0)
    row = rows.sigma_mass_cm2_g
    return float((1.0 - t) * row[k] + t * row[k + 1])


def _recorded_candidates(
    ep: Any, tr: np.ndarray, extra_blocks: int = 24
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Candidate inputs ``(u3, Sigma_el, Sigma_tot, recorded decision)`` re-derived from the
    elastic trace (module docstring), for every recorded elastic event (decision 1) and for
    ``extra_blocks`` further blocks of the same stream at the same energy (decision = the python
    twin's, evaluated by the caller)."""
    assert len(ep.geometry.materials) == 1, "single-material configuration expected"
    key = key_from_seed(ep.requested.run.seed)
    nuc_rows, el_rows = ep.nuclear.rows[0], ep.nuclear.elastic.rows[0]
    u3s, sels, stots, dec = [], [], [], []
    for r in tr:
        h, gid, t1 = int(r[0]), int(r[1]), float(r[4])
        s_el, s_n = _row_at(el_rows, t1), _row_at(nuc_rows, t1)
        blocks: list[list[float]] = []
        nc = -1
        for b in range(20000):
            words = draw_block(key, h, gid, b, PURPOSE_NUCLEAR)
            blocks.append([u01_py(x, "float64") for x in words])
            if b >= 1 and blocks[b - 1][2] == r[5] and blocks[b][:2] == [r[6], r[7]]:
                nc = b - 1
                break
        assert nc >= 0, "candidate block not recovered from the trace"
        for b, d in [(nc, 1)] + [(b, -1) for b in range(min(extra_blocks, nc))]:
            u3s.append(blocks[b][3])
            sels.append(s_el)
            stots.append(s_n + s_el)
            dec.append(d)
    return np.array(u3s), np.array(sels), np.array(stots), np.array(dec)


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_channel_choice_parity_recorded_candidates(
    pair: tuple[Any, Any, Any], precision: str
) -> None:
    """P5-ext ``channel_is_elastic`` on recorded candidate inputs of the python run: every recorded
    elastic event is elastic in the python reference's float64 evaluation and in the Warp float64
    twin; in either precision the Warp decision equals the python twin's on all recovered
    candidates, including the other (mostly non-elastic) blocks of the same Philox streams."""
    ep, py, _ = pair
    tr = py.meta["elastic_trace"]
    assert len(tr) > 0
    u3, se, st, dec = _recorded_candidates(ep, tr)
    got, _ = _warp_channel(precision, u3, se, st)
    want = _twin_channel(precision, u3, se, st)
    np.testing.assert_array_equal(got, want)
    rec = dec == 1
    if precision == "float64":
        assert np.all(got[rec] == 1)  # the recorded events were chosen elastic
    assert np.mean(got[rec]) >= 0.999  # float32 may flip a measure-zero boundary case
    assert int(got[~rec].sum()) > 0  # the same streams also give elastic decisions elsewhere
    n_other = int((~rec).sum())
    print(f"P5-ext channel {precision}: {int(rec.sum())} recorded elastic, {n_other} other "
          f"candidate blocks (elastic {int(got[~rec].sum())})")  # fmt: skip


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_faster_is_recoil_parity_recorded(pair: tuple[Any, Any, Any], precision: str) -> None:
    """P5-ext ``faster_is_recoil`` on the recorded outgoing energies (primary, other) of the
    elastic trace, both argument orders, equals the python twin."""
    _, py, _ = pair
    tr = py.meta["elastic_trace"]
    for a, b in ((tr[:, 13], tr[:, 17]), (tr[:, 17], tr[:, 13])):
        _, got = _warp_channel(precision, a, a, a, a, b)
        npr = _REAL[precision][1]
        want = np.array(
            [EL_PY.faster_is_recoil(npr(x), npr(y)) for x, y in zip(a, b, strict=True)],
            dtype=np.int32,
        )
        np.testing.assert_array_equal(got, want)


def test_elastic_event_counts_equal_python(pair: tuple[Any, Any, Any]) -> None:
    """Run level: the elastic event counters of warp-cpu float64 equal the python reference's
    (p-p and p + A events), and the recorded event count equals their sum."""
    _, py, wr = pair
    _, cp = _sums(py)
    _, cw = _sums(wr)
    sp, _ = _sums(py)
    sw, _ = _sums(wr)
    assert sp[-2] == sw[-2] and sp[-1] == sw[-1] and cp == cw
    assert int(sp[-1] + sp[-2]) == len(py.meta["elastic_trace"])


# ---- P7: Warp samplers ------------------------------------------------------------------------
def _mu_kernel(real: Any) -> Any:
    EL = make_elastic(real)

    def mu_k(
        u: wp.array(dtype=real),
        edges: wp.array(dtype=real),
        base_lo: int,
        base_hi: int,
        t: real,
        n_q: int,
        out: wp.array(dtype=wp.float64),
    ):
        i = wp.tid()
        uu = real(u[i])
        out[i] = wp.float64(EL.sample_mu_edges(uu, edges, base_lo, base_hi, t, n_q))

    mu_k.__name__ = f"el_mu_sampler_{'f32' if real is wp.float32 else 'f64'}"
    mu_k.__qualname__ = mu_k.__name__
    return wp.kernel(enable_backward=False, module="unique")(mu_k)


@pytest.mark.parametrize("real", [wp.float64, wp.float32])
@pytest.mark.parametrize("target,name", [(0, "H-1"), (3, "O-16")])
@pytest.mark.parametrize("energy", [20.0, 100.0, 200.0])
def test_p7_warp_sampler_chi2(
    table: ElasticTable, real: Any, target: int, name: str, energy: float
) -> None:
    assert table.target_names[target] == name
    n = 20000
    u = np.random.default_rng(20515 + target * 1000 + int(energy)).random(n)
    a = table.arrays
    grid = a["grid_e_mev"]
    k = int(NU.grid_locate(energy, grid, grid.size))
    t = min(max((energy - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
    n_g, n_e = a["edges_mu"].shape[1:]
    np_real = np.float64 if real is wp.float64 else np.float32
    u_d = wp.array(u.astype(np_real), dtype=real, device="cpu")
    ed_d = wp.array(a["edges_mu"].reshape(-1).astype(np_real), dtype=real, device="cpu")
    out = wp.zeros(n, dtype=wp.float64, device="cpu")
    wp.launch(
        _mu_kernel(real), dim=n,
        inputs=[u_d, ed_d, (target * n_g + k) * n_e, (target * n_g + k + 1) * n_e, real(t), n_e - 1,
                out],
        device="cpu",
    )  # fmt: skip
    mu = out.numpy()
    if real is wp.float64:
        np.testing.assert_array_equal(mu, table.sample_mu(target, energy, u))
    chi2, p, dof = chi2_equiprobable(mu, table.edges_at(target, energy), 64)
    print(f"P7 warp-cpu {np_real.__name__} {name} {energy} MeV chi2 {chi2:.1f} dof {dof} p {p:.4f}")
    assert p > 0.001


# ---- bit identity of the Warp paths -------------------------------------------------------------
@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_sigma_el_zero_warp_is_bit_identical_to_elastic_false(
    tid: str, precision: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The compiled elastic variant with Sigma_el = 0 (additive majorant, u3 draw, domain bounds)
    equals the elastic=False kernel in every non-elastic column, counter and grid, in the same
    process; the elastic columns are zero except the diagnostics domain counters."""
    n = 24
    off = run_warp_range(
        Simulation(_warp_cfg(tid, 100.0, n, precision, elastic=False)).effective, 0, n, "cpu"
    )
    with monkeypatch.context() as mp:
        _scale_elastic(mp, 0.0)
        on = run_warp_range(Simulation(_warp_cfg(tid, 100.0, n, precision)).effective, 0, n, "cpu")
    n_t, n_c = len(off.tally_components), len(off.counter_sums)
    assert [math.fsum(c) for c in on.tally_components[:n_t]] == [
        math.fsum(c) for c in off.tally_components
    ]
    assert on.counter_sums[:n_c] == off.counter_sums
    for a, b in zip(on.edep, off.edep, strict=True):
        np.testing.assert_array_equal(a, b)
    assert all(math.fsum(c) == 0.0 for c in on.tally_components[n_t:])


def test_elastic_false_warp_matches_python_elastic_false(tid: str) -> None:
    """Same-process anchor of the pre-D2 path: the elastic=False warp-cpu float64 run follows the
    python elastic=False run (V8 parity, <= 1e-10) after the kernel change."""
    n = 8
    cfg = _warp_cfg(tid, 150.0, n, "float64", trace=n, elastic=False)
    pcfg = replace(cfg, run=replace(cfg.run, backend="python"))
    py = run_reference_range(Simulation(pcfg).effective, 0, n)
    wr = run_warp_range(Simulation(cfg).effective, 0, n, "cpu")
    m = _compare(py, wr)
    assert max(m.values()) <= TOL, m


# ---- config, device cache, fail-closed ----------------------------------------------------------
def test_warp_backends_accept_elastic_and_pack_device(tid: str) -> None:
    eff = Simulation(_warp_cfg(tid, 100.0, 4, "float32")).effective
    el = eff.nuclear.elastic
    assert el is not None
    s = eff.summary()["elastic"]
    host = ed.pack_elastic(el.table, eff.geometry.materials, rows=el.rows)
    assert s["device_sha256"] == ed.host_sha256(host, "float32")
    assert (
        "device_sha256"
        not in Simulation(replace(_cfg(tid, 100.0, 4))).effective.summary()["elastic"]
    )  # python backend: no uploaded arrays


def test_device_domain_bounds_equal_elasticsetup(tid: str) -> None:
    """The kernel's per-material e_min_pa / e_min_pp, derived on device from the packed
    ``e_min_shape`` / ``e_min_pp`` / ``mat_target`` fields, equal ``ElasticSetup``."""
    eff = Simulation(_warp_cfg(tid, 100.0, 4, "float64")).effective
    el = eff.nuclear.elastic
    host = ed.pack_elastic(el.table, eff.geometry.materials, rows=el.rows)
    a = host.arrays
    for m in range(host.n_materials):
        pa, pp = 0.0, 0.0
        for j in range(int(a["mat_ntargets"][m])):
            t = int(a["mat_target"][m * host.kmax + j])
            if t == 0:
                pp = float(a["e_min_pp"][0])
            else:
                pa = max(pa, float(a["e_min_shape"][t]))
        assert pa == el.e_min_pa[m] and pp == el.e_min_pp[m]


def test_device_cached_once_andverified(tid: str) -> None:
    ed.clear_elastic_device_cache()
    eff = Simulation(_warp_cfg(tid, 100.0, 4, "float64")).effective
    from ionmc.transport.warp_driver import _elastic_device

    a = _elastic_device(eff, wp.float64, "cpu")
    assert _elastic_device(eff, wp.float64, "cpu") is a and a.verified
    eff32 = Simulation(_warp_cfg(tid, 100.0, 4, "float32")).effective
    assert _elastic_device(eff32, wp.float32, "cpu") is not a
    with pytest.raises(ValueError, match="precision"):  # precision of the device differs from eff
        _elastic_device(eff, wp.float32, "cpu")
    # a corrupted upload fails closed
    ed.clear_elastic_device_cache()
    b = _elastic_device(eff, wp.float64, "cpu")
    b.verified = False
    b.sigma.fill_(1.0)
    with pytest.raises(ValueError, match="sha256"):
        _elastic_device(eff, wp.float64, "cpu")
    ed.clear_elastic_device_cache()


def test_elastic_only_and_balance_warp_cpu(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    _scale_elastic(monkeypatch, 30.0)
    n = 24
    res = Simulation(
        replace(_warp_cfg(tid, 150.0, n, "float64"), diagnostics=DiagnosticsOptions())
    ).run()
    b = res.energy_balance
    assert res.valid and b.relative_residual <= 1e-12 and b.grid_relative_residual(0) <= 1e-12
    assert b.elastic_mev["elastic_pp_events"] > 0 and b.elastic_mev["elastic_pa_events"] > 0
    only = Simulation(
        replace(
            _warp_cfg(tid, 100.0, 16, "float64", elastic_only=True),
            diagnostics=DiagnosticsOptions(),
        )
    ).run()
    assert only.valid and only.energy_balance.relative_residual <= 1e-12
    assert sum(r["events"] for r in only.diagnostics.get("nuclear", {}).values()) == 0
    assert only.energy_balance.nuclear_mev["nuclear_local"] == 0.0


def test_float32_elastic_closes_balance(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    _scale_elastic(monkeypatch, 30.0)
    res = Simulation(
        replace(_warp_cfg(tid, 150.0, 24, "float32"), diagnostics=DiagnosticsOptions())
    ).run()
    assert res.valid and res.counters.nuclear_conservation == 0
    assert res.energy_balance.elastic_mev["elastic_pa_events"] > 0
    assert res.energy_balance.relative_residual <= 1e-9


def test_device_mass_array_is_cached_hashed_and_verified(tid: str) -> None:
    """V3-005C C4 review: the float64 target masses ``el_mass`` the kernel consumes are a member of
    the cached device (uploaded once with it), part of the recorded sha256 and of the readback
    verification: a corrupted device-resident mass array is detected; the hash changes with a
    mass value; the f32 / f64 hashes differ only through the precision-cast arrays."""
    from ionmc.transport.warp_driver import _elastic_device

    ed.clear_elastic_device_cache()
    eff = Simulation(_warp_cfg(tid, 100.0, 4, "float64")).effective
    el = eff.nuclear.elastic
    dev = _elastic_device(eff, wp.float64, "cpu")
    host = ed.pack_elastic(el.table, eff.geometry.materials, rows=el.rows)
    assert dev.el_mass.dtype == wp.float64 and "el_mass" in dev.readback()
    np.testing.assert_array_equal(dev.el_mass.numpy(), host.arrays["target_mass"])
    # the summary hash changes when one mass value changes (both precisions)
    h2 = ed.ElasticHost(**{**host.__dict__, "arrays": {**host.arrays}})
    h2.arrays["target_mass"] = host.arrays["target_mass"].copy()
    h2.arrays["target_mass"][1] = np.nextafter(h2.arrays["target_mass"][1], np.inf)
    for prec in ("float32", "float64"):
        assert ed.host_sha256(h2, prec) != ed.host_sha256(host, prec)
    # f32 vs f64: the float64 mass bytes are identical, the precision-cast arrays differ
    c32, c64 = ed._cast(host, np.float32), ed._cast(host, np.float64)
    assert (
        c32["el_mass"].tobytes() == c64["el_mass"].tobytes() and c32["el_mass"].dtype == np.float64
    )
    assert ed.host_sha256(host, "float32") != ed.host_sha256(host, "float64")
    assert all(c32[k].tobytes() != c64[k].tobytes() for k in ("grid", "target_mass"))
    # corruption of the cached device array after caching is detected on re-verification
    ed.clear_elastic_device_cache()
    bad = _elastic_device(eff, wp.float64, "cpu")
    bad.verified = False
    arr = bad.el_mass.numpy()
    arr[1] *= 1.0000001
    bad.el_mass.assign(arr)
    with pytest.raises(ValueError, match="sha256"):
        _elastic_device(eff, wp.float64, "cpu")
    ed.clear_elastic_device_cache()


def test_no_table_array_is_uploaded_per_run(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """The elastic device (all arrays, the float64 masses included) is built once per process and
    key: a second ``run_warp_range`` constructs no ``ElasticDevice`` and hands the kernel the very
    same device arrays."""
    import ionmc.transport.warp_driver as wd

    ed.clear_elastic_device_cache()
    built: list[int] = []
    orig = ed.ElasticDevice.__init__

    def counting(self: Any, *a: Any, **k: Any) -> None:
        built.append(1)
        orig(self, *a, **k)

    monkeypatch.setattr(ed.ElasticDevice, "__init__", counting)
    seen: list[Any] = []
    orig_in = wd._nuclear_inputs

    def spy(*a: Any, **k: Any) -> Any:
        out = orig_in(*a, **k)
        seen.append(out[1])
        return out

    monkeypatch.setattr(wd, "_nuclear_inputs", spy)
    eff = Simulation(_warp_cfg(tid, 100.0, 4, "float64")).effective
    run_warp_range(eff, 0, 4, "cpu")
    run_warp_range(eff, 0, 4, "cpu")
    assert len(built) == 1 and len(seen) == 2
    dev = wd._elastic_device(eff, wp.float64, "cpu")
    for nd in seen:
        assert nd.el_mass is dev.el_mass and nd.el_grid is dev.grid and nd.el_cum is dev.cum_sigma
    ed.clear_elastic_device_cache()


# ---- C1-ext forced counters on warp-cpu -------------------------------------------------------
def _forced(tid: str, monkeypatch: pytest.MonkeyPatch) -> Any:
    _scale_elastic(monkeypatch, 30.0)
    cfg = _warp_cfg(tid, 150.0, 32, "float64", allow_invalid=True)
    return Simulation(replace(cfg, diagnostics=DiagnosticsOptions())).run()


@pytest.mark.parametrize(
    "counter", ["queue_overflow", "genealogy_overflow", "nuclear_conservation"]
)
def test_forced_runtime_limits_fail_closed_with_elastic_warp_cpu(
    tid: str, monkeypatch: pytest.MonkeyPatch, counter: str
) -> None:
    if counter == "queue_overflow":
        monkeypatch.setattr(kn, "STACK_CAPACITY", 0)
    elif counter == "genealogy_overflow":
        monkeypatch.setattr(kn, "CHILD_LIMIT", 0)
    else:
        monkeypatch.setattr(kn, "LEDGER_TOL", -1.0)
    res = _forced(tid, monkeypatch)
    assert res.counters.as_dict()[counter] > 0 and not res.valid
    assert res.energy_balance.relative_residual <= 1e-12


def test_forced_elastic_majorant_violation_warp_cpu(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    orig = ElasticTable.material_rows

    def low(self: ElasticTable, material: Any, f_e: float = 0.02) -> Any:
        m = orig(self, material, f_e)
        return replace(
            m, sigma_hat_window=m.sigma_hat_window * 0.0, sigma_hat_end=m.sigma_hat_end * 0.0
        )

    monkeypatch.setattr(ElasticTable, "material_rows", low)
    cfg = _warp_cfg(tid, 150.0, 8, "float64", allow_invalid=True)
    res = Simulation(replace(cfg, diagnostics=DiagnosticsOptions())).run()
    assert not res.valid and res.counters.majorant_violation > 0


# ---- CUDA (skipped without a device; qualified by hr5c on the host runner) -----------------------
@pytest.mark.cuda
@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_cuda_matches_warp_cpu_elastic(
    tid: str, precision: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    _scale_elastic(monkeypatch, 30.0)
    n = 32
    eff = Simulation(_warp_cfg(tid, 150.0, n, precision)).effective
    cpu = run_warp_range(eff, 0, n, "cpu")
    cu = run_warp_range(eff, 0, n, "cuda:0")
    sc, cc = _sums(cpu)
    sg, cg = _sums(cu)
    assert cc == cg or precision == "float32"
    assert max(abs(a - b) for a, b in zip(sc, sg, strict=True)) <= (
        1e-9 if precision == "float64" else 1e-2
    )
    assert cu.meta["elastic_device_sha256"] == cpu.meta["elastic_device_sha256"]
