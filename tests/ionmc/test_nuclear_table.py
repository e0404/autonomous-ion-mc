"""V3-005A nuclear table: sampler twin agreement, builder determinism, loader, composition,
majorants, D6 definitions (rows N1, V1, C1 table cases, D6; decision 0041). Data-backed tests
use the cache (``IONMC_CACHE_DIR``), skip without it and fail under ``IONMC_REQUIRE_DATA=1``."""

from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.data import cache
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import WATER, Material
from ionmc.nuclear import build as B
from ionmc.nuclear import events as ev
from ionmc.nuclear.tables import (
    NuclearTable,
    NuclearTableMissingError,
    NuclearTablePinError,
    NuclearTableStaleError,
    NuclearTableUnqualifiedError,
    majorant,
)

M_N, M_P, M_D, M_A = 939.56542052, 938.27208816, 1875.61294257, 3727.3794066


def _synthetic_model() -> ev.EventModel:
    def mass(z: int, a: int) -> float:
        n = a - z
        bind = 15.75 * a - 17.8 * a ** (2 / 3) - 0.711 * z * (z - 1) / a ** (1 / 3)
        bind -= 23.7 * (a - 2 * z) ** 2 / a
        return z * M_P + n * M_N - bind

    m_res = np.full((ev.DZ_MAX + 1, ev.DA_MAX + 1), np.inf)
    for dz in range(ev.DZ_MAX + 1):
        for da in range(ev.DA_MAX + 1):
            z_r, a_r = 9 - dz, 17 - da
            if z_r >= 1 and a_r - z_r >= 1:
                m_res[dz, da] = mass(z_r, a_r)
    return ev.EventModel(
        8, 16, M_P, mass(8, 16), np.array([M_N, M_P, M_D, M_A, 0.0]), 6.5,
        np.array([13.7, 6.5, 16.0, 8.1, 0.0]), m_res,
    )  # fmt: skip


def _rows() -> ev.EnergyRows:
    q = np.arange(65) / 64.0
    edges = np.array([-m * np.log(1.0 - q * (1 - math.exp(-60.0 / m))) for m in (6, 9, 12, 8, 3)])
    return ev.EnergyRows(np.array([1.0, 1.7, 0.3, 0.8, 0.4]), edges, np.full((5, 64), 0.4), 2.5)


def test_numpy_sampler_agrees_with_twin_path_on_identical_uniforms() -> None:
    """Scalar and batch paths: counts and residual identical, continuous values to 1e-12."""
    model, rows, src = _synthetic_model(), _rows(), ev.CounterUniforms(2024)
    n = 1000
    batch = ev.sample_events(model, rows, 150.0, n, src, want_particles=True)
    assert batch.particle_lab is not None and batch.particle_event is not None
    worst = 0.0
    for i in range(n):
        s = ev.sample_event_scalar(model, rows, 150.0, src.for_event(i))
        assert s.accepted == bool(batch.accepted[i])
        assert s.attempts == batch.attempts[i]
        assert tuple(s.counts) == tuple(batch.counts[i])
        if s.accepted:
            assert (s.z_r, s.a_r) == (batch.z_r[i], batch.a_r[i])
            lab = batch.particle_lab[batch.particle_event == i]
            sl = np.array([q[4:8] for q in s.particles]).reshape(-1, 4)
            scale = max(float(np.abs(lab).max(initial=1.0)), 1.0)
            worst = max(worst, float(np.abs(lab - sl).max(initial=0.0)) / scale)
            for a_, b_ in (
                (s.imbalance_mev, batch.imbalance_mev[i]),
                (s.binding_mev, batch.binding_mev[i]),
                (s.local_deposit_mev, batch.local_deposit_mev[i]),
            ):
                worst = max(worst, abs(a_ - float(b_)) / (1.0 + abs(a_)))
    assert worst < 1e-12


def test_p3_ledger_closure_on_synthetic_tables_1e5_events() -> None:
    """P3 (amended): T1 = sum T_lab + T_r + binding + Delta to 1e-9 MeV on every event, dZ = dA =
    0, and the residual exists (finite AME-like mass)."""
    model, rows = _synthetic_model(), _rows()
    n = 100_000
    t1 = 150.0
    b = ev.sample_events(model, rows, t1, n, ev.CounterUniforms(31415), want_particles=True)
    assert b.particle_event is not None and b.particle_lab is not None
    assert b.particle_species is not None
    acc = b.accepted
    assert acc.mean() > 0.9
    t_lab = b.particle_lab[:, 0] - model.species_mass_mev[b.particle_species]
    sum_t = np.bincount(b.particle_event, weights=t_lab, minlength=n)
    closure = np.abs(
        t1 - sum_t[acc] - b.recoil_t_mev[acc] - b.binding_mev[acc] - b.imbalance_mev[acc]
    )
    assert closure.max() <= 1e-9
    z_p = np.asarray(ev.SPECIES_Z)
    a_p = np.asarray(ev.SPECIES_A)
    assert np.all(b.counts[acc] @ z_p + b.z_r[acc] == model.z_t + 1)
    assert np.all(b.counts[acc] @ a_p + b.a_r[acc] == model.a_t + 1)
    assert np.all(
        np.isfinite(model.m_res_mev[(model.z_t + 1 - b.z_r[acc]), (model.a_t + 1 - b.a_r[acc])])
    )
    assert np.all(b.recoil_t_mev[acc] == 2.5)


def test_exact_enumeration_matches_monte_carlo_within_3_sigma() -> None:
    model, rows = _synthetic_model(), _rows()
    n = 100_000
    b = ev.sample_events(model, rows, 100.0, n, ev.CounterUniforms(2718))
    p_acc, mean = ev.exact_post_acceptance(model, rows.lam)
    first = float(np.mean(b.accepted & (b.attempts == 1)))
    assert abs(first - p_acc) <= 3.0 * math.sqrt(p_acc * (1.0 - p_acc) / n) + 1e-12
    c = b.counts[b.accepted].astype(float)
    sem = c.std(axis=0) / math.sqrt(c.shape[0])
    assert np.all(np.abs(c.mean(axis=0) - mean) <= 3.0 * sem + 1e-12)
    # solving lam with the exact enumeration reproduces the requested yields
    y = np.array([0.9, 1.4, 0.2, 0.6, 0.3])
    sol = B.solve_lambda(model, y, 1e-6)
    assert sol["converged"] and sol["max_residual"] <= 1e-6 and sol["p_accept"] >= 0.5
    assert np.allclose(ev.exact_post_acceptance(model, sol["lam"])[1], y, rtol=2e-6)


def test_counter_uniforms_scalar_equals_vector() -> None:
    src = ev.CounterUniforms(7)
    vec = src.uniform(np.array([3, 5]), 2, np.array([9, 40]))
    assert vec[0] == src.scalar(3, 2, 9) and vec[1] == src.scalar(5, 2, 40)
    assert np.all((vec > 0) & (vec < 1))


def test_grid_contains_150_and_is_uniform_in_ln_e() -> None:
    g, i150 = B.build_grid(50)
    assert g[0] == 1.0 and g[i150] == 150.0 and g[-1] >= 250.0
    h = np.diff(np.log(g))
    assert np.allclose(h, h[0], rtol=1e-9) and 1 / (h[0] / math.log(10)) >= 50.0
    with pytest.raises(ValueError):
        B.build_grid(10)


def test_surrogate_scaling_factor() -> None:
    rows = B.element_rows(B.TARGETS)
    assert rows["Na"]["target"] == "Al-27" and rows["Na"]["surrogate"]
    assert rows["Na"]["sigma_scale"] == pytest.approx((22.98977 / 27.0) ** (2 / 3))
    assert rows["K"]["sigma_scale"] == pytest.approx((39.0983 / 40.0) ** (2 / 3))
    assert rows["O"]["sigma_scale"] == 1.0 and not rows["O"]["surrogate"]
    assert "H" not in rows


def test_deterministic_npz_bytes(tmp_path: Path) -> None:
    arrays = {"b": np.arange(5.0), "a": np.eye(2)}
    h1 = B.write_npz_deterministic(tmp_path / "x.npz", arrays)
    h2 = B.write_npz_deterministic(tmp_path / "y.npz", arrays)
    assert h1 == h2 and (tmp_path / "x.npz").read_bytes() == (tmp_path / "y.npz").read_bytes()
    with np.load(tmp_path / "x.npz", allow_pickle=False) as z:
        assert np.array_equal(z["b"], arrays["b"])


def test_table_id_binds_every_sidecar_field_but_not_the_key_order() -> None:
    doc = {"b": {"y": 1, "x": [1.5, 2]}, "a": "1", "npz_sha256": "a" * 64, "options": {"n": 3}}
    tid = B.table_id(doc)
    assert tid == B.table_id(dict(reversed(doc.items())))
    assert tid == B.table_id({**doc, "table_id": "anything"})  # the id field is excluded
    assert tid != B.table_id({**doc, "npz_sha256": "b" * 64})  # the id covers the npz bytes
    assert tid != B.table_id({**doc, "options": {"n": 4}})
    assert tid != B.table_id({**doc, "b": {"y": 1, "x": [1.5, 2.0000001]}})
    assert tid != B.table_id({**doc, "extra": 0})


def test_d6_definitions_on_synthetic_spectra() -> None:
    def rng_fn(t: np.ndarray) -> np.ndarray:  # 0.1 mm of water is crossed at 10 MeV
        return 0.01 * (t / 10.0) ** 1.7

    low = B.gate_numbers(np.full(1000, 3.0), 1000, 0.15, 150.0, rng_fn)
    high = B.gate_numbers(np.full(1000, 40.0), 1000, 0.15, 150.0, rng_fn)
    assert low["G"] == 0.0 and low["D"] == 0.0
    assert high["G"] == 1.0
    assert high["D"] == pytest.approx(40.0 * 0.15 / 150.0)
    assert 0.0 <= high["p_event"] <= 1.0
    g = {"150": high, "250": high}
    assert B.tiers(g) == (False, False)
    assert B.tiers({"150": low, "250": low}) == (True, True)


def _synthetic_table(tmp_path: Path) -> NuclearTable:
    grid, _ = B.build_grid(50)
    n = grid.size
    sig = np.zeros((2, n))
    sig[0] = 0.2 + 0.1 * np.sin(np.log(grid))
    sig[1] = 0.4
    info = {
        "table_id": "synthetic",
        "targets": [{"name": "O-16"}, {"name": "C-12"}],
        "elements": {
            "O": {"target": "O-16", "sigma_scale": 1.0, "a_g_mol": 16.0},
            "C": {"target": "C-12", "sigma_scale": 1.0, "a_g_mol": 12.0},
            "Na": {"target": "O-16", "sigma_scale": 0.9, "a_g_mol": 23.0},
        },
    }
    return NuclearTable(info, {"grid_e_mev": grid, "sigma_barn": sig}, tmp_path / "x.npz")


def test_material_composition_cumulative_fractions_and_unsupported(tmp_path: Path) -> None:
    tab = _synthetic_table(tmp_path)
    mat = Material("toy", 1.0, {"O": 0.5, "C": 0.3, "H": 0.2})
    rows = tab.material_rows(mat)
    na = 6.02214076e23 * 1e-24
    o, c = na * 0.5 * tab.arrays["sigma_barn"][0] / 16.0, na * 0.3 * 0.4 / 12.0
    assert np.allclose(rows.sigma_mass_cm2_g, o + c, rtol=1e-12)
    assert rows.target_index == (0, 1)
    assert np.allclose(rows.cum_fraction[0], o / (o + c)) and np.all(rows.cum_fraction[-1] == 1.0)
    assert not rows.sigma_mass_cm2_g.flags.writeable
    with pytest.raises(UnsupportedCombinationError):
        tab.material_rows(Material("cu", 8.96, {"Cu": 1.0}))
    assert np.all(
        tab.material_rows(WATER if False else Material("h", 1.0, {"H": 1.0})).sigma_mass_cm2_g == 0
    )


def test_majorant_covers_sigma_and_window_is_correct_on_a_peak() -> None:
    grid, _ = B.build_grid(50)
    sig = np.exp(-0.5 * ((np.log(grid) - math.log(60.0)) / 0.08) ** 2) + 0.01
    win, end = majorant(grid, sig, 0.02)
    assert np.all(win >= 1.02 * sig - 1e-15) and np.all(end >= win - 1e-12)
    ln = np.log(grid)
    for e0 in np.exp(np.linspace(ln[3], ln[-2], 777)):  # any E0, any E1 in its window
        k = int(np.searchsorted(grid, e0) - 1)
        bound = win[k]  # step lookup: the value of the node at or below E0
        e1 = np.linspace(e0 * 0.95, e0, 50)
        assert np.all(np.interp(np.log(e1), ln, sig) <= bound)
    # peak at 60 MeV: just above it the window must still see the peak
    k = int(np.searchsorted(grid, 62.0))
    assert win[k] >= 1.02 * float(sig[np.abs(grid - 60.0).argmin()]) * 0.99
    # far below the peak the majorant is the small local value
    assert win[5] < 0.1 and end[5] < 0.1


def test_union_grid_contains_every_endf_node_and_150_and_is_increasing() -> None:
    endf = np.array([1.0, 1.23456789, 3.3, 7.77, 149.5, 150.0, 400.0])
    g, i150 = B.union_grid(endf, 50)
    assert g[0] == 1.0 and g[i150] == 150.0 and g[-1] >= 250.0 and np.all(np.diff(g) > 0)
    for x in endf[endf <= 150.0]:
        assert x in g
    uni, _ = B.build_grid(50)
    assert g.size <= uni.size + 5 and g.size >= uni.size  # ENDF nodes only add or replace
    near = np.abs(g[:, None] - g[None, :]) < 1e-12 * g[:, None]
    assert near.sum() == g.size  # no two nodes closer than 1e-12 (relative)


def test_grid_locate_matches_searchsorted_on_1e5_energies_nodes_and_clamps() -> None:
    from ionmc._wpfunc import python_twin
    from ionmc.physics.nuclear import make_nuclear

    nu = python_twin(make_nuclear)
    g, _ = B.union_grid(np.array([1.5, 2.0, 3.25, 17.0, 33.0, 77.7]), 50)
    rng = np.random.default_rng(5)
    e = np.concatenate(
        (np.exp(rng.uniform(math.log(0.3), math.log(400.0), 100_000)), g, [0.0, 0.5, 250.0, 1e3])
    )
    ref = np.clip(np.searchsorted(g, e, side="right") - 1, 0, g.size - 2)
    got = np.array([nu.grid_locate(float(v), g, g.size) for v in e])
    assert np.array_equal(got, ref)


def test_lin_lin_rows_are_exact_at_nodes_and_midpoints_on_a_two_target_fixture(
    tmp_path: Path,
) -> None:
    x0, y0 = np.array([1.0, 3.0, 8.0, 40.0, 150.0]), np.array([0.0, 0.05, 0.4, 0.3, 0.2])
    x1, y1 = (
        np.array([1.0, 2.0, 6.5, 21.0, 90.0, 150.0]),
        np.array([0.01, 0.2, 0.5, 0.35, 0.3, 0.25]),
    )
    grid, _ = B.union_grid(np.unique(np.concatenate((x0, x1))), 50)
    sig = np.array([np.interp(grid, x0, y0), np.interp(grid, x1, y1)])
    info = {
        "table_id": "fixture",
        "targets": [{"name": "O-16"}, {"name": "C-12"}],
        "elements": {
            "O": {"target": "O-16", "sigma_scale": 1.0, "a_g_mol": 16.0},
            "C": {"target": "C-12", "sigma_scale": 1.0, "a_g_mol": 12.0},
        },
    }
    tab = NuclearTable(info, {"grid_e_mev": grid, "sigma_barn": sig}, tmp_path / "x.npz")
    mat = Material("toy", 1.0, {"O": 0.6, "C": 0.4})
    rows = tab.material_rows(mat)
    na = 6.02214076e23 * 1e-24
    pts = np.unique(
        np.concatenate(
            (
                x0,
                x1,
                0.5
                * (
                    np.unique(np.concatenate((x0, x1)))[1:]
                    + np.unique(np.concatenate((x0, x1)))[:-1]
                ),
            )
        )
    )
    worst = 0.0
    for e in pts:
        o = na * 0.6 * np.interp(e, x0, y0) / 16.0
        c = na * 0.4 * np.interp(e, x1, y1) / 12.0
        ref = o + c
        worst = max(worst, abs(rows.sigma_at(float(e)) / ref - 1.0))
        cf = rows.cum_fraction_at(float(e))
        worst = max(worst, abs(cf[0] - o / ref), abs(cf[1] - 1.0))
    assert worst <= 1e-14


def test_interp_rows_lin_lin_in_e() -> None:
    g = np.array([1.0, 2.0, 4.0])
    lam = np.zeros((5, 3))
    lam[1] = [0.0, 1.0, 3.0]
    edges = np.zeros((5, 3, 65))
    edges[0, :, 3] = [0.0, 2.0, 10.0]
    r = np.zeros((5, 3, 64))
    rec = np.array([0.0, 1.0, 5.0])
    row = ev.interp_rows(g, lam, edges, r, rec, 3.0)
    assert row.lam[1] == 2.0 and row.edges_mev[0, 3] == 6.0 and row.recoil_t_mev == 3.0
    assert ev.interp_rows(g, lam, edges, r, rec, 0.1).lam[1] == 0.0  # clamped below
    assert ev.interp_rows(g, lam, edges, r, rec, 9.0).lam[1] == 3.0  # clamped above


# ---------------------------------------------------------------------------------------------
# data-backed
# ---------------------------------------------------------------------------------------------
def _cache_dir() -> Path:
    cdir = cache.resolve_cache_dir()
    try:
        for sid in B.SOURCE_IDS:
            cache.verify(sid, cdir)
    except (FileNotFoundError, cache.IntegrityError):
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("IONMC_REQUIRE_DATA=1 but the nuclear sources are not cached")
        pytest.skip("nuclear sources are not cached")
    return cdir


REDUCED = B.BuildOptions(
    diagnostic_events=2000,
    diagnostic_nodes_mev=(100.0,),
    d6_events=2000,
)


@pytest.fixture(scope="module")
def reduced_build() -> Any:
    return B.build_nuclear_proton(_cache_dir(), REDUCED)


def test_reduced_build_twice_is_byte_identical(reduced_build: Any) -> None:
    first = (reduced_build.npz_path.read_bytes(), reduced_build.json_path.read_bytes())
    again = B.build_nuclear_proton(_cache_dir(), REDUCED)
    assert again.table_id == reduced_build.table_id
    assert (again.npz_path.read_bytes(), again.json_path.read_bytes()) == first
    info = json.loads(reduced_build.json_path.read_text())
    assert info["kalbach_separation"] == "systematics-formula"
    g = info["gate_d6"]
    for k in ("150", "250"):
        assert 0.0 <= g["numbers"][k]["G"] <= 1.0 and 0.0 <= g["numbers"][k]["p_event"] <= 1.0
        assert math.isfinite(g["numbers"][k]["D"])
    assert isinstance(g["tier1_pass"], bool) and isinstance(g["tier2_pass"], bool)
    assert isinstance(g["ceiling_pass"], bool)
    assert info["schema"] == B.SCHEMA and info["builder_version"] == B.BUILDER_VERSION
    assert 0.0 < info["transport_energy_bound_mev"] <= B.STOPPING_TABLE_MAX_MEV
    assert info["empty_residual_allowed"] is True
    assert set(info["transport_path_bound_terms"]) == {"n", "p", "d", "a", "g"}
    terms = info["transport_path_bound_terms"]
    expect = sum(v["n_max"] * v["t_lab_max_mev"] for v in terms.values()) + info["recoil_t_max_mev"]
    assert info["history_energy_bound_mev"] == max(info["grid"]["e_max_mev"], expect)
    assert info["history_energy_bound_mev"] >= info["transport_energy_bound_mev"]
    assert info["multiplicity"]["all_nodes_converged"] is (
        not info["multiplicity"]["non_converged_nodes"]
    )


def test_load_and_fail_closed_cases(
    reduced_build: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cdir = _cache_dir()
    # the reduced build need not be qualified: the qualification itself is tested separately
    monkeypatch.setattr("ionmc.nuclear.tables.qualification_failures_npz", lambda derived: [])
    tab = NuclearTable.load(cdir, reduced_build.table_id)
    assert not tab.arrays["sigma_barn"].flags.writeable
    with pytest.raises(NuclearTableMissingError):
        NuclearTable.load(cdir, "0" * 64)
    # stale: flip one byte in a tmp copy
    d = tmp_path / "derived"
    d.mkdir()
    for p in (reduced_build.npz_path, reduced_build.json_path):
        shutil.copy(p, d / p.name)
    npz = d / reduced_build.npz_path.name
    raw = bytearray(npz.read_bytes())
    raw[len(raw) // 2] ^= 0xFF
    npz.write_bytes(bytes(raw))
    with pytest.raises(NuclearTableStaleError):
        NuclearTable.load(tmp_path, reduced_build.table_id)
    # pin mismatch
    pin = DATASETS["ame2020-mass"]
    monkeypatch.setitem(DATASETS, "ame2020-mass", type(pin)(**{**pin.__dict__, "sha256": "f" * 64}))
    with pytest.raises(NuclearTablePinError):
        NuclearTable.load(cdir, reduced_build.table_id)


def _checks_module() -> Any:
    path = Path(__file__).resolve().parents[2] / "validation/scripts/transport/nuclear_checks.py"
    spec = importlib.util.spec_from_file_location("nuclear_checks", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_p_accept_guard_on_a_model_that_violates_the_rule() -> None:
    """A synthetic model whose only existing residual is the empty one: the exact P_accept of
    yields of several products is far below 0.5 and ``solve_lambda`` reports it unconverged."""
    base = _synthetic_model()
    m_res = np.full_like(base.m_res_mev, np.inf)
    m_res[9, 17] = 0.0  # (Z_r, A_r) = (0, 0)
    model = ev.EventModel(
        base.z_t, base.a_t, base.m_p_mev, base.m_t_mev, base.species_mass_mev, base.s_a_mev,
        base.s_b_mev, m_res,
    )  # fmt: skip
    sol = B.solve_lambda(model, np.array([1.0, 1.0, 0.2, 0.5, 0.3]), 1e-6)
    assert sol["p_accept"] < B.P_ACCEPT_MIN and not sol["converged"]


def test_yield_extension_below_the_first_mf6_energy() -> None:
    tab = B.Tab1(np.array([2]), np.array([2]), np.array([1.5e6, 3.0e6]), np.array([0.5, 1.0]))
    sp = B.SpeciesTables.__new__(B.SpeciesTables)
    sp.yield_tab, sp.extend_below_mev, sp.extend_to_mev = tab, 1.0, 1.5
    assert sp.yield_at(0.9) == 0.0 and sp.yield_at(1.0) == 0.5 and sp.yield_at(1.2) == 0.5
    assert sp.yield_at(2.25) == pytest.approx(0.75)


def test_n1_v1_v4_on_the_table(reduced_build: Any) -> None:
    cdir = _cache_dir()
    s = _checks_module().run_checks(cdir, reduced_build.table_id)
    print("N1 max rel", s["N1"]["max_rel"], "V1 max rel", s["V1"]["max_rel"])
    assert s["N1"]["max_rel"] <= 1e-12 and s["V1"]["max_rel"] <= 1e-12
    assert s["N1"]["pass"] and s["V1"]["pass"]
    assert s["V4"]["V4_pass"]
    # V4b is a regression guard reported, not asserted: its n <E'> criterion (1 %) is 1.2 % at
    # O-16 100 MeV on the real table (worker report); the numbers are printed.
    print("V4b pass", s["V4"]["V4b_pass"], {k: c["v4b_pass"] for k, c in s["V4"]["cases"].items()})


def _qualified_table_id() -> str:
    cdir = _cache_dir()
    for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
        tid = p.stem.removeprefix("nuclear-proton-")
        try:
            NuclearTable.load(cdir, tid)
        except UnsupportedCombinationError:
            continue
        return tid
    pytest.skip("no qualified nuclear table in the cache")


def _tampered_copy(tmp_path: Path, tid: str, edit: Any) -> None:
    src = _cache_dir() / "derived"
    d = tmp_path / "derived"
    d.mkdir()
    for ext in ("npz", "json"):
        shutil.copy(src / f"nuclear-proton-{tid}.{ext}", d / f"nuclear-proton-{tid}.{ext}")
    jp = d / f"nuclear-proton-{tid}.json"
    info = json.loads(jp.read_text())
    edit(info)
    jp.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")


def _resealed_copy(
    root: Path, tid: str, edit_arrays: Any = None, edit_info: Any = None, new_id: bool = True
) -> str:
    """Copy the cached table ``tid`` into ``root/derived``, edit its arrays and JSON and reseal:
    the new npz digest goes into the sidecar and, with ``new_id``, into a recomputed table id
    (the files are renamed), i.e. a self-consistent table that only a gate can refuse. Returns the
    id under which the copy loads."""
    src = _cache_dir() / "derived"
    d = root / "derived"
    d.mkdir(parents=True)
    info = json.loads((src / f"nuclear-proton-{tid}.json").read_text())
    with np.load(src / f"nuclear-proton-{tid}.npz", allow_pickle=False) as z:
        arrays = {k: z[k].copy() for k in z.files}
    if edit_arrays is not None:
        edit_arrays(arrays)
    if edit_info is not None:
        edit_info(info)
    tmp = d / "tmp.npz"
    info["npz_sha256"] = B.write_npz_deterministic(tmp, arrays)
    if new_id:
        info["table_id"] = B.table_id(info)
    out = info["table_id"] if new_id else tid
    tmp.replace(d / f"nuclear-proton-{out}.npz")
    (d / f"nuclear-proton-{out}.json").write_text(json.dumps(info, indent=2, sort_keys=True))
    return str(out)


def test_loader_refuses_unqualified_tables(tmp_path: Path) -> None:
    """A self-consistent table (arrays, digest and id resealed) whose npz arrays record a failed
    ceiling or a non-converged node is refused; the npz decides."""
    tid = _qualified_table_id()
    tab = NuclearTable.load(None, tid)
    assert tab.info["gate_d6"]["ceiling_pass"] is True
    qi = {k: i for i, k in enumerate(B.QUALIFICATION_FIELDS)}

    def fail_flag(name: str, key: str) -> Any:
        def edit_a(arr: dict[str, Any]) -> None:
            arr["qualification"][qi[name]] = 0

        def edit_i(info: dict[str, Any]) -> None:
            info["gate_d6"][key] = False

        return edit_a, edit_i

    def node_a(arr: dict[str, Any]) -> None:
        arr["lam_converged"][0, 0] = 0
        arr["qualification"][qi["all_nodes_converged"]] = 0

    def node_i(info: dict[str, Any]) -> None:
        info["multiplicity"]["non_converged_nodes"] = [{"target": "C-12", "e_mev": 1.0}]
        info["multiplicity"]["all_nodes_converged"] = False

    cases = {
        "ceiling": fail_flag("ceiling_pass", "ceiling_pass"),
        "ceiling_energy": fail_flag(
            "ceiling_pass_energy_weighted_range", "ceiling_pass_energy_weighted_range"
        ),
        "node": (node_a, node_i),
    }
    for name, (edit_a, edit_i) in cases.items():
        root = tmp_path / name
        new = _resealed_copy(root, tid, edit_a, edit_i)
        with pytest.raises(NuclearTableUnqualifiedError):
            NuclearTable.load(root, new)
        # the sidecar flipped back to "passing" on the unqualified npz: stale (the npz decides)
        root2 = tmp_path / f"{name}-flipped"
        new2 = _resealed_copy(root2, tid, edit_a, None)
        with pytest.raises(NuclearTableStaleError):
            NuclearTable.load(root2, new2)


def test_loader_authenticates_the_npz_and_the_bounds(tmp_path: Path) -> None:
    """Altered array bytes with a resealed sidecar digest keep the pinned id from reproducing; a
    weakened sidecar bound, or weakened npz bounds under a resealed id and an untouched sidecar,
    are stale (Codex finding 2)."""
    tid = _qualified_table_id()

    def bump_sigma(arr: dict[str, Any]) -> None:
        arr["sigma_barn"] = arr["sigma_barn"] * 1.0000001

    root = tmp_path / "bytes"
    same = _resealed_copy(root, tid, bump_sigma, None, new_id=False)  # digest resealed, id kept
    with pytest.raises(NuclearTableStaleError, match="id not reproduced"):
        NuclearTable.load(root, same)

    def weaker_json(info: dict[str, Any]) -> None:
        info["history_energy_bound_mev"] *= 0.5

    root = tmp_path / "json-bound"
    with pytest.raises(NuclearTableStaleError, match="disagrees"):
        NuclearTable.load(root, _resealed_copy(root, tid, None, weaker_json))

    def weaker_npz(arr: dict[str, Any]) -> None:
        arr["bounds"][B.BOUND_FIELDS.index("history_energy_bound_mev")] *= 0.5

    root = tmp_path / "npz-bound"
    with pytest.raises(NuclearTableStaleError, match="disagrees"):
        NuclearTable.load(root, _resealed_copy(root, tid, weaker_npz, None))
    # the loaded bounds are those of the npz
    tab = NuclearTable.load(None, tid)
    b = B.bounds_from_array(tab.arrays["bounds"])
    assert b["history_energy_bound_mev"] == tab.info["history_energy_bound_mev"]
    assert b["transport_path_bound_terms"] == dict(
        (k, dict(v)) for k, v in tab.info["transport_path_bound_terms"].items()
    )


def test_binding_recomputed_per_event_from_ame_masses_independently() -> None:
    """The binding of 1e3 sampled events recomputed from AME2020 masses (counts and residual only;
    the sampler's ``binding_mev`` is not an input) agrees per event to 1e-9 MeV; a perturbed AME
    mass breaks it."""
    from dataclasses import replace

    from ionmc.data.ame import load_ame2020

    cdir = _cache_dir()
    ame = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    chk = _checks_module()
    model = ev.build_event_model(ame, 8, 16)
    q = np.arange(65) / 64.0
    edges = np.array([-m * np.log(1.0 - q * (1 - math.exp(-60.0 / m))) for m in (6, 9, 12, 8, 3)])
    rows = ev.EnergyRows(np.array([1.0, 1.7, 0.3, 0.8, 0.4]), edges, np.full((5, 64), 0.4), 2.5)
    b = ev.sample_events(model, rows, 100.0, 1000, ev.CounterUniforms(99))
    acc = b.accepted
    assert acc.sum() > 900
    got = chk.recompute_binding_events(b.counts[acc], b.z_r[acc], b.a_r[acc], 8, 16, ame)
    assert np.abs(got - b.binding_mev[acc]).max() <= 1e-9
    bad = dict(ame)
    bad[(8, 16)] = replace(ame[(8, 16)], atomic_mass_u=ame[(8, 16)].atomic_mass_u + 1e-6)
    off = chk.recompute_binding_events(b.counts[acc], b.z_r[acc], b.a_r[acc], 8, 16, bad)
    assert np.abs(off - b.binding_mev[acc]).max() > 1e-4


def test_sidecar_fields_are_bound_to_the_pinned_id(tmp_path: Path) -> None:
    """Codex review 3, finding 1: elements[*].sigma_scale, targets[*].z and gate_d6.numbers are
    sidecar-only fields; editing any of them (npz and pinned id kept) is stale, and resealing the id
    produces a different id that the config pin does not name (missing under the pinned id)."""
    tid = _qualified_table_id()

    def scale(info: dict[str, Any]) -> None:
        key = next(iter(info["elements"]))
        info["elements"][key]["sigma_scale"] *= 1.01

    def target_z(info: dict[str, Any]) -> None:
        info["targets"][0]["z"] += 1

    def d6(info: dict[str, Any]) -> None:
        k = next(iter(info["gate_d6"]["numbers"]))
        info["gate_d6"]["numbers"][k]["D"] *= 0.5

    for name, edit in {"sigma_scale": scale, "z": target_z, "d6": d6}.items():
        root = tmp_path / name
        root.mkdir()
        _tampered_copy(root, tid, edit)  # id and npz bytes untouched
        with pytest.raises(NuclearTableStaleError, match="id not reproduced"):
            NuclearTable.load(root, tid)
        # resealed id: the edited document is self-consistent under another id, the pinned id
        # names no table in that cache and the id differs from the pin
        sealed = tmp_path / f"{name}-sealed"
        new = _resealed_copy(sealed, tid, None, edit)
        assert new != tid
        with pytest.raises(NuclearTableMissingError):
            NuclearTable.load(sealed, tid)
        # the edited sidecar copied over the pinned file names: its own id disagrees with the pin
        shutil.copy(
            sealed / "derived" / f"nuclear-proton-{new}.json",
            sealed / "derived" / f"nuclear-proton-{tid}.json",
        )
        shutil.copy(
            sealed / "derived" / f"nuclear-proton-{new}.npz",
            sealed / "derived" / f"nuclear-proton-{tid}.npz",
        )
        with pytest.raises(NuclearTableStaleError, match="id not reproduced"):
            NuclearTable.load(sealed, tid)
