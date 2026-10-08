"""V3-005A work packages C8 and C9: deuteron transport tables and the species-aware water row
(row V9, R1), ``producible(nuclear=True)``, the per-species scoring hook, ``validate()`` with
``nuclear=True`` (every configuration case of row C1), the capacity bounds of decision 0041
section 5 (amended 2026-10-07), the effective-config record and the capability report.

The data-backed tests use a built nuclear table of the cache (``IONMC_CACHE_DIR``) and skip
without one (``IONMC_REQUIRE_DATA=1`` makes that a failure)."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import ionmc
from ionmc.config import SimulationConfig, validate
from ionmc.data import cache
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import ADIPOSE_TISSUE_ICRP, AIR, LEAD, WATER
from ionmc.nuclear.tables import (
    NuclearTable,
    NuclearTableError,
    NuclearTableMissingError,
    NuclearTablePinError,
    NuclearTableStaleError,
)
from ionmc.physics.projectiles import ALPHA, DEUTERON, HELIUM3, PROTON
from ionmc.physics.stopping import BetheStoppingSource, StoppingTable
from ionmc.scoring import TallyRequest
from ionmc.simulation import Simulation
from ionmc.species import producible
from ionmc.transport.channels import MAX_PARTICLES, STRAGGLING_MARGIN, mixed_path_bound_mm
from ionmc.transport.scoring_ref import ReferenceChannelScorer
from ionmc.transport.tables import _ARRAY_FIELDS, TransportTables

MakeConfig = Callable[..., SimulationConfig]
MATS = (WATER, ADIPOSE_TISSUE_ICRP, AIR)

# ---- R1: proton tables bit-identical (hashes recorded on 9913ddc6, before C8) -----------------
PROTON_ARRAY_SHA16 = {
    "e_min_mev": "cc143326a2646c60",
    "e_max_mev": "2afc22467f0bcfd9",
    "ln_e0": "9d908ecfb6b256de",
    "inv_dln_e": "756fb1275ae6f813",
    "ln_s_mass": "ebe4c95c2a5261f5",
    "ln_r_mass": "f696c4373325ddab",
    "r_mass": "45e9aa8eda4f9246",
    "f_mass": "db5b7d524a663e90",
    "d_f": "b42c90e3930e2ec3",
    "r_min_g_cm2": "46ec7968f8b1a50d",
    "r_max_g_cm2": "c3c192a77c35b6de",
    "ln_r0": "f5667b8bad6ddbff",
    "inv_dln_r": "5a23211f69c0ec9b",
    "ln_e_of_r": "d441ade2b9202723",
    "z_over_a": "2453ee0ed2cc80bc",
    "inv_rho_xs_cm2_g": "473b7fd190782455",
    "nominal_density_g_cm3": "e72d5af20206197b",
    "water_ln_s_mass": "a3e6ec3771f3f74c",
}
PROTON_TABLE_SHA256 = "dfbacd1b680eacdc52eac7a48c761bf5433d2c6298960e75cbb7e4a0c70d154a"
_LIBM_CANARY = "51c37367c89be6f4"


def _libm_matches_recording() -> bool:
    x = np.exp(np.linspace(0.0, 5.0, 1000)) + np.log1p(np.linspace(0.0, 3.0, 1000))
    return hashlib.sha256(x.tobytes()).hexdigest()[:16] == _LIBM_CANARY


class _NistStarStandIn:
    """Carries the name of the NIST PSTAR/ASTAR source (which has no deuteron table)."""

    name = "nist-star"

    def table(self, material: Any, projectile: Any) -> StoppingTable:
        raise ValueError("no deuteron table")


@pytest.fixture(scope="module")
def bethe_src() -> BetheStoppingSource:
    return BetheStoppingSource()


@pytest.fixture(scope="module")
def p_tables(bethe_src: BetheStoppingSource) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [bethe_src.table(m, PROTON) for m in MATS], water=bethe_src.table(WATER, PROTON)
    )


@pytest.fixture(scope="module")
def d_tables(bethe_src: BetheStoppingSource) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [bethe_src.table(m, DEUTERON) for m in MATS],
        water=bethe_src.table(WATER, DEUTERON),
        projectile=DEUTERON,
    )


def test_proton_tables_are_bit_identical_to_the_pre_v3_005a_build(
    p_tables: TransportTables,
) -> None:
    if not _libm_matches_recording():
        pytest.skip("this libm differs from the machine that recorded the hashes")
    for name in (*_ARRAY_FIELDS, "water_ln_s_mass"):
        h = hashlib.sha256(np.ascontiguousarray(getattr(p_tables, name)).tobytes()).hexdigest()
        assert h[:16] == PROTON_ARRAY_SHA16[name], name
    assert p_tables.sha256 == PROTON_TABLE_SHA256
    assert all("energy_axis" not in d for d in p_tables.identity)  # proton identity unchanged
    assert p_tables.deuteron is None


def test_transport_tables_reject_other_projectiles(bethe_src: BetheStoppingSource) -> None:
    for proj in (HELIUM3, ALPHA):  # z >= 2: effective charge (decision 0038), not implemented here
        with pytest.raises(ValueError, match="z = 1"):
            TransportTables.from_stopping_tables([bethe_src.table(WATER, proj)], projectile=proj)
    with pytest.raises(ValueError, match="energy_axis"):
        TransportTables.from_stopping_tables(
            [bethe_src.table(WATER, DEUTERON)], projectile=DEUTERON, energy_axis="per_nucleon"
        )
    with pytest.raises(ValueError, match="describe"):  # mismatched projectile
        TransportTables.from_stopping_tables([bethe_src.table(WATER, PROTON)], projectile=DEUTERON)


def test_v9_deuteron_tables_obey_the_velocity_scaling_identity(
    p_tables: TransportTables, d_tables: TransportTables
) -> None:
    """S_d(E) = S_p(E m_p/m_d) and R_d(E) = (m_d/m_p) R_p(E m_p/m_d) at every deuteron grid node
    with E >= 2 MeV/u (water and tissue): relative <= 1e-3 (observed max reported by the run)."""
    ratio = DEUTERON.mass_mev / PROTON.mass_mev
    worst_s = worst_r = 0.0
    for m in (0, 1):  # water, adipose tissue
        e = np.exp(d_tables.ln_e0[m] + np.arange(d_tables.n_e) / d_tables.inv_dln_e[m])
        for en in e[e >= 2.0 * DEUTERON.a]:
            ep = float(en) / ratio
            s_d = d_tables.stopping_mass(m, float(en))
            s_p = p_tables.stopping_mass(m, ep)
            r_d = d_tables.range_g_cm2(m, float(en))
            r_p = p_tables.range_g_cm2(m, ep)
            worst_s = max(worst_s, abs(s_d / s_p - 1.0))
            worst_r = max(worst_r, abs(r_d / (ratio * r_p) - 1.0))
    print(f"V9 max relative deviation: stopping {worst_s:.3e}, range {worst_r:.3e}")
    assert worst_s <= 1.0e-3 and worst_r <= 1.0e-3
    assert d_tables.projectile == DEUTERON
    assert all(d["energy_axis"] == "total_kinetic_mev" for d in d_tables.identity)
    assert float(d_tables.e_min_mev.max()) == 2.0 and float(d_tables.e_max_mev.min()) == 1000.0
    # range construction unchanged: the exact closed form integrates dE / S on the total axis
    e1 = 20.0
    assert d_tables.energy_from_range(0, d_tables.range_g_cm2(0, e1)) == pytest.approx(e1, rel=1e-6)


def test_s_water_species_rows(p_tables: TransportTables, d_tables: TransportTables) -> None:
    both = replace(p_tables, deuteron=d_tables)
    assert both.s_water(30.0, "proton") == p_tables.s_water(30.0)
    s_d = both.s_water(60.0, "deuteron")
    assert s_d == d_tables.s_water(60.0, "deuteron")
    # same velocity: S_d(2 e) = S_p(e) to the Bethe mass-ratio accuracy
    assert s_d == pytest.approx(p_tables.s_water(30.0), rel=1e-3)
    for sp in ("triton", "helium3", "alpha", "nuclear_local"):
        with pytest.raises(ValueError):
            both.s_water(30.0, sp)
    with pytest.raises(ValueError):  # proton-only tables have no deuteron row
        p_tables.s_water(30.0, "deuteron")


def test_producible_sets() -> None:
    assert producible(PROTON) == frozenset({("proton", "primary")}) == producible(PROTON, False)
    assert producible(PROTON, nuclear=True) == frozenset(
        {
            ("proton", "primary"),
            ("proton", "secondary"),
            ("deuteron", "secondary"),
            ("nuclear_local", "secondary"),
        }
    )
    with pytest.raises(ValueError):
        producible(DEUTERON, nuclear=True)


def test_scorer_selects_the_water_row_and_mass_number_per_species(
    p_tables: TransportTables, d_tables: TransportTables
) -> None:
    from ionmc.transport.channels import ChannelPlan

    plan = ChannelPlan(
        channels=(), quantities=(), species_match=np.zeros((0, 65), dtype=np.int8),
        ch_begin=(0,), ch_end=(0,), total_size=0, n_residual=0, lookups=(), bounds={},
        memory_bytes=0, histories_per_batch=1, path_bound_mm=1.0,
    )  # fmt: skip
    tabs = replace(p_tables, deuteron=d_tables)
    sc = ReferenceChannelScorer(plan, [], 1, tables=tabs)
    sp, gp = sc.water_state(40.0, 0)
    sd, gd = sc.water_state(80.0, 1)
    assert sp == pytest.approx(p_tables.s_water(40.0)) and sd == pytest.approx(
        d_tables.s_water(80.0, "deuteron")
    )
    assert sd == pytest.approx(sp, rel=1e-3) and gd == pytest.approx(gp, rel=2e-3)
    sc.begin_step(80.0, 0.5, 0.1, 1)
    assert sc.a_step == 2 and sc.s_mid == pytest.approx(sd)
    sc.begin_step(40.0, 0.5, 0.1)  # default species 0
    assert sc.a_step == 1
    for bad in (2, 3, 4):
        with pytest.raises(ValueError):
            sc.water_state(40.0, bad)
    # without companion tables species 1 has no row (and nothing else changes)
    with pytest.raises(ValueError):
        ReferenceChannelScorer(plan, [], 1, tables=p_tables).water_state(40.0, 1)


# ---- C9: validate() with nuclear=True ---------------------------------------------------------
def _table_id() -> str:
    """A built nuclear table of the cache that the current loader accepts."""
    try:
        cdir = cache.resolve_cache_dir(None)
    except Exception:  # pragma: no cover - no cache configured
        cdir = None
    found: list[tuple[bool, str]] = []
    if cdir is not None:
        for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
            tid = p.stem.removeprefix("nuclear-proton-")
            try:
                info = NuclearTable.load(None, tid).info
            except NuclearTableError:
                continue
            if "transport_path_bound_terms" in info:  # C7b-complete tables only
                found.append((info["options"]["diagnostic_nodes_mev"] is None, tid))
    if found:
        return max(found)[1]  # a full (not reduced-diagnostics) build first
    if os.environ.get("IONMC_REQUIRE_DATA") == "1":
        pytest.fail("no loadable nuclear table in the cache")
    pytest.skip("no built nuclear table in the cache (IONMC_CACHE_DIR)")


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


def _nuc(cfg: SimulationConfig, tid: str | None, **phys: Any) -> SimulationConfig:
    return replace(cfg, physics=replace(cfg.physics, nuclear=True, nuclear_table_id=tid, **phys))


def test_physics_options_defaults_and_validation(make_config: MakeConfig) -> None:
    ph = make_config().physics
    assert ph.e_cut_deuteron_mev == 4.0 and ph.nuclear_table_id is None
    with pytest.raises(UnsupportedCombinationError):
        replace(ph, e_cut_deuteron_mev=0.0)
    with pytest.raises(UnsupportedCombinationError):
        replace(ph, nuclear_table_id=3)  # type: ignore[arg-type]


def test_c1_nuclear_configuration_cases_raise_before_transport(
    make_config: MakeConfig,
    tid: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ok = _nuc(make_config(energy=100.0), tid)
    validate(ok)  # the reference case is valid
    # the Warp backends accept nuclear=True since V3-005B (C11); the cases were rejected before
    validate(_nuc(make_config(energy=100.0, backend="warp-cpu", precision="float32"), tid))
    cases: dict[str, tuple[SimulationConfig, str]] = {
        "E0 + 6 sigma_E above 250 MeV": (
            _nuc(make_config(energy=245.0, energy_sigma=1.0), tid),
            "250",
        ),
        "E0 above 250 MeV": (_nuc(make_config(energy=260.0), tid), "250"),
        "no table id": (_nuc(make_config(energy=100.0), None), "nuclear_table_id"),
        "unsupported element": (
            _nuc(
                make_config(
                    energy=100.0,
                    geometry=BoxPhantom((-30.0, -30.0, 0.0), (60, 60, 60), LEAD),
                ),
                tid,
            ),
            "Pb",
        ),
        "unproducible species (alpha)": (
            replace(
                ok,
                scoring=ok.scoring,
                tallies=(TallyRequest("t", ok.scoring[0].name, "dose", species=("alpha",)),),
            ),
            "producible",
        ),
        "unproducible species (triton)": (
            replace(
                ok,
                tallies=(TallyRequest("t", ok.scoring[0].name, "fluence", species=("triton",)),),
            ),
            "producible",
        ),
    }
    for label, (cfg, text) in cases.items():
        with pytest.raises(UnsupportedCombinationError, match=text):
            validate(cfg)
            pytest.fail(f"{label}: no exception")
        with pytest.raises(UnsupportedCombinationError):
            Simulation(cfg)  # nothing runs: construction already fails
    # non-proton source
    d_src = replace(ok.source, projectile=DEUTERON)
    with pytest.raises(UnsupportedCombinationError):
        validate(replace(ok, source=d_src))
    # nist-star stopping (explicit message)
    nist = replace(ok, physics=replace(ok.physics, stopping=_NistStarStandIn()))
    with pytest.raises(UnsupportedCombinationError, match="deuteron"):
        validate(nist)
    # missing table
    with pytest.raises(NuclearTableMissingError) as miss:
        validate(_nuc(make_config(energy=100.0), "0" * 64))
    assert isinstance(miss.value, UnsupportedCombinationError)
    # stale table (npz bytes changed) in a tmp cache
    d = tmp_path / "derived"
    d.mkdir()
    src_dir = cache.resolve_cache_dir(None) / "derived"
    for suffix in (".npz", ".json"):
        shutil.copy(src_dir / f"nuclear-proton-{tid}{suffix}", d / f"nuclear-proton-{tid}{suffix}")
    npz = d / f"nuclear-proton-{tid}.npz"
    raw = bytearray(npz.read_bytes())
    raw[len(raw) // 2] ^= 0xFF
    npz.write_bytes(bytes(raw))
    with monkeypatch.context() as mp:
        mp.setenv(cache.CACHE_ENV, str(tmp_path))
        with pytest.raises(NuclearTableStaleError):
            validate(ok)
    # mis-pinned source
    pin = DATASETS["ame2020-mass"]
    with monkeypatch.context() as mp:
        mp.setitem(DATASETS, "ame2020-mass", type(pin)(**{**pin.__dict__, "sha256": "f" * 64}))
        with pytest.raises(NuclearTablePinError):
            validate(ok)


def test_nuclear_effective_config_and_capabilities(make_config: MakeConfig, tid: str) -> None:
    cfg = _nuc(make_config(energy=100.0), tid)
    eff = validate(cfg)
    assert eff.nuclear is not None and eff.tables.deuteron is not None
    s = eff.summary()
    n = s["nuclear"]
    info = NuclearTable.load(None, tid).info
    assert n["table_id"] == tid and n["npz_sha256"] == info["npz_sha256"]
    assert n["source_sha256"] == {k: v["sha256"] for k, v in sorted(info["sources"].items())}
    assert set(n["elements"]) == {"H", "O"} & set(info["elements"]) or set(n["elements"]) == {"O"}
    assert n["elements"]["O"]["surrogate"] is False
    assert n["q_plus_table_mev"] == info["q_plus_table_mev"]
    assert n["transport_energy_bound_mev"] == info["transport_energy_bound_mev"]
    assert n["history_energy_bound_mev"] == info["history_energy_bound_mev"]
    assert n["grid_size"] == info["grid"]["n_points"] == 604
    assert n["multiplicity_model"] == info["multiplicity"]["model"]
    assert n["majorant_factor"] == 1.02
    assert n["e_cut_deuteron_mev"] == 4.0 == s["physics"]["e_cut_deuteron_mev"]
    assert {(p["species"], p["generation"]) for p in n["producible"]} == producible(PROTON, True)
    assert n["deuteron_tables_sha256"] == eff.tables.deuteron.sha256
    json.dumps(s)  # JSON-serialisable
    # nuclear=False: no nuclear section, summary keys unchanged
    off = validate(make_config(energy=100.0)).summary()
    assert "nuclear" not in off and "e_cut_deuteron_mev" not in off["physics"]
    # capability report: nuclear section only on request; the default report is the pre-V3-005A one
    cap = ionmc.capabilities(nuclear=True)
    assert (
        cap["nuclear"]["backends"] == ["python", "warp-cpu", "warp-cuda"]
        and "nuclear" not in ionmc.capabilities()
    )
    assert cap["nuclear"]["deuteron_cutoff_mev_default"] == cfg.physics.e_cut_deuteron_mev
    assert cap["nuclear"]["max_particles_per_history"] == MAX_PARTICLES
    base = ionmc.capabilities()
    base["backends"].pop("warp-cuda")  # host dependent
    digest = hashlib.sha256(json.dumps(base, sort_keys=True).encode()).hexdigest()
    assert digest == CAPABILITIES_SHA256  # the report of 9913ddc6 (nuclear=False)
    # the tally section is consistent with the nuclear one (secondaries transported and accepted)
    assert {(p["species"], p["generation"]) for p in cap["tallies"]["producible"]} == producible(
        PROTON, True
    )
    assert "secondary" in cap["tallies"]["generations"]["accepted"]
    assert "secondary" not in cap["tallies"]["generations"]["rejected"]
    assert "not transported yet" not in cap["tallies"]["generations"]["note"]
    # nuclear=True: physics.nuclear true, tally backends python only with warp deferred, the
    # fail-closed list consistent with the accepted secondaries (Codex finding 4)
    assert cap["physics"]["nuclear"] is True and base["physics"]["nuclear"] is False
    assert list(cap["tallies"]["backends"]) == ["python", "warp-cpu", "warp-cuda"]
    assert cap["tallies"]["untested_in_sandbox"] == ["warp-cuda"]  # CUDA: no GPU in the sandbox
    assert "untested_in_sandbox" not in base["tallies"]
    fc = cap["tallies"]["fail_closed"]
    assert fc[0] == "unknown or unproducible species"  # no sentence lists 'secondary' as an error
    assert not any("secondary" in e or "not an error" in e for e in fc)  # accepted only above
    assert any("secondary" in e for e in base["tallies"]["fail_closed"])  # nuclear=False: pinned
    assert base["tallies"]["fail_closed"][0] == (
        "unknown or unproducible species, generation 'secondary'"
    )
    del cap["nuclear"]
    cap["backends"].pop("warp-cuda")
    assert cap["tallies"] != base["tallies"]
    cap["tallies"] = base["tallies"]
    cap["physics"] = base["physics"]
    assert cap == base


CAPABILITIES_SHA256 = "a3d3052e5a894302a5fcd90b63d710d49fe8b8191dbf7055e517004ad279472e"


def test_capacity_bounds_follow_decision_0041_section_5(
    make_config: MakeConfig, tid: str, bethe_src: BetheStoppingSource
) -> None:
    base = make_config(energy=100.0)
    grid = base.scoring[0].name
    tallies = (
        TallyRequest("let", grid, "let_d"),
        TallyRequest("dose", grid, "dose"),
        TallyRequest("fl", grid, "fluence"),
    )
    off = validate(replace(base, tallies=tallies)).channels
    on_cfg = replace(_nuc(base, tid), tallies=tallies)
    eff = validate(on_cfg)
    assert eff.channels is not None and off is not None
    b = eff.channels.bounds
    info = NuclearTable.load(None, tid).info
    t_bound = info["transport_energy_bound_mev"]
    h_bound = info["history_energy_bound_mev"]
    terms = info["transport_path_bound_terms"]
    assert set(terms) == {"n", "p", "d", "a", "g"}
    assert (
        h_bound
        == sum(v["n_max"] * v["t_lab_max_mev"] for v in terms.values()) + info["recoil_t_max_mev"]
    )
    assert h_bound >= t_bound and h_bound > 2.0 * t_bound
    e_hi = 100.0
    rho = {0: float(eff.geometry.densities_g_cm3().min())}
    mpb_p = mixed_path_bound_mm(eff.tables, rho, e_hi)
    expect = mpb_p
    assert eff.tables.deuteron is not None
    for key, tab in (("p", eff.tables), ("d", eff.tables.deuteron)):  # transported species only
        expect += terms[key]["n_max"] * mixed_path_bound_mm(tab, rho, terms[key]["t_lab_max_mev"])
    assert eff.channels.path_bound_mm == pytest.approx(STRAGGLING_MARGIN * expect, rel=1e-12)
    assert b["E_bound_mev"] == max(e_hi, h_bound) and b["E_hi_mev"] == e_hi
    assert eff.channels.path_bound_mm > off.path_bound_mm
    assert b["max_particles"] == MAX_PARTICLES
    # the piece-count capacity is 32 times larger and the E bound is the larger one
    n_ch = [c for c in eff.channels.channels if c.kind == "N"][0]
    n_off = [c for c in off.channels if c.kind == "N"][0]
    assert n_ch.bound_per_history == pytest.approx(MAX_PARTICLES * n_off.bound_per_history)
    e_ch = [c for c in eff.channels.channels if c.kind == "E"][0]
    assert e_ch.bound_per_history == max(e_hi, h_bound)
    # a nuclear=False run is untouched by the new code path
    assert "E_bound_mev" not in off.bounds


def test_bounds_at_250_mev_in_water(make_config: MakeConfig, tid: str) -> None:
    """The recorded numbers: B_L and the E bound for 250 MeV protons in a water box."""
    base = make_config(energy=250.0, n=2, n_batches=2)
    grid = base.scoring[0].name
    cfg = replace(_nuc(base, tid), tallies=(TallyRequest("let", grid, "let_d"),))
    ch = validate(cfg).channels
    assert ch is not None
    print(
        f"250 MeV water: B_L = {ch.path_bound_mm:.6f} mm, E bound = {ch.bounds['E_bound_mev']:.6f} "
        f"MeV, terms {ch.bounds}"
    )
    assert ch.bounds["E_bound_mev"] >= 250.0
