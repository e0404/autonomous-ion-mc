"""V3-004 (decision 0040): species registry, tally requests, lookup tables, channel compiler and
the fail-closed rules of acceptance row A12 (plan validation/plans/v3-004-acceptance.md)."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import ionmc.config as config_module
from ionmc.config import SimulationConfig, validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.lookup import LookupTable, resample_uniform
from ionmc.materials import AIR
from ionmc.physics.projectiles import DEUTERON, PROTON
from ionmc.scoring import TALLY_QUANTITIES, ScoringGrid, TallyRequest
from ionmc.species import (
    N_SPECIES_SLOTS,
    REGISTRY,
    producible,
    species_by_id,
    species_by_name,
    species_of_projectile,
    transported_ids,
)
from ionmc.transport.channels import (
    CLASS_LOCAL,
    CLASS_STEP,
    EXCLUDED_CHANNEL_NAME,
    LOCAL_PIECE_COUNT_NAME,
    MAX_CHANNELS,
    PIECE_COUNT_NAME,
    ChannelPlan,
    adaptive_exponent,
    fe_exponent,
)

MakeConfig = Callable[..., SimulationConfig]
SYNTHETIC = Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"
LOG_EDGES = (2.0, 4.0, 8.0, 16.0)
LIN_EDGES = (1.0, 2.0, 3.0)
GRID = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 30), name="dose")


@pytest.fixture
def channels_enabled(monkeypatch: pytest.MonkeyPatch) -> None:
    """All backends implement channels (kept so the tests state their dependence explicitly)."""
    monkeypatch.setattr(
        config_module, "CHANNEL_BACKENDS", ("python", "warp-cpu", "warp-cuda"), raising=True
    )


def _plan(cfg: SimulationConfig) -> ChannelPlan:
    plan = validate(cfg).channels
    assert plan is not None
    return plan


def _with(cfg: SimulationConfig, *tallies: TallyRequest, **kw: Any) -> SimulationConfig:
    return replace(cfg, tallies=tuple(tallies), **kw)


def _req(name: str, quantity: str, **kw: Any) -> TallyRequest:
    return TallyRequest(name, "dose", quantity, **kw)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- species registry


def test_registry_ids_are_frozen_and_append_only() -> None:
    assert [(s.id, s.name) for s in REGISTRY if s.transported] == [
        (0, "proton"),
        (1, "deuteron"),
        (2, "triton"),
        (3, "helium3"),
        (4, "alpha"),
    ]
    assert species_by_name("nuclear_local").transported is False
    assert species_by_name("nuclear_local").id == 64 < N_SPECIES_SLOTS
    assert species_by_id(0).name == "proton"
    assert transported_ids() == (0, 1, 2, 3, 4)
    with pytest.raises(ValueError):
        species_by_name("carbon12")  # not appended yet (V3-009)
    assert species_of_projectile(PROTON).id == 0


def test_producible_is_the_source_primary_only() -> None:
    assert producible(PROTON) == frozenset({("proton", "primary")})
    assert producible(DEUTERON) == frozenset({("deuteron", "primary")})


# --------------------------------------------------------------------------- TallyRequest


def test_tally_request_rules() -> None:
    assert _req("a", "let_d").species is None
    for kwargs in (
        {"quantity": "bogus"},
        {"quantity": "let_d", "generation": "tertiary"},
        {"quantity": "let_d", "species": ()},
        {"quantity": "let_d", "species": ("proton", "proton")},
        {"quantity": "let_d", "lookup": "x"},  # lookup only for lookup_*
        {"quantity": "lookup_sum"},  # needs a lookup
        {"quantity": "fluence_spectrum"},  # needs edges
        {"quantity": "let_d", "energy_edges_mev_per_u": (1.0, 2.0)},
        {"quantity": "fluence_spectrum", "energy_edges_mev_per_u": (2.0, 1.0)},
        {"quantity": "fluence_spectrum", "energy_edges_mev_per_u": (1.0,)},
    ):
        with pytest.raises(UnsupportedCombinationError):
            TallyRequest("t", "dose", **kwargs)


# --------------------------------------------------------------------------- LookupTable


def _doc(**over: Any) -> dict[str, Any]:
    doc = json.loads(SYNTHETIC.read_text())
    doc.update(over)
    return doc


def _write(tmp_path: Path, doc: dict[str, Any]) -> Path:
    p = tmp_path / "t.json"
    p.write_text(json.dumps(doc))
    return p


def test_synthetic_fixture_loads_with_provenance() -> None:
    raw = SYNTHETIC.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    lk = LookupTable.from_file(SYNTHETIC, expected_sha256=sha)
    assert lk.synthetic is True and lk.file_sha256 == sha
    assert lk.axis == "let_water_kev_um" and lk.axis_range == (0.0, 100.0)
    x = np.array([0.0, 2.5, 100.0, 120.0])
    f, inside = lk.evaluate("proton", x)
    np.testing.assert_allclose(f[:3], 1.0 + 0.1 * x[:3], rtol=1e-15)
    assert inside.tolist() == [True, True, True, False]
    prov = lk.provenance()
    assert prov["file_sha256"] == sha and prov["synthetic"] is True
    assert prov["content_sha256"] == lk.content_sha256 != sha
    assert prov["resampling"] is None and prov["license"] and prov["citation"]
    json.dumps(prov)  # serialisable
    with pytest.raises(UnsupportedCombinationError, match="sha256"):
        LookupTable.from_file(SYNTHETIC, expected_sha256="0" * 64)
    with pytest.raises((ValueError, TypeError)):
        lk.axis_values[0] = 1.0  # read-only


@pytest.mark.parametrize(
    "mutate, match",
    [
        (
            lambda d: d.update(
                axis_values=[0.0, 1.0, 3.0] + [4.0] * 0, species={"proton": [1, 2, 3]}
            ),
            "uniform",
        ),
        (lambda d: d["species"].update(proton=[-1.0] + d["species"]["proton"][1:]), "non-negative"),
        (
            lambda d: d["species"].update(proton=[float("nan")] + d["species"]["proton"][1:]),
            "finite",
        ),
        (lambda d: d["species"].update(proton=[1.0, 2.0]), "values for"),
        (lambda d: d["species"].update(unobtainium=d["species"]["proton"]), "unknown species"),
        (lambda d: d.pop("citation"), "required fields"),
        (lambda d: d.pop("synthetic"), "required fields"),
        (lambda d: d.update(license=""), "non-empty"),
        (lambda d: d.update(axis="depth_mm"), "axis must be one of"),
        (lambda d: d.update(axis_spacing="log"), "positive"),  # axis starts at 0
    ],
)
def test_lookup_file_rejections(
    tmp_path: Path, mutate: Callable[[dict[str, Any]], Any], match: str
) -> None:
    doc = _doc()
    mutate(doc)
    with pytest.raises(UnsupportedCombinationError, match=match):
        LookupTable.from_file(_write(tmp_path, doc))


def test_lookup_log_axis_uniform_in_ln_and_resampling_is_recorded() -> None:
    x = np.array([1.0, 2.0, 5.0, 10.0, 40.0])  # non-uniform source data
    v = {"proton": 2.0 * x}
    with pytest.raises(UnsupportedCombinationError, match="uniform"):
        LookupTable("t", "q", "1", "energy_per_nucleon_mev", "linear", x, v, "c", "l", "s", False)
    new, vals, record = resample_uniform(x, v, n_points=9, spacing="log")
    lk = LookupTable(
        "t", "q", "1", "energy_per_nucleon_mev", "log", new, vals, "c", "l", "s", False,
        resampling=record,
    )  # fmt: skip
    assert lk.axis_range == (1.0, 40.0)
    assert lk.provenance()["resampling"]["n_points_in"] == 5
    assert lk.provenance()["resampling"]["source_axis_sha256"]
    # linear interpolation in x of a linear function is exact
    np.testing.assert_allclose(vals["proton"], 2.0 * new, rtol=1e-12)


# --------------------------------------------------------------------------- compiler


def test_compile_dedup_offsets_and_automatic_channels(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg = _with(
        make_config(energy=100.0, n=400, n_batches=20),
        _req("let_t", "let_t"),
        _req("let_d", "let_d"),
        _req("dose", "dose"),
        _req("edep", "edep"),  # same channel as dose
        _req("flu", "fluence"),  # same L channel as the denominator of let_t
        _req("let_e", "let_d_eps"),
    )
    plan = _plan(cfg)
    kinds = [c.kind for c in plan.channels]
    # L (shared by let_t and fluence), LS (shared by let_t and let_d), LS2, ES, E_step, E (edep =
    # dose), N (step), N_local, and the excluded-energy channel: 9 channels, no duplicates
    assert sorted(kinds) == sorted(["L", "LS", "LS2", "ES", "E", "E", "N", "N", "E"])
    assert len(plan.channels) == 9
    n_step = plan.channels[plan.count_channel(0, CLASS_STEP)]
    n_loc = plan.channels[plan.count_channel(0, CLASS_LOCAL)]
    assert n_step.k == n_loc.k == 0 and n_step.residual_column == n_loc.residual_column == -1
    assert (
        plan.species_match[plan.count_channel(0, CLASS_LOCAL), species_by_name("nuclear_local").id]
        == 1
    )
    q = {d.name: d for d in plan.quantities}
    assert q["dose"].numerator == q["edep"].numerator
    assert q["flu"].numerator == q["let_t"].denominator
    assert q["let_t"].numerator == q["let_d"].denominator  # LS
    assert q["let_d"].kind == "ratio" and q["flu"].kind == "linear"
    # edep covers both classes and the pseudo-species; E_step only steps; excluded only local
    e_edep = plan.channels[q["edep"].numerator]
    e_step = plan.channels[q["let_e"].denominator]
    assert e_edep.class_mask == CLASS_STEP | CLASS_LOCAL and e_step.class_mask == CLASS_STEP
    local = [c for c in plan.channels if c.class_mask == CLASS_LOCAL and c.kind == "E"]
    assert len(local) == 1
    assert (
        plan.species_match[plan.channels.index(local[0]), species_by_name("nuclear_local").id] == 1
    )
    assert plan.species_match[plan.channels.index(e_step), species_by_name("nuclear_local").id] == 0
    assert plan.species_match.dtype == np.int8
    assert plan.species_match.shape == (9, N_SPECIES_SLOTS)
    # offsets are contiguous, sizes are the grid voxels, residual column for every channel but N
    off = 0
    for c in plan.channels:
        assert c.offset == off and c.size == GRID.n_voxels
        off += c.size
    assert plan.total_size == off
    assert plan.ch_begin == (0,) and plan.ch_end == (9,)
    res = sorted(c.residual_column for c in plan.channels)
    assert res == [-1, -1, 0, 1, 2, 3, 4, 5, 6] and plan.n_residual == 7
    assert EXCLUDED_CHANNEL_NAME == "edep_excluded_from_let"
    assert plan.memory_bytes == 20 * (GRID.n_voxels + plan.total_size) * 8
    json.dumps(validate(cfg).summary())  # the effective configuration stays serialisable


def test_compile_orders_channels_by_grid_and_spectrum_sizes(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    g2 = ScoringGrid((-30.0, -30.0, 0.0), (5.0, 5.0, 5.0), (6, 6, 12), name="coarse")
    cfg = make_config(energy=100.0, n=400, n_batches=20, scoring=(GRID, g2))
    reqs = (
        TallyRequest("c", "coarse", "let_t"),
        TallyRequest("f", "dose", "fluence"),
        TallyRequest("spec", "coarse", "fluence_spectrum", energy_edges_mev_per_u=LOG_EDGES),
        TallyRequest("lin", "coarse", "fluence_spectrum", energy_edges_mev_per_u=LIN_EDGES),
    )
    plan = _plan(_with(cfg, *reqs))
    assert [c.grid for c in plan.channels] == sorted(c.grid for c in plan.channels)
    assert plan.ch_begin[0] == 0 and plan.ch_end[0] == plan.ch_begin[1]
    assert plan.ch_end[1] == len(plan.channels)
    fl = [c for c in plan.channels if c.kind == "FL"]
    sizes = sorted(c.size for c in fl)
    assert sizes == [g2.n_voxels * 4, g2.n_voxels * 5]  # n_bins + 2
    assert {c.spectrum.log for c in fl if c.spectrum} == {True, False}
    l_chan = next(c for c in plan.channels if c.kind == "L")
    assert all(c.k == l_chan.k for c in fl)  # spectra share the exponent of L


def test_species_filter_gives_separate_channels(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    plan = _plan(
        _with(
            cfg,
            _req("all", "fluence"),
            _req("p", "fluence", species=("proton",)),
            _req("prim", "fluence", generation="primary"),
        )
    )
    ls = [c for c in plan.channels if c.kind == "L"]
    assert len(ls) == 3
    prim = [c for c in ls if (c.gen_lo, c.gen_hi) == (0, 0)]
    assert len(prim) == 1
    only_p = [c for c in ls if sum(c.species) == 1]
    assert len(only_p) == 1 and only_p[0].species[species_by_name("proton").id] == 1


# --------------------------------------------------------------------------- quanta


@pytest.mark.parametrize(
    "hpb, kl, kls, kls2",
    [(10**6, 34, 33, 28), (10**7, 31, 30, 24), (28_600_000, 29, 28, 23)],
)
def test_worked_example_quanta(
    make_config: MakeConfig, channels_enabled: None, hpb: int, kl: int, kls: int, kls2: int
) -> None:
    """Delta section (ii): water, protons, 150 MeV, E_cut 2 MeV (ramp envelope
    S_bar_max 45.1 MeV/mm, B_L about 197 mm, B_LS about 339 MeV, B_LS2 about 1.53e4)."""
    cfg = make_config(energy=150.0, n=20 * hpb, n_batches=20, e_cut=2.0)
    plan = _plan(_with(cfg, _req("ld", "let_d"), _req("let", "let_d_eps"), _req("flu", "fluence")))
    k = {c.kind: c.k for c in plan.channels}
    assert k["L"] == kl and k["LS"] == kls and k["LS2"] == kls2
    assert k["E"] == 30 and k["N"] == 0
    b = plan.bounds
    s_bar_max = b["S_w_max_mev_per_mm"]  # ramp envelope, E_mid >= E_cut/2
    assert s_bar_max == pytest.approx(45.1, rel=0.02)
    assert 0.0 < b["S_w_min_mev_per_mm"] < 0.5 and b["gamma_max"] == pytest.approx(0.81, rel=0.02)
    assert b["B_L_mm"] == pytest.approx(197.0, rel=0.02)
    assert b["r_max"] == pytest.approx(1.81, rel=0.01)  # (1 + |gamma|) over [E_cut/2, E_hi]
    assert b["B_LS_mev"] == pytest.approx(1.25 * b["r_max"] * 150.0, rel=1e-9)
    assert b["B_LS2"] == pytest.approx(1.53e4, rel=0.02)
    # capacity: hpb * B_c * 2^k < 2^62 for every adaptive kind, and one more bit would not fit
    for c in plan.channels:
        if c.kind in ("L", "LS", "LS2", "ES"):
            assert plan.histories_per_batch * c.bound_per_history * 2.0**c.k < 2.0**62


def test_adaptive_exponent_rule_and_determinism() -> None:
    assert adaptive_exponent(10**6, 197.0) == 34
    assert adaptive_exponent(1, 1.0) == 40  # capped
    assert adaptive_exponent(2**10, 2.0**20) == 31  # 2^62 itself is excluded
    assert adaptive_exponent(10**6, 197.0) == adaptive_exponent(10**6, 197.0)
    # k is monotone decreasing in hpb*B and never lets the total reach 2^62
    for hpb in (3, 10**3, 10**9):
        for b in (0.5, 7.0, 1e5):
            k = adaptive_exponent(hpb, b)
            assert hpb * b * 2.0**k < 2.0**62
            assert k == 40 or hpb * b * 2.0 ** (k + 1) >= 2.0**62
    with pytest.raises(UnsupportedCombinationError):
        adaptive_exponent(10, 0.0)
    assert fe_exponent(1.0) == 30 and fe_exponent(2.0) == 29 and fe_exponent(1.5) == 29
    assert fe_exponent(0.25) == 32


def test_quanta_are_a_function_of_config_only(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg = _with(make_config(energy=100.0, n=2000, n_batches=20), _req("ld", "let_d"))
    a, b = _plan(cfg), _plan(cfg)
    assert [c.k for c in a.channels] == [c.k for c in b.channels]
    assert a.summary() == b.summary()


def test_precision_floor_fails_closed(make_config: MakeConfig, channels_enabled: None) -> None:
    """In a low-density medium the per-history bounds r_max and B_L grow by the density ratio, so
    at a high hpb the adaptive quantum exceeds 2^-16 of a 1 mm entrance piece: fail closed."""
    air = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 60.0), AIR)
    fine = make_config(energy=100.0, n=20 * 10**5, n_batches=20, geometry=air)
    _plan(_with(fine, _req("ld", "let_d")))  # a small hpb fits the floor
    heavy = make_config(energy=100.0, n=20 * 4 * 10**7, n_batches=20, geometry=air)
    with pytest.raises(UnsupportedCombinationError, match="precision floor"):
        validate(_with(heavy, _req("ld", "let_d")))
    # the same hpb is fine in water (the E capacity limit is 4.3e7 histories per batch at 100 MeV)
    ok = make_config(energy=100.0, n=20 * 4 * 10**7, n_batches=20)
    _plan(_with(ok, _req("ld", "let_d")))


def test_energy_capacity_rule_still_applies_with_tallies(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg2 = _with(make_config(energy=150.0, n=20 * 10**8, n_batches=20), _req("ld", "let_d"))
    with pytest.raises(UnsupportedCombinationError, match="capacity"):
        validate(cfg2)


def test_memory_guard_counts_channels(make_config: MakeConfig, channels_enabled: None) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    base = _with(cfg, _req("ld", "let_d"))
    plan = _plan(base)
    edep_only = 20 * GRID.n_voxels * 8
    assert plan.memory_bytes > edep_only
    tight = replace(base, run=replace(base.run, memory_budget_bytes=plan.memory_bytes - 1))
    with pytest.raises(UnsupportedCombinationError, match="memory budget"):
        validate(tight)
    ok = replace(base, run=replace(base.run, memory_budget_bytes=plan.memory_bytes))
    validate(ok)
    two_workers = replace(
        base, run=replace(base.run, cpu_workers=2, memory_budget_bytes=plan.memory_bytes)
    )
    with pytest.raises(UnsupportedCombinationError, match="memory budget"):
        validate(two_workers)


# --------------------------------------------------------------------------- A12


def test_no_tallies_is_unchanged_and_tallies_are_never_ignored(
    make_config: MakeConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = make_config(energy=20.0)
    assert validate(cfg).channels is None
    assert validate(cfg).summary()["tallies"] is None
    # every backend has channels now; a backend that lacks them (simulated) rejects the request
    monkeypatch.setattr(config_module, "CHANNEL_BACKENDS", ("python",), raising=True)
    for backend, prec in (("warp-cpu", "float32"), ("warp-cpu", "float64")):
        c = _with(make_config(energy=20.0, backend=backend, precision=prec), _req("ld", "let_d"))
        with pytest.raises(UnsupportedCombinationError, match="does not implement scoring"):
            validate(c)


def _lookup(
    species: str = "proton",
    lo: float = 0.0,
    hi: float = 100.0,
    axis: str = "let_water_kev_um",
    spacing: str = "linear",
) -> LookupTable:
    x = np.linspace(lo, hi, 11) if spacing == "linear" else np.geomspace(lo, hi, 11)
    return LookupTable(
        "lk", "q", "1", axis, spacing, x, {species: 1.0 + 0.1 * x}, "c", "l", "s", True
    )


_EN = "energy_per_nucleon_mev"
_LK = {"lookup": "lk"}
_BAD_EDGES = {"energy_edges_mev_per_u": (1.0, 2.0, 5.0, 6.0)}
A12_CASES = [
    ("unknown species", "let_d", {"species": ("unobtainium",)}, [], "unknown species"),
    ("unproducible species", "let_d", {"species": ("alpha",)}, [], "not producible"),
    ("unproducible deuteron", "fluence", {"species": ("deuteron",)}, [], "not producible"),
    ("secondary generation", "let_d", {"generation": "secondary"}, [], "not producible"),
    (
        "secondary of p",
        "fluence",
        {"species": ("proton",), "generation": "secondary"},
        [],
        "not producible",
    ),
    ("let_medium local", "let_d", {"let_medium": "local"}, [], "let_medium"),
    ("dose to water", "dose", {"dose_reference": "water"}, [], "dose_reference"),
    ("unknown grid", "let_d", {"grid": "nogrid"}, [], "unknown scoring grid"),
    ("lookup missing", "lookup_sum", {"lookup": "nope"}, [], "unknown lookup"),
    ("lookup species gap", "lookup_sum", _LK, [_lookup("deuteron")], "no values for"),
    ("lookup coverage let", "lookup_sum", _LK, [_lookup(hi=5.0)], "covers"),
    ("lookup coverage energy", "lookup_sum", _LK, [_lookup(lo=3.0, hi=60.0, axis=_EN)], "covers"),
    ("unused lookup", "let_d", {}, [_lookup()], "not used"),
    ("non-uniform spectrum", "fluence_spectrum", _BAD_EDGES, [], "uniform"),
    ("pseudo species LET", "let_d", {"species": ("nuclear_local",)}, [], "not transported"),
    # review 67f03e06 (a): an explicit nuclear_local edep/dose request is not producible yet
    ("nuclear_local edep", "edep", {"species": ("nuclear_local",)}, [], "not producible"),
    ("nuclear_local dose", "dose", {"species": ("nuclear_local",)}, [], "not producible"),
]


@pytest.mark.parametrize("name, quantity, kwargs, lookups, match", A12_CASES)
def test_a12_fail_closed(
    make_config: MakeConfig,
    channels_enabled: None,
    name: str,
    quantity: str,
    kwargs: dict[str, Any],
    lookups: list[LookupTable],
    match: str,
) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    kw = dict(kwargs)
    grid = kw.pop("grid", "dose")
    req = TallyRequest("x", grid, quantity, **kw)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError, match=match):
        validate(_with(cfg, req, lookups=tuple(lookups)))


def test_a12_duplicate_request_names(make_config: MakeConfig, channels_enabled: None) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    with pytest.raises(UnsupportedCombinationError, match="unique"):
        validate(_with(cfg, _req("x", "let_d"), _req("x", "let_t")))


def test_a12_valid_lookup_is_accepted(make_config: MakeConfig, channels_enabled: None) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    lk = LookupTable.from_file(SYNTHETIC)
    # LET axis must cover the water-S interval of 1..100 MeV (up to about 26 keV/um): ok (0..100)
    plan = _plan(
        _with(
            cfg,
            _req("lk", "lookup_dose_avg", lookup="synthetic_let_linear"),
            _req("lks", "lookup_sum", lookup="synthetic_let_linear", species=("proton",)),
            lookups=(lk,),
        )
    )
    fe = [c for c in plan.channels if c.kind == "FE"]
    assert len(fe) == 2 and all(c.lookup == 0 for c in fe)
    assert fe[0].k == fe[1].k == fe_exponent(float(lk.values["proton"].max()))
    assert plan.summary()["lookups"][0]["synthetic"] is True


def test_a12_energy_axis_coverage_uses_per_nucleon_range(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    ok = _lookup(lo=0.9, hi=120.0, axis="energy_per_nucleon_mev", spacing="log")
    _plan(_with(cfg, _req("x", "lookup_sum", lookup="lk"), lookups=(ok,)))
    short = _lookup(lo=0.9, hi=99.0, axis="energy_per_nucleon_mev", spacing="log")
    with pytest.raises(UnsupportedCombinationError, match="covers"):
        validate(_with(cfg, _req("x", "lookup_sum", lookup="lk"), lookups=(short,)))


def test_a12_lookup_provenance_in_effective_summary(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    lk = LookupTable.from_file(SYNTHETIC)
    eff = validate(_with(cfg, _req("x", "lookup_sum", lookup=lk.name), lookups=(lk,)))
    prov = eff.summary()["tallies"]["lookups"][0]
    assert prov["file_sha256"] == hashlib.sha256(SYNTHETIC.read_bytes()).hexdigest()
    assert prov["citation"] and prov["license"] and prov["source"]
    assert not math.isnan(sum(c["k"] for c in eff.summary()["tallies"]["channels"]))


def test_fe_quantum_bound_and_floor_are_per_channel(
    make_config: MakeConfig, channels_enabled: None
) -> None:
    """Review 67f03e06 (b): each FE channel uses the largest value of its own table (and species),
    not one global f_max over every table."""
    x = np.linspace(0.0, 100.0, 11)
    small = LookupTable("small", "q", "1", "let_water_kev_um", "linear", x,
                        {"proton": 1.0e-3 * (1.0 + x)}, "c", "l", "s", True)  # fmt: skip
    big = LookupTable("big", "q", "1", "let_water_kev_um", "linear", x,
                      {"proton": 1.0e3 * (1.0 + x)}, "c", "l", "s", True)  # fmt: skip
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    plan = _plan(
        _with(
            cfg,
            _req("a", "lookup_sum", lookup="small"),
            _req("b", "lookup_sum", lookup="big"),
            lookups=(small, big),
        )
    )
    fe = {plan.lookups[c.lookup].name: c for c in plan.channels if c.kind == "FE"}
    f_small = float(small.values["proton"].max())
    f_big = float(big.values["proton"].max())
    assert fe["small"].k == fe_exponent(f_small) and fe["big"].k == fe_exponent(f_big)
    assert fe["small"].k > fe["big"].k + 15  # not one shared exponent
    e_hi = plan.bounds["E_hi_mev"]
    assert fe["small"].bound_per_history == pytest.approx(f_small * e_hi)
    assert fe["big"].bound_per_history == pytest.approx(f_big * e_hi)
    # alone, the small table gives the same quantum as next to the big one
    alone = _plan(_with(cfg, _req("a", "lookup_sum", lookup="small"), lookups=(small,)))
    assert [c.k for c in alone.channels if c.kind == "FE"] == [fe["small"].k]


def test_lookup_table_is_deeply_immutable() -> None:
    """Review 67f03e06 (c): the values mapping and the provenance structures are frozen, not only
    the arrays; the caller's inputs are not aliased."""
    x = np.linspace(1.0, 100.0, 5)
    src_values = {"proton": 1.0 + 0.1 * x}
    rec = {"method": "m", "nested": {"a": [1, 2]}}
    lk = LookupTable(
        "lk", "q", "1", "let_water_kev_um", "linear", x, src_values, "c", "l", "s", True,
        resampling=rec,
    )  # fmt: skip
    before = (lk.content_sha256, lk.provenance())
    with pytest.raises(TypeError):
        lk.values["proton"] = np.zeros(5)  # type: ignore[index]  # mapping replacement
    with pytest.raises(TypeError):
        lk.values["alpha"] = np.zeros(5)  # type: ignore[index]
    with pytest.raises(TypeError):
        del lk.values["proton"]  # type: ignore[attr-defined]
    with pytest.raises(AttributeError):
        lk.values.update({})  # type: ignore[attr-defined]
    with pytest.raises(ValueError):
        lk.values["proton"][0] = 5.0  # array element
    with pytest.raises(ValueError):
        lk.values["proton"].setflags(write=True)
    with pytest.raises(TypeError):
        lk.resampling["method"] = "x"  # type: ignore[index]
    with pytest.raises(TypeError):
        lk.resampling["nested"]["a"] = 3  # type: ignore[index]
    assert isinstance(lk.resampling["nested"]["a"], tuple)  # type: ignore[index]
    # inputs are copied, not aliased
    src_values["proton"][0] = 99.0
    src_values["alpha"] = np.zeros(5)
    rec["nested"]["a"].append(3)
    assert "alpha" not in lk.values and lk.values["proton"][0] == pytest.approx(1.1)
    assert lk._content_hash() == lk.content_sha256 == before[0]
    assert lk.provenance() == before[1]
    prov = lk.provenance()
    prov["resampling"]["nested"]["a"].append(4)  # provenance() hands out a fresh plain copy
    assert lk.provenance() == before[1]
    json.dumps(prov)


def test_capability_report_matches_the_enforced_scoring_contract() -> None:
    """The public capability report is derived from, and checked against, the objects that
    enforce the scoring contract (it cannot drift again)."""
    import typing

    import ionmc
    from ionmc.config import CHANNEL_BACKENDS
    from ionmc.lookup import AXES
    from ionmc.physics.projectiles import PROTON
    from ionmc.species import producible

    cap = ionmc.capabilities()
    t = cap["tallies"]
    literal = list(typing.get_args(typing.get_type_hints(TallyRequest)["quantity"]))
    assert t["quantities"] == literal == list(TALLY_QUANTITIES)
    assert set(t["backends"]) == set(CHANNEL_BACKENDS) <= set(cap["backend_names"])
    for b, entry in t["backends"].items():
        assert entry["quantities"] == literal, b
        assert isinstance(entry["available"], bool)
    assert t["backends"]["python"]["available"] is True
    assert {(d["species"], d["generation"]) for d in t["producible"]} == set(producible(PROTON))
    assert t["producible"] == [{"species": "proton", "generation": "primary"}]
    assert t["generations"]["accepted"] == ["all", "primary"]
    assert t["generations"]["rejected"] == ["secondary"]
    assert t["let_medium"] == ["water"] and t["dose_reference"] == ["medium"]
    assert t["lookup"]["axes"] == list(AXES) and t["lookup"]["uniform_axis_required"] is True
    from ionmc.transport.channels import FE_F_MAX_EXPONENT, FE_F_MIN_EXPONENT

    assert t["lookup"]["max_value_range"] == [2.0**FE_F_MIN_EXPONENT, 2.0**FE_F_MAX_EXPONENT]
    assert t["max_channels"] == MAX_CHANNELS
    assert set(t["automatic_channels"]) == {
        PIECE_COUNT_NAME, LOCAL_PIECE_COUNT_NAME, EXCLUDED_CHANNEL_NAME
    }  # fmt: skip
    assert t["fail_closed"] and any("secondary" in r for r in t["fail_closed"])
    assert cap["species"] == ["proton"]
    # the report is JSON-serialisable
    json.dumps(cap)


def _fe_plan_with_scale(make_config: MakeConfig, scale: float) -> ChannelPlan:
    x = np.linspace(0.0, 100.0, 11)
    lk = LookupTable("lk", "q", "1", "let_water_kev_um", "linear", x,
                     {"proton": scale * (1.0 + x) / 101.0}, "c", "l", "s", True)  # fmt: skip
    cfg = make_config(energy=100.0, n=400, n_batches=20)
    return _plan(_with(cfg, _req("a", "lookup_sum", lookup="lk"), lookups=(lk,)))


@pytest.mark.parametrize("scale", [1e-300, 1e-30, 2.0**-61, 2.0**61, 1e30, 1e300])
def test_fe_extreme_lookup_magnitudes_fail_closed_in_validate(
    make_config: MakeConfig, channels_enabled: None, scale: float
) -> None:
    """Review 63d3ef7b (1): a finite lookup whose largest value is outside [2^-60, 2^60] cannot
    give a representable quantum and scale; validate() rejects it (it never reaches a backend)."""
    with pytest.raises(UnsupportedCombinationError, match="admissible magnitude range"):
        _fe_plan_with_scale(make_config, scale)


@pytest.mark.parametrize("scale", [2.0**-60, 2.0**-40, 1e-3, 1.0, 1e6, 2.0**60])
def test_fe_admissible_boundary_values_give_representable_scales(
    make_config: MakeConfig, channels_enabled: None, scale: float
) -> None:
    """The boundary values are either accepted with a finite power-of-two scale and a finite
    bound, or rejected by another fail-closed rule (never by an overflow later)."""
    try:
        plan = _fe_plan_with_scale(make_config, scale)
    except UnsupportedCombinationError as exc:
        assert "admissible magnitude range" not in str(exc)
        return
    for c in plan.channels:
        assert -64 <= c.k <= 100
        assert math.isfinite(math.ldexp(1.0, c.k)) and math.isfinite(c.bound_per_history)
        assert math.isfinite(c.bound_per_history * math.ldexp(1.0, c.k))
