"""Model ``icru90-water-v1`` (decision 0044; V3-013 step S4).

Rows ID/CC-HASH, DOM-1, DOM-3, DOM-4, DOM-5, DOM-6, ROLE-1, PER-2, PER-3 (ICRU case) and the
consistency checks CC-NODE, CC-INT, CC-I2 and CC-SPAN (section A.5) of
``validation/plans/v3-013-acceptance.md``. The A.5 checks are gating for implementation
correctness and are never evidence. The tests that need the ICRU 90 object (and, for CC-I2 and the
committed aggregates, the PSTAR/ASTAR objects) find the cache like the other cache-dependent tests
(``cache.resolve_cache_dir()``: ``IONMC_CACHE_DIR`` or the default) and skip with a message when an
object is absent; ``IONMC_REQUIRE_DATA=1`` makes that a failure. The tests that need no data run on
synthetic power-law arrays. No table value is stored in Git.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.data import cache
from ionmc.data.icru90 import Icru90Water, load_icru90_water
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import PMMA, WATER, Material, water
from ionmc.physics import stopping_icru90, stopping_models
from ionmc.physics.projectiles import ALPHA, CARBON12, DEUTERON, OXYGEN16, PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.physics.stopping_evidence import STOPPING_CIRCULAR_ROWS, STOPPING_EVIDENCE
from ionmc.physics.stopping_icru90 import (
    Icru90WaterSource,
    build_icru90_source,
    icru90_start_range,
)
from ionmc.physics.stopping_models import STOPPING_MODELS, resolve, spec_of, spec_sha256
from ionmc.physics.stopping_selection import (
    check_window,
    select,
    table_model_id,
    validate_identity_block,
)
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource
from ionmc.transport.tables import TransportTables

MID = "icru90-water-v1"
ICRU_ID = "geant4-icru90-stopping-11.4.2"
PSTAR_ID, ASTAR_ID = "nist-pstar-water-2005", "nist-astar-water-2005"
ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "validation/scripts/stopping/v13_consistency.py"
RESULTS = ROOT / "validation/results/stopping/v3-013"
REL = 1.0e-12


def _need(ids: tuple[str, ...]) -> Path:
    d = cache.resolve_cache_dir()
    try:
        for ds in ids:
            cache.verify(ds, d)
    except (FileNotFoundError, cache.IntegrityError) as exc:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(f"IONMC_REQUIRE_DATA=1 but {ids} are not cached: {exc}")
        pytest.skip(f"{ids} not cached in {d}")
    return d


@pytest.fixture(scope="module")
def cdir() -> Path:
    return _need((ICRU_ID,))


@pytest.fixture(scope="module")
def cdir_all() -> Path:
    return _need((ICRU_ID, PSTAR_ID, ASTAR_ID))


@pytest.fixture(scope="module")
def gen() -> Any:
    spec = importlib.util.spec_from_file_location("v13_consistency", SCRIPT)
    assert spec is not None and spec.loader is not None
    m = importlib.util.module_from_spec(spec)
    sys.modules["v13_consistency"] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def data(cdir: Path) -> Icru90Water:
    return load_icru90_water(cache.verify(ICRU_ID, cdir))


def _synthetic() -> Icru90Water:
    """Power-law arrays S = c E_u^-p per projectile (no ICRU value): total energies span the
    windows of the model on both sides, the alpha energies are total (4 E_u)."""
    e_p = np.geomspace(0.1, 1000.0, 41)
    e_a = np.geomspace(0.4, 1000.0, 41)
    return Icru90Water(
        e_p,
        50.0 * e_p**-0.8,
        e_a,
        150.0 * (e_a / 4.0) ** -0.7,
        content_sha256="synthetic",
        dataset_id=None,
    )


def _cfg(proj: Any = PROTON, e_u: float | None = None, **physics: Any) -> SimulationConfig:
    e_u = e_u if e_u is not None else (20.0 if proj.z == 1 else 50.0)
    sigma = physics.pop("sigma_mev", 0.0)
    geo = physics.pop("geometry", None) or BoxPhantom(
        (-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER
    )
    physics.setdefault("stopping_model", MID)
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e_u * proj.a, sigma, 0.0),
        geometry=geo,
        scoring=(ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 2.0), (10, 10, 20)),),
        physics=PhysicsOptions(nuclear=physics.pop("nuclear", False), **physics),
        run=RunOptions(backend="python", precision="float64", seed=1, n_histories=2, n_batches=2),
    )


# ---- registered spec (cache-independent) ----


def test_registered_spec_is_cache_independent_and_carries_one_construction_pin() -> None:
    from ionmc.data.registry import DATASETS

    spec = STOPPING_MODELS[MID]
    assert spec.kind == "tabulated" and spec.version == "1" and spec.decision == "0044"
    (pin,) = spec.data
    ds = DATASETS[ICRU_ID]
    assert pin.dataset_id == ICRU_ID and pin.archive_sha256.startswith("afd14794")
    assert (pin.archive_sha256, pin.parser, pin.licence) == (ds.sha256, "icru90", ds.license)
    assert pin.licence == "Geant4 Software License"
    assert pin.member is None and pin.member_sha256 is None
    assert pin.use_role == "construction-permitted"
    assert ds.role == "evaluation"  # the dataset's own role is unchanged by the model
    assert spec.interpolation == "log-log of S onto 200/decade from 1 MeV/u"
    assert spec.range_start == "icru90-integrated-from-lowest-node-v1"
    assert spec.domain.species == frozenset({"proton", "alpha"})
    assert spec.domain.material_rule == stopping_models.MATERIAL_RULE_WATER_I78
    assert spec.domain.e_min_per_u == {"proton": 1.0, "alpha": 1.0}
    assert spec.domain.e_max_per_u == {"proton": 500.0, "alpha": 250.0}
    assert spec.options["points_per_decade"] == 200
    digest = spec_sha256(spec)
    assert digest not in {spec_sha256(s) for k, s in STOPPING_MODELS.items() if k != MID}
    assert [r.status for r in STOPPING_EVIDENCE[digest]] == ["consistency-only"]
    # the circular rows: the dataset and the rows named by plan row PER-2
    assert STOPPING_CIRCULAR_ROWS[digest] == (
        ICRU_ID,
        "I2 (ICRU 90 part)",
        "decision-0038 ICRU 90 comparisons",
    )


def test_capabilities_list_the_model_with_pin_licence_role_and_record() -> None:
    import ionmc

    cap = ionmc.capabilities()["physics"]["stopping_models"]["models"][MID]
    assert cap["data_pins"] == [
        {
            "dataset_id": ICRU_ID,
            "licence": "Geant4 Software License",
            "use_role": "construction-permitted",
        }
    ]
    assert [r["status"] for r in cap["evidence"]] == ["consistency-only"]
    assert cap["domain"]["species"] == ["alpha", "proton"]


@pytest.mark.parametrize("role", ["evaluation-only", "lineage-only"])
def test_role1_a_non_construction_pin_cannot_be_resolved(
    role: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = STOPPING_MODELS[MID]
    weak = replace(spec, data=(replace(spec.data[0], use_role=role),))
    monkeypatch.setitem(stopping_models._REGISTRY, MID, weak)
    with pytest.raises(ValueError, match="use_role"):
        resolve(MID, cache_dir="/nonexistent")  # refused before any cache access


# ---- start range and synthetic construction (cache-independent) ----


def test_start_range_is_exact_on_a_power_law_with_constant_s_only_below_the_lowest_node() -> None:
    c, p, a, e0 = 60.0, 0.75, 4.0, 0.2
    e = np.geomspace(e0, 100.0, 30)
    s = c * e**-p
    got = icru90_start_range(e, s, a)
    exact = a * e0 ** (p + 1.0) / c + a / (c * (p + 1.0)) * (1.0 - e0 ** (p + 1.0))
    assert got == pytest.approx(exact, rel=REL)
    # a node exactly at the start energy adds no extra node; a lowest node at the start gives
    # the constant-S term only
    e1 = np.append(e[e < 1.0], [1.0, 3.0])
    s1 = c * e1**-p
    assert icru90_start_range(e1, s1, a) == pytest.approx(exact, rel=REL)
    assert icru90_start_range(np.array([1.0, 2.0]), np.array([10.0, 8.0]), 4.0) == pytest.approx(
        0.4, rel=REL
    )
    with pytest.raises(ValueError, match="starts at"):
        icru90_start_range(np.array([2.0, 3.0]), np.array([1.0, 1.0]), 1.0)


def test_synthetic_model_tables_windows_and_identity() -> None:
    model = build_icru90_source(_synthetic())
    assert isinstance(model, Icru90WaterSource)
    assert model.window_per_u(PROTON) == (1.0, 500.0) and model.window_per_u(ALPHA) == (1.0, 250.0)
    for proj, c, p in ((PROTON, 50.0, 0.8), (ALPHA, 150.0, 0.7)):
        t = model.table(WATER, proj)
        assert t.energy_per_u[0] == 1.0
        np.testing.assert_allclose(t.s_el_mass, c * t.energy_per_u**-p, rtol=1e-12)
        assert t.metadata["model_id"] == MID and t.metadata["version"] == "1"
        assert t.metadata["decision"] == "0044" and t.metadata["range_construction"]
        assert table_model_id({"metadata": dict(t.metadata)}, proj) == MID
        # range: exact power-law integral from the lowest node plus the constant-S term
        e_lo = 0.1  # lowest node per u of the synthetic arrays
        a = float(proj.a)
        r_start = a * e_lo ** (p + 1.0) / c + a / (c * (p + 1.0)) * (1.0 - e_lo ** (p + 1.0))
        assert t.csda_range_g_cm2[0] == pytest.approx(r_start, rel=REL)
    assert model.content_sha256_for(PROTON) != model.content_sha256_for(ALPHA)
    assert len({model.content_sha256, model.content_sha256_for(PROTON)}) == 2
    # an unverified (synthetic) object is not the registered model: no pin, spec_of unregistered
    assert spec_of(model)[0] == "unregistered" and spec_of(model)[1].data == ()


def test_synthetic_dom1_and_dom4_raise_naming_model_species_material_and_domain() -> None:
    model = build_icru90_source(_synthetic())
    for proj in (DEUTERON, CARBON12, OXYGEN16):
        with pytest.raises(UnsupportedCombinationError) as ei:
            model.table(WATER, proj)
        for part in (MID, f"'{proj.name}'", "'water'", "energy any", "domain"):
            assert part in str(ei.value)
    for i_ev in (75.0, 78.5, 80.0):
        with pytest.raises(UnsupportedCombinationError) as ei:
            model.table(water(i_ev), PROTON)
        msg = str(ei.value)
        for part in (MID, "'proton'", "'water'", f"I = {i_ev} eV", "domain", "78.0 eV only"):
            assert part in msg
    model.table(water(78.0), PROTON)
    with pytest.raises(UnsupportedCombinationError, match="'pmma'.*domain"):
        model.table(PMMA, PROTON)  # structural rule (DOM-2), named by model and domain


def test_synthetic_source_must_cover_the_windows() -> None:
    d = _synthetic()
    short = replace(
        d, proton_energy_mev=d.proton_energy_mev[:-8], proton_stopping=d.proton_stopping[:-8]
    )
    # the proton window is capped at the array end (no extrapolation)
    m = build_icru90_source(short)
    assert m.window_per_u(PROTON)[1] < 500.0
    high = replace(
        d, alpha_energy_mev=d.alpha_energy_mev[20:], alpha_stopping=d.alpha_stopping[20:]
    )
    with pytest.raises(ValueError, match="outside the array span"):  # starts above 1 MeV/u
        build_icru90_source(high)


# ---- cache-dependent: the loaded model ----


def test_resolve_builds_the_registered_model_from_the_verified_cache(
    cdir: Path, data: Icru90Water
) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    assert isinstance(model, Icru90WaterSource)
    assert spec_of(model) == (MID, STOPPING_MODELS[MID])
    changed = Icru90WaterSource(replace(model.proton, points_per_decade=100), model.alpha)
    assert spec_of(changed)[0] == "unregistered"
    unverified = build_icru90_source(replace(data, sha256=None, dataset_id=None))
    assert spec_of(unverified)[0] == "unregistered"
    for proj in (PROTON, ALPHA):
        t = model.table(WATER, proj)
        meta = t.metadata
        assert meta["model_id"] == MID and meta["version"] == "1" and meta["decision"] == "0044"
        assert meta["dataset_id"] == ICRU_ID and meta["dataset_sha256"] == data.sha256
        assert meta["content_sha256"] == model.content_sha256_for(proj)
        assert meta["I_eV"] == 78.0 and meta["source"] == "tabulated"


def test_cc_hash_new_model_metadata_carry_the_id_and_the_hash_differs(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    nist_like = BetheStoppingSource()
    for proj in (PROTON, ALPHA):
        t = model.table(WATER, proj)
        b = BetheStoppingSource().table(WATER, proj)
        assert "model_id" in t.metadata and "model_id" not in b.metadata
        h_icru = TransportTables.from_stopping_tables([t], projectile=proj).sha256
        h_bethe = TransportTables.from_stopping_tables([b], projectile=proj).sha256
        assert h_icru != h_bethe
        # same arrays but another model id give another hash (the id is in the identity)
        other = replace(t, metadata={**t.metadata, "model_id": "icru90-water-v2"})
        assert TransportTables.from_stopping_tables([other], projectile=proj).sha256 != h_icru
    assert nist_like.table(WATER, PROTON).metadata.get("model_id") is None


def test_cc_span_loader_counts_and_coverage(cdir: Path, data: Icru90Water) -> None:
    assert data.proton_energy_mev.size == 57 and data.proton_stopping.size == 57
    assert data.alpha_energy_mev.size == 49 and data.alpha_stopping.size == 49
    model = resolve(MID, cache_dir=str(cdir))
    for proj, e_tot in ((PROTON, data.proton_energy_mev), (ALPHA, data.alpha_energy_mev)):
        lo, hi = model.window_per_u(proj)
        e_u = e_tot / proj.a
        assert e_u[0] <= lo and hi <= e_u[-1] * (1 + 1e-12)
    assert model.window_per_u(PROTON) == (1.0, 500.0)
    assert model.window_per_u(ALPHA) == (1.0, 250.0)


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=["proton", "alpha"])
def test_cc_node_fine_grid_equals_the_loglog_interpolant_of_the_array(
    cdir: Path, data: Icru90Water, gen: Any, proj: Any
) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    assert gen.cc_node_fine(model, data, proj) <= REL
    chord = gen.cc_node_chord(model, data, proj)  # report-only: aggregates, no assertion on size
    assert chord["n_points"] > 0 and math.isfinite(chord["max_abs_relative_deviation"])
    assert chord["rms_relative_deviation"] <= chord["max_abs_relative_deviation"]


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=["proton", "alpha"])
def test_cc_int_aggregates_per_decade_are_finite_and_cover_the_window(
    cdir: Path, data: Icru90Water, gen: Any, proj: Any
) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    out = gen.cc_int(data, proj, model.window_per_u(proj))
    assert [d["decade_per_u_mev"] for d in out] == [[1.0, 10.0], [10.0, 100.0], [100.0, 1000.0]]
    for d in out:
        assert d["n_points"] > 100 and set(d) == {
            "decade_per_u_mev",
            "window_per_u_mev",
            "n_points",
            "max_abs_relative_deviation",
            "rms_relative_deviation",
        }
        assert 0.0 < d["rms_relative_deviation"] <= d["max_abs_relative_deviation"] < 0.05


def test_cc_int_monotone_cubic_helper_reproduces_nodes_and_preserves_monotonicity(
    gen: Any,
) -> None:
    x = np.linspace(0.0, 5.0, 12)
    y = -0.8 * x + 0.05 * np.sin(2 * x)
    y = np.minimum.accumulate(y)  # monotone decreasing data
    d = gen.pchip_slopes(x, y)
    np.testing.assert_allclose(gen.pchip_eval(x, y, d, x), y, rtol=0, atol=1e-13)
    xq = np.linspace(0.0, 5.0, 400)
    assert np.all(np.diff(gen.pchip_eval(x, y, d, xq)) <= 1e-12)


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=["proton", "alpha"])
def test_cc_i2_node_arrays_equal_the_i2_construction(cdir_all: Path, gen: Any, proj: Any) -> None:
    model, data, pstar, astar = gen.load_inputs(cdir_all)
    star = pstar if proj is PROTON else astar
    r_m, r_i, e = gen.cc_i2_arrays(model, data, star, proj)
    assert e[0] == 1.0 and r_m[0] == 0.0 and r_i[0] == 0.0  # exact at the first node
    assert np.all(np.abs(r_m[1:] / r_i[1:] - 1.0) <= REL)
    assert gen.cc_i2_max_rel(model, data, star, proj) <= REL
    offset = gen.cc_i2_start_offset(data, star, proj)  # report-only: one number per projectile
    assert math.isfinite(offset)


def test_committed_aggregates_are_reproduced_and_hold_no_table_values(
    cdir_all: Path, gen: Any
) -> None:
    docs = gen.compute(cdir_all)
    text_all = ""
    for name, fresh in docs.items():
        path = RESULTS / f"{name}.json"
        committed = json.loads(path.read_text())
        _same(fresh, committed, name)
        text_all += path.read_text().replace(" ", "")
        assert path.read_text().endswith("}\n") and not path.read_text().endswith("\n\n")
    model, data, pstar, astar = gen.load_inputs(cdir_all)
    for arr in (
        data.proton_stopping,
        data.alpha_stopping,
        data.proton_energy_mev,
        data.alpha_energy_mev,
        pstar.s_electronic,
        pstar.csda_range,
        astar.s_electronic,
        astar.csda_range,
    ):
        for v in arr:
            s = f"{v:.8g}"
            assert len(s) < 7 or s not in text_all, "a table value appears in a committed JSON"
    cc = json.loads((RESULTS / "cc-int.json").read_text())
    assert cc["species"]["proton"]["cc_span"]["n_nodes"] == 57
    assert cc["species"]["alpha"]["cc_span"]["n_nodes"] == 49
    assert cc["species"]["proton"]["cc_node_max_rel_dev_at_fine_nodes"] <= REL
    assert cc["species"]["alpha"]["cc_i2_max_rel_dev_node_arrays"] <= REL


def _same(a: Any, b: Any, path: str = "") -> None:
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            _same(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _same(x, y, f"{path}[{i}]")
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), path
    else:
        assert a == b, path


# ---- domains (DOM-1, DOM-3, DOM-4, DOM-5) ----


def _refused(cfg: SimulationConfig, *parts: str) -> None:
    with pytest.raises(UnsupportedCombinationError) as ei:
        validate(cfg)
    msg = str(ei.value)
    for part in parts:
        assert part in msg, (part, msg)


@pytest.mark.parametrize("proj", [CARBON12, OXYGEN16], ids=lambda p: p.name)
def test_dom1_other_species_fail_closed(cdir: Path, proj: Any) -> None:
    _refused(_cfg(proj, 100.0), MID, proj.name, "'water'", "energy any", "domain")


def test_dom1_deuteron_fails_closed_at_the_table_and_at_validate(cdir: Path) -> None:
    # the deuteron is not a source projectile of the engine at all (earlier, pre-existing rule);
    # the model itself has no deuteron table
    with pytest.raises(UnsupportedCombinationError):
        validate(_cfg(DEUTERON, 100.0))
    with pytest.raises(UnsupportedCombinationError, match="deuteron.*domain"):
        resolve(MID, cache_dir=str(cdir)).table(WATER, DEUTERON)


def test_dom3_energy_above_the_window_fails_closed_for_both_species(cdir: Path) -> None:
    # The alpha window end (250 MeV/u) is also the pre-existing source rule of decision 0043
    # section 5, which `validate` applies first; the model's own window check (the interval
    # [E_cut, E0 + 6 sigma_E] against the per-species window) is called directly as well.
    for cfg in (_cfg(ALPHA, 251.0), _cfg(ALPHA, 245.0, sigma_mev=10.0)):  # 980 + 60 MeV
        with pytest.raises(UnsupportedCombinationError, match="250"):
            validate(cfg)
    validate(_cfg(ALPHA, 200.0, sigma_mev=1.0))
    sel = select(None, MID)
    check_window(sel, ALPHA, "water", "source interval", 200.0, 1000.0)  # the exact end passes
    with pytest.raises(UnsupportedCombinationError) as ei:
        check_window(sel, ALPHA, "water", "source interval", 200.0, 1004.0)  # E0 + 6 sigma_E
    for part in (MID, "'alpha'", "'water'", "251 MeV/u", "domain"):
        assert part in str(ei.value)
    _refused(_cfg(PROTON, 520.0), MID, "'proton'", "520", "domain")
    _refused(_cfg(PROTON, 480.0, sigma_mev=10.0), MID, "'proton'", "domain")  # 540 MeV
    validate(_cfg(PROTON, 480.0, sigma_mev=1.0))
    blk = validate(_cfg(ALPHA, 100.0)).stopping_identity()
    assert blk["domain_resolved"] == {"alpha": {"e_min_per_u": 1.0, "e_max_per_u": 250.0}}


def test_dom4_requested_water_i_other_than_78_ev_fails_closed(cdir: Path) -> None:
    for i_ev in (75.0, 77.9, 78.1, 80.0):
        box = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), water(i_ev))
        for proj in (PROTON, ALPHA):
            _refused(_cfg(proj, geometry=box), MID, proj.name, f"I = {i_ev} eV", "domain")
    ok = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), water(78.0))
    validate(_cfg(PROTON, geometry=ok))
    # a structurally non-water material of any I is refused by the structural rule (DOM-2)
    pmma = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), PMMA)
    _refused(_cfg(geometry=pmma), MID)
    # named otherwise but structurally water with I = 78 eV: accepted
    real = Material("tank liquid", 1.0, {"H": 0.111894, "O": 0.888106}, 78.0)
    validate(_cfg(geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), real)))


def test_dom5_nuclear_with_the_icru90_model_is_refused_naming_the_model(cdir: Path) -> None:
    _refused(_cfg(nuclear=True, nuclear_table_id="x"), MID, "deuteron", "nuclear=True")


# ---- DOM-6 ----


def test_dom6_missing_object_hash_mismatch_and_parse_failure_raise(
    cdir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        resolve(MID, cache_dir=str(empty))
    monkeypatch.setenv("IONMC_CACHE_DIR", str(empty))
    with pytest.raises(UnsupportedCombinationError, match="FileNotFoundError"):
        select(None, MID)  # the selection point fails closed with the loader's message
    monkeypatch.undo()
    # hash mismatch: a copy of the ICRU object with other bytes under the same name
    bad = tmp_path / "bad"
    (bad / "objects").mkdir(parents=True)
    (bad / "manifests").mkdir()
    sha = stopping_models.DATASETS[ICRU_ID].sha256
    shutil.copy(cache.manifest_path(ICRU_ID, cdir), cache.manifest_path(ICRU_ID, bad))
    obj = cache.object_path(sha, bad)
    shutil.copy(cache.object_path(sha, cdir), obj)
    os.chmod(obj, 0o644)  # cache objects may be read-only (host runner: hard-linked store)
    obj.write_bytes(obj.read_bytes() + b"x")
    with pytest.raises(cache.IntegrityError):
        resolve(MID, cache_dir=str(bad))

    # parse failure of an object that passed the hash check
    def broken(_path: Any, **_k: Any) -> None:
        raise ValueError("array T0_proton not found in source")

    monkeypatch.setattr(stopping_icru90, "load_icru90_water", broken)
    with pytest.raises(ValueError, match="T0_proton"):
        resolve(MID, cache_dir=str(cdir))


# ---- PER-2 / PER-3 (ICRU case) ----


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=lambda p: p.name)
def test_per3_data_used_lists_the_single_icru_pin_with_the_projectile_table(
    cdir: Path, proj: Any
) -> None:
    spec = STOPPING_MODELS[MID]
    pin = stopping_models._canonical(spec.data[0])
    eff = validate(_cfg(proj))
    blk = eff.summary()["stopping"]
    validate_identity_block(blk)
    assert blk["model_id_requested"] == blk["model_id_effective"] == MID
    assert blk["data_declared"] == [pin]
    assert blk["data_used"] == [
        {"pin": pin, "role_in_this_run": "construction", "table": [proj.name, "water"]}
    ]
    assert [r["status"] for r in blk["evidence_at_run"]] == ["consistency-only"]
    # PER-2: circular rows, and the MCS proton path of an ion run
    assert blk["circular_rows"] == [
        ICRU_ID,
        "I2 (ICRU 90 part)",
        "decision-0038 ICRU 90 comparisons",
    ]
    assert ("mcs_proton_path_stopping" in blk) == (proj.z >= 2)
    if proj.z >= 2:
        assert blk["mcs_proton_path_stopping"] == "bethe-ionmc-v1"
    model = resolve(MID, cache_dir=str(cdir))
    assert blk["content_sha256"] == model.content_sha256_for(proj)
    assert json.loads(json.dumps(blk)) == blk
    # negative control: a pin of a projectile that was not transported is rejected
    other = "alpha" if proj is PROTON else "proton"
    wrong = {**blk, "data_used": [{**blk["data_used"][0], "table": [other, "water"]}]}
    with pytest.raises(ValueError, match="not transported"):
        validate_identity_block(wrong)


def test_single_model_assertion_names_the_model_for_both_species(cdir: Path) -> None:
    for proj in (PROTON, ALPHA):
        eff = validate(_cfg(proj))
        assert eff.stopping_identity()["model_id_effective"] == MID
        t = eff.summary()["stopping"]
        assert t["kind"] == "tabulated" and t["version"] == "1" and t["decision"] == "0044"
