"""Model ``nist-pstar-astar-water-2005`` (decision 0044; V3-013 step S3).

Rows ID-NIST, DOM-1 (species), DOM-2, DOM-5, DOM-6, ROLE-1 and the NIST cases of PER-3 of
``validation/plans/v3-013-acceptance.md``. The tests that need the PSTAR/ASTAR objects find the
cache like the other cache-dependent tests (``cache.resolve_cache_dir()``: ``IONMC_CACHE_DIR`` or
the default) and skip with a message when an object is absent; ``IONMC_REQUIRE_DATA=1`` makes that
a failure. No table value is stored in Git: the fixture holds identities and hashes only.
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.data import cache
from ionmc.data.nist_star import load_star_table
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import PMMA, WATER, Material
from ionmc.physics import stopping_models
from ionmc.physics.projectiles import ALPHA, CARBON12, DEUTERON, PROTON
from ionmc.physics.stopping import NistStarStoppingSource
from ionmc.physics.stopping_evidence import STOPPING_CIRCULAR_ROWS, STOPPING_EVIDENCE
from ionmc.physics.stopping_models import STOPPING_MODELS, resolve, spec_of, spec_sha256
from ionmc.physics.stopping_nist import NistStarCombinedSource
from ionmc.physics.stopping_selection import select, validate_identity_block
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource
from ionmc.transport.tables import TransportTables, thaw

MID = "nist-pstar-astar-water-2005"
PSTAR_ID, ASTAR_ID = "nist-pstar-water-2005", "nist-astar-water-2005"
FIXTURE = Path(__file__).parent / "data" / "nist_identity_afe42b3.json"


@pytest.fixture(scope="module")
def cdir() -> Path:
    d = cache.resolve_cache_dir()
    try:
        for ds in (PSTAR_ID, ASTAR_ID):
            cache.verify(ds, d)
    except (FileNotFoundError, cache.IntegrityError) as exc:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(
                f"IONMC_REQUIRE_DATA=1 but the NIST PSTAR/ASTAR objects are not cached: {exc}"
            )
        pytest.skip(f"NIST PSTAR/ASTAR objects are not cached in {d}")
    return d


def _old_way(cdir: Path, program_id: str) -> NistStarStoppingSource:
    return NistStarStoppingSource(load_star_table(cache.verify(program_id, cdir)))


def _cfg(proj: Any = PROTON, e_u: float | None = None, **physics: Any) -> SimulationConfig:
    e_u = e_u if e_u is not None else (20.0 if proj.z == 1 else 50.0)
    geo = physics.pop("geometry", None) or BoxPhantom(
        (-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER
    )
    physics.setdefault("stopping_model", MID)
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e_u * proj.a, 0.0, 0.0),
        geometry=geo,
        scoring=(ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 2.0), (10, 10, 20)),),
        physics=PhysicsOptions(nuclear=physics.pop("nuclear", False), **physics),
        run=RunOptions(backend="python", precision="float64", seed=1, n_histories=2, n_batches=2),
    )


def _strip(x: Any, key: str = "retrieved_at") -> Any:
    if isinstance(x, dict):
        return {k: _strip(v, key) for k, v in x.items() if k != key}
    if isinstance(x, list):
        return [_strip(v, key) for v in x]
    return x


# ---- ID-NIST ----


@pytest.mark.parametrize(
    ("name", "proj", "ds"), [("proton", PROTON, PSTAR_ID), ("alpha", ALPHA, ASTAR_ID)]
)
def test_id_nist_tables_identity_and_hash_equal_the_single_program_source(
    cdir: Path, name: str, proj: Any, ds: str
) -> None:
    old = _old_way(cdir, ds)
    model = resolve(MID, cache_dir=str(cdir))
    assert isinstance(model, NistStarCombinedSource)
    t_old, t_new = old.table(WATER, proj), model.table(WATER, proj)
    assert dict(t_new.metadata) == dict(t_old.metadata)  # byte-identical per-table metadata
    assert "model_id" not in t_new.metadata  # the id goes into the identity block, not the table
    a = TransportTables.from_stopping_tables([t_old], projectile=proj)
    b = TransportTables.from_stopping_tables([t_new], projectile=proj)
    assert thaw(b.identity) == thaw(a.identity) and b.sha256 == a.sha256
    # against the values computed at the S2b commit before any S3 edit (afe42b3); retrieved_at is
    # the cache manifest's timestamp and ties the hash to this cache
    fx = json.loads(FIXTURE.read_text())
    assert fx["base_commit"] == "afe42b3"
    assert _strip(json.loads(json.dumps(thaw(b.identity), default=str))) == _strip(
        fx["identity"][name]
    )
    if fx["metadata"][name]["retrieved_at"] == t_new.metadata["retrieved_at"]:
        assert b.sha256 == fx["sha256"][name]


def test_id_nist_single_model_assertion_names_the_model_for_both_species(cdir: Path) -> None:
    for proj in (PROTON, ALPHA):
        eff = validate(_cfg(proj))
        assert eff.stopping_identity()["model_id_effective"] == MID


# ---- spec_of ----


def test_spec_of_bare_program_is_unregistered_and_the_combined_source_is_the_model(
    cdir: Path,
) -> None:
    assert spec_of(_old_way(cdir, PSTAR_ID))[0] == "unregistered"
    assert spec_of(_old_way(cdir, ASTAR_ID))[0] == "unregistered"
    model = resolve(MID, cache_dir=str(cdir))
    assert spec_of(model) == (MID, STOPPING_MODELS[MID])
    changed = NistStarCombinedSource(replace(model.pstar, points_per_decade=100), model.astar)
    assert spec_of(changed)[0] == "unregistered"
    capped = NistStarCombinedSource(replace(model.pstar, e_max_per_u=100.0), model.astar)
    assert spec_of(capped)[0] == "unregistered"


def test_registered_spec_is_cache_independent_and_carries_two_construction_pins() -> None:
    spec = STOPPING_MODELS[MID]
    assert spec.kind == "tabulated" and spec.version == "1" and spec.decision == "0044"
    assert [p.dataset_id for p in spec.data] == [PSTAR_ID, ASTAR_ID]
    from ionmc.data.registry import DATASETS

    for pin in spec.data:
        ds = DATASETS[pin.dataset_id]
        assert (pin.archive_sha256, pin.parser, pin.licence) == (ds.sha256, "nist_star", ds.license)
        assert pin.member is None and pin.use_role == "construction-permitted"
        assert ds.role == "evaluation"  # the dataset's own role is unchanged by the model
    assert spec.domain.species == frozenset({"proton", "alpha"})
    assert spec.domain.material_rule == stopping_models.MATERIAL_RULE_WATER
    assert spec.domain.e_min_per_u == {"proton": 1.0, "alpha": 1.0}
    digest = spec_sha256(spec)
    assert [r.status for r in STOPPING_EVIDENCE[digest]] == ["consistency-only"]
    assert STOPPING_CIRCULAR_ROWS[digest] == (PSTAR_ID, ASTAR_ID)


# ---- DOM-1 (species), DOM-2, DOM-5 ----


def test_dom1_other_species_fail_closed(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    for proj in (DEUTERON, CARBON12):
        with pytest.raises(UnsupportedCombinationError, match="no table for projectile"):
            model.table(WATER, proj)


def test_dom1_carbon_run_is_refused(cdir: Path) -> None:
    with pytest.raises(UnsupportedCombinationError, match="carbon12"):
        validate(_cfg(CARBON12, 100.0))


def test_dom2_non_water_material_fails_by_the_structural_rule_not_the_name(cdir: Path) -> None:
    pmma = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), PMMA)
    with pytest.raises(UnsupportedCombinationError):
        validate(_cfg(geometry=pmma))
    # named "water" but not water: refused; named otherwise but structurally water: accepted
    fake = Material("water", 1.2, {"H": 0.111894, "O": 0.888106}, 78.0)
    with pytest.raises(UnsupportedCombinationError):
        validate(_cfg(geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), fake)))
    real = Material("tank liquid", 1.0, {"H": 0.111894, "O": 0.888106}, 78.0)
    validate(_cfg(geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), real)))


def test_dom5_nuclear_with_the_nist_model_is_refused(cdir: Path) -> None:
    with pytest.raises(UnsupportedCombinationError, match="nist-star"):
        validate(_cfg(nuclear=True, nuclear_table_id="x"))


def test_window_is_the_resolved_table_window_not_the_spec_cap(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    lo_p, hi_p = model.window_per_u(PROTON)
    lo_a, hi_a = model.window_per_u(ALPHA)
    assert lo_p == lo_a == 1.0 and hi_p == 500.0 and hi_a < 500.0
    blk = validate(_cfg(ALPHA)).stopping_identity()
    assert blk["domain_resolved"] == {"alpha": {"e_min_per_u": lo_a, "e_max_per_u": hi_a}}
    with pytest.raises(UnsupportedCombinationError, match="MeV/u"):
        validate(_cfg(ALPHA, hi_a * 1.01))


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
    # hash mismatch: a copy of the cache whose ASTAR object has other bytes under the same name
    bad = tmp_path / "bad"
    shutil.copytree(cdir, bad, ignore=shutil.ignore_patterns("*icru*"))
    obj = bad / "objects" / stopping_models.DATASETS[ASTAR_ID].sha256
    os.chmod(obj, 0o644)  # cache objects may be read-only (host runner: hard-linked store)
    obj.write_bytes(obj.read_bytes() + b"x")
    with pytest.raises(cache.IntegrityError):
        resolve(MID, cache_dir=str(bad))
    # parse failure of an object that passed the hash check
    from ionmc.physics import stopping_nist

    def broken(_path: Any, **_k: Any) -> None:
        raise ValueError("no data rows found")

    monkeypatch.setattr(stopping_nist, "load_star_table", broken)
    with pytest.raises(ValueError, match="no data rows"):
        resolve(MID, cache_dir=str(cdir))
    # a table of the wrong program under the right pin
    real = load_star_table(cache.verify(PSTAR_ID, cdir))
    monkeypatch.setattr(stopping_nist, "load_star_table", lambda _p, **_k: real)
    with pytest.raises((ValueError, cache.IntegrityError)):
        resolve(MID, cache_dir=str(cdir))


# ---- ROLE-1 ----


@pytest.mark.parametrize("role", ["evaluation-only", "lineage-only"])
def test_role1_a_non_construction_pin_of_the_nist_model_cannot_be_resolved(
    role: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = STOPPING_MODELS[MID]
    weak = replace(spec, data=(replace(spec.data[0], use_role=role), spec.data[1]))
    monkeypatch.setitem(stopping_models._REGISTRY, MID, weak)
    with pytest.raises(ValueError, match="use_role"):
        resolve(MID, cache_dir="/nonexistent")  # refused before any cache access


# ---- PER-3 (NIST cases) ----


@pytest.mark.parametrize(
    ("proj", "used", "name"), [(PROTON, PSTAR_ID, "proton"), (ALPHA, ASTAR_ID, "alpha")]
)
def test_per3_data_used_lists_only_the_pin_read_for_the_transported_projectile(
    cdir: Path, proj: Any, used: str, name: str
) -> None:
    spec = STOPPING_MODELS[MID]
    pins = {p.dataset_id: stopping_models._canonical(p) for p in spec.data}
    blk = validate(_cfg(proj)).summary()["stopping"]
    validate_identity_block(blk)
    assert blk["data_declared"] == [pins[PSTAR_ID], pins[ASTAR_ID]]
    assert blk["data_used"] == [
        {
            "pin": pins[used],
            "role_in_this_run": "construction",
            "table": [name, "water"],
            "cache_dir": str(cdir),
        }
    ]
    assert blk["cache_dir"] == str(cdir)
    assert blk["circular_rows"] == [PSTAR_ID, ASTAR_ID]
    assert [r["status"] for r in blk["evidence_at_run"]] == ["consistency-only"]
    # content digest: of the parsed arrays of the program read, differs between the programs
    model = resolve(MID, cache_dir=str(cdir))
    assert blk["content_sha256"] == model.content_sha256_for(proj)
    assert model.content_sha256_for(PROTON) != model.content_sha256_for(ALPHA)
    # negative control: the other program's pin was not read in this run
    other = ASTAR_ID if used == PSTAR_ID else PSTAR_ID
    wrong = {
        **blk,
        "data_used": [
            {
                **blk["data_used"][0],
                "pin": pins[other],
                "table": ["alpha" if name == "proton" else "proton", "water"],
            }
        ],
    }
    with pytest.raises(ValueError, match="not transported"):
        validate_identity_block(wrong)


def test_per3_block_is_json_round_trippable_and_outside_the_table_hash(cdir: Path) -> None:
    eff = validate(_cfg(PROTON))
    blk = eff.summary()["stopping"]
    assert json.loads(json.dumps(blk)) == blk
    assert (
        eff.tables.sha256
        == TransportTables.from_stopping_tables(
            [_old_way(cdir, PSTAR_ID).table(WATER, PROTON)], projectile=PROTON
        ).sha256
    )


# ---- identity from the live arrays (review finding 1) ----


def _combined_with(
    src: NistStarCombinedSource, program: str, **arrays: Any
) -> NistStarCombinedSource:
    """A new combined source whose ``program`` table has other arrays (the metadata unchanged)."""
    base = {"PSTAR": src.pstar, "ASTAR": src.astar}
    st = base[program].star_table
    new = NistStarStoppingSource(replace(st, **arrays))
    if program == "PSTAR":
        return NistStarCombinedSource(new, src.astar)
    return NistStarCombinedSource(src.pstar, new)


def test_loaded_arrays_are_read_only(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    for src in (model.pstar, model.astar):
        for name in ("energy_mev", "s_electronic", "s_nuclear", "csda_range", "detour"):
            arr = getattr(src.star_table, name)
            assert not arr.flags.writeable
            with pytest.raises(ValueError, match="read-only"):
                arr[0] = 1.0


def test_registered_digests_match_the_loaded_content_and_the_spec(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    digests = STOPPING_MODELS[MID].options["content_sha256"]
    assert digests == {
        "proton": model.content_sha256_for(PROTON),
        "alpha": model.content_sha256_for(ALPHA),
    }
    assert spec_of(model)[0] == MID


@pytest.mark.parametrize("program", ["PSTAR", "ASTAR"])
@pytest.mark.parametrize("array", ["s_electronic", "csda_range", "energy_mev"])
def test_modified_content_with_the_pinned_metadata_is_unregistered(
    cdir: Path, program: str, array: str
) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    st = (model.pstar if program == "PSTAR" else model.astar).star_table
    changed = getattr(st, array).copy()  # a writable copy; the loaded array cannot be written
    changed[3] *= 1.001
    other = _combined_with(model, program, **{array: changed})
    assert other.pstar.star_table.dataset_id == PSTAR_ID  # the pinned metadata are intact
    mid, spec = spec_of(other)
    assert mid == "unregistered" and spec.decision == "user"
    modified = PSTAR_ID if program == "PSTAR" else ASTAR_ID
    # no registered pin is claimed for the modified program (the unmodified one keeps its own)
    assert modified not in {pin.dataset_id for pin in spec.data}
    assert spec_sha256(spec) not in STOPPING_EVIDENCE  # the registered evidence is not attached


def test_constructed_object_claiming_the_pinned_metadata_is_unregistered(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))

    def fake(st: Any) -> NistStarStoppingSource:
        return NistStarStoppingSource(
            replace(st, s_electronic=st.s_electronic * 2.0, csda_range=st.csda_range * 2.0)
        )

    other = NistStarCombinedSource(fake(model.pstar.star_table), fake(model.astar.star_table))
    assert other.pstar.star_table.sha256 == model.pstar.star_table.sha256
    mid, spec = spec_of(other)
    assert mid == "unregistered" and spec_sha256(spec) not in STOPPING_EVIDENCE
    assert spec_of(model)[0] == MID and spec_sha256(spec_of(model)[1]) in STOPPING_EVIDENCE


def test_subclasses_never_keep_the_registered_nist_identity(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))
    assert model.is_exact_implementation() and spec_of(model)[0] == MID

    class SubComposite(NistStarCombinedSource):
        def table(self, material: Any, projectile: Any) -> Any:
            return super().table(material, projectile)

    class SubStar(NistStarStoppingSource):
        pass

    assert spec_of(SubComposite(model.pstar, model.astar))[0] == "unregistered"
    star = model.pstar
    inner = SubStar(star.star_table, star.e_min_per_u, star.e_max_per_u, star.points_per_decade)
    mid, spec = spec_of(NistStarCombinedSource(inner, model.astar))
    assert mid == "unregistered" and spec_sha256(spec) not in STOPPING_EVIDENCE
    assert spec_of(model)[0] == MID


def test_subclass_overriding_table_and_the_exactness_method_is_unregistered(cdir: Path) -> None:
    model = resolve(MID, cache_dir=str(cdir))

    class Liar(NistStarCombinedSource):
        def is_exact_implementation(self) -> bool:
            return True

        def instance_spec(self) -> Any:
            return STOPPING_MODELS[MID]

        def table(self, material: Any, projectile: Any) -> Any:
            return super().table(material, projectile)

    mid, spec = spec_of(Liar(model.pstar, model.astar))
    assert mid == "unregistered" and spec.data == () and spec.decision == "user"
    assert spec_sha256(spec) not in STOPPING_EVIDENCE


# ---- provider type and consumed-content binding (REVIEW-bcd2acbc) ----


class _StatefulStarTable:
    """A look-alike of a StarTable: pinned arrays everywhere except inside ``table()``."""

    def __init__(self, real: Any) -> None:
        object.__setattr__(self, "_real", real)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._real, name)

    def _arr(self, name: str, factor: float) -> Any:
        import sys

        a = getattr(self._real, name)
        if sys._getframe(2).f_code.co_name == "table":
            b = a * factor
            b.setflags(write=False)
            return b
        return a

    @property
    def energy_mev(self) -> Any:
        return self._arr("energy_mev", 1.0)

    @property
    def s_electronic(self) -> Any:
        return self._arr("s_electronic", 1.02)

    @property
    def csda_range(self) -> Any:
        return self._arr("csda_range", 1.0)


def _lookalike_model(cdir: Path) -> NistStarCombinedSource:
    model = resolve(MID, cache_dir=str(cdir))

    def wrap(src: NistStarStoppingSource) -> NistStarStoppingSource:
        return NistStarStoppingSource(
            _StatefulStarTable(src.star_table),  # type: ignore[arg-type]
            src.e_min_per_u,
            src.e_max_per_u,
            src.points_per_decade,
        )

    return NistStarCombinedSource(wrap(model.pstar), wrap(model.astar), cache_dir=model.cache_dir)


def test_stateful_provider_is_unregistered_and_the_consumed_digest_fails_validate(
    cdir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = _lookalike_model(cdir)
    assert spec_of(fake)[0] == "unregistered"  # the provider is not an exact StarTable
    # bypass only the provider type check: the per-table consumed digest still fails the run
    monkeypatch.setattr(stopping_models, "_nist_provider_exact", lambda t: True)
    assert spec_of(fake)[0] == MID  # pinned on every read outside table()
    cfg = SimulationConfig(
        source=_cfg().source,
        geometry=_cfg().geometry,
        scoring=_cfg().scoring,
        physics=PhysicsOptions(nuclear=False, stopping=fake),
        run=_cfg().run,
    )
    assert cfg.physics.selection.observed[0] == MID
    with pytest.raises(UnsupportedCombinationError, match="identity changed"):
        validate(cfg)


def test_untouched_nist_model_identity_is_stable_and_tables_carry_the_consumed_digest(
    cdir: Path,
) -> None:
    cfg = _cfg()
    assert cfg.physics.selection.observed[0] == MID
    validate(cfg)
    model = cfg.physics.stopping_source
    t = model.table(WATER, PROTON)
    assert t.consumed_content_sha256 == STOPPING_MODELS[MID].options["content_sha256"]["proton"]
    assert "consumed_content_sha256" not in t.metadata


# ---- exact primitive settings and the rebuilt-table closure (REVIEW-dd6310ee) ----


def _in_frame(*names: str) -> bool:
    import sys

    return sys._getframe(2).f_code.co_name in names


class _CondBool:
    """Truthy on inspection, falsy inside the physics kernel."""

    def __bool__(self) -> bool:
        return not _in_frame("bethe_mass_stopping")


class _CondInt:
    """An int-like setting that is 200 everywhere except inside ``log_grid`` (100 there)."""

    def __int__(self) -> int:
        return 200

    def __index__(self) -> int:
        return 200

    def __float__(self) -> float:
        return 200.0

    def __eq__(self, other: object) -> bool:
        return other == 200

    __hash__ = None  # type: ignore[assignment]

    def __mul__(self, other: float) -> float:
        return other * (100.0 if _in_frame("log_grid") else 200.0)

    __rmul__ = __mul__


def test_conditional_bethe_setting_is_unregistered_and_fails_on_the_rebuilt_tables(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from ionmc.physics.stopping import BetheOptions, BetheStoppingSource

    src = BetheStoppingSource(options=BetheOptions(shell=_CondBool()))  # type: ignore[arg-type]
    assert spec_of(src)[0] == "unregistered"  # exact-type rule: shell is not a bool
    cfg = _cfg(stopping=src, stopping_model=None)
    assert cfg.physics.selection.requested_id == "unregistered"
    monkeypatch.setattr(stopping_models, "_bethe_exact", lambda s: True)
    cfg = _cfg(stopping=src, stopping_model=None)
    assert cfg.physics.selection.requested_id == "bethe-ionmc-v1"
    with pytest.raises(UnsupportedCombinationError, match="identity changed"):
        validate(cfg)


def test_stateful_grid_setting_is_unregistered_and_fails_on_the_rebuilt_tables(
    cdir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model = resolve(MID, cache_dir=str(cdir))

    def wrap(src: NistStarStoppingSource) -> NistStarStoppingSource:
        return NistStarStoppingSource(
            src.star_table,
            src.e_min_per_u,
            src.e_max_per_u,
            _CondInt(),  # type: ignore[arg-type]
        )

    fake = NistStarCombinedSource(wrap(model.pstar), wrap(model.astar), cache_dir=model.cache_dir)
    assert spec_of(fake)[0] == "unregistered"
    monkeypatch.setattr(stopping_models, "_exact_types", lambda *a: True)
    assert spec_of(fake)[0] == MID
    cfg = _cfg(stopping=fake, stopping_model=None)
    with pytest.raises(UnsupportedCombinationError, match="identity changed"):
        validate(cfg)


@pytest.mark.parametrize("model_id", ["bethe-ionmc-v1", MID])
def test_untouched_registered_models_validate_with_identical_rebuilt_tables(
    cdir: Path, model_id: str
) -> None:
    cfg = _cfg(stopping_model=model_id)
    assert cfg.physics.selection.observed[0] == model_id
    validate(cfg)
    cfg_obj = _cfg(stopping=resolve(model_id, cache_dir=str(cdir)), stopping_model=None)
    assert cfg_obj.physics.selection.observed[0] == model_id
    validate(cfg_obj)
