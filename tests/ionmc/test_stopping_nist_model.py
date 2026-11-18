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
        {"pin": pins[used], "role_in_this_run": "construction", "table": [name, "water"]}
    ]
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
