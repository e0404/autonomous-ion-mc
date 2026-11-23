"""Row OUT-12 (plan v3-006, amendment A-2): tier-2 ``verify_effective``.

Needs the data cache of this environment (CI and local tiers); every test rebuilds the configuration
of a saved bundle and compares ``EffectiveConfig.summary()`` with the persisted ``effective``.
"""

from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_v3_006_results_format import (  # noqa: E402
    _doc,
    _ion_config,
    _mutate,
    _nuclear_table_id,
    _voxel_result,
    lookup_bundle,  # noqa: F401  (fixture)
    make_result,
)

from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, PROTON
from ionmc.results import VerificationError, load, save, verify_effective
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

LOOKUP_FILE = Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"


def _tree_digest(b: Path) -> dict[str, str]:
    return {
        str(p.relative_to(b)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(b.rglob("*"))
        if p.is_file()
    }


@pytest.fixture(scope="module")
def alpha_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    cfg = SimulationConfig(
        PencilBeamSource(ALPHA, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 400.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 4.0), (4, 4, 10)),),
        PhysicsOptions(nuclear=False, stopping_model="icru90-water-v1", multiple_scattering=False),
        RunOptions("warp-cpu", "float64", 17, 16, n_batches=4),
    )
    return save(Simulation(cfg).run(), tmp_path_factory.mktemp("al") / "b")


def test_synthetic_default_and_voxel_bundles_are_reproduced(tmp_path: Path) -> None:
    for i, r in enumerate((make_result(tally=True), _voxel_result(30.0))):
        b = save(r, tmp_path / f"b{i}")
        before = _tree_digest(b)
        rep = verify_effective(b)
        assert rep.ok and rep.status == "reproduced" and rep.differences == ()
        assert any("evidence_at_run" in x for x in rep.exclusions)
        assert _tree_digest(b) == before  # never modifies the bundle


def test_carbon_and_icru90_alpha_are_reproduced(tmp_path: Path, alpha_bundle: Path) -> None:
    b = save(Simulation(_ion_config()).run(), tmp_path / "c")
    assert verify_effective(b).ok
    assert verify_effective(load(alpha_bundle)).ok  # a LoadedResult is accepted too


def test_nuclear_bundle_is_reproduced(tmp_path: Path) -> None:
    cfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 150.0, 0.0, 0.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 200.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 10.0), (4, 4, 20)),),
        PhysicsOptions(nuclear=True, nuclear_table_id=_nuclear_table_id(),
                       stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 16, 40, n_batches=4),
    )  # fmt: skip
    b = save(Simulation(cfg).run(), tmp_path / "n")
    rep = verify_effective(b)
    assert rep.ok, rep.differences


def test_lookup_with_input_files_has_no_exclusion_and_placeholders_report_them(
    lookup_bundle: Path,  # noqa: F811
) -> None:
    role = _doc(lookup_bundle)["input_manifest"]["files"][-1]["role"]
    assert role.startswith("lookup:")
    rep = verify_effective(lookup_bundle, input_files={role: LOOKUP_FILE})
    assert rep.ok and rep.lookups_from == {role: "input_file"}
    assert not any("bound_per_history" in x for x in rep.exclusions)
    ph = verify_effective(lookup_bundle)
    assert ph.ok and ph.lookups_from == {role: "placeholder"}
    assert any("bound_per_history" in x for x in ph.exclusions)


def test_mismatching_input_file_and_unknown_role_are_tier2_failures(
    lookup_bundle: Path,  # noqa: F811
    tmp_path: Path,
) -> None:
    role = _doc(lookup_bundle)["input_manifest"]["files"][-1]["role"]
    bad = tmp_path / "t.json"
    bad.write_bytes(LOOKUP_FILE.read_bytes().replace(b"\n", b" \n", 1))  # same table, other bytes
    with pytest.raises(VerificationError, match="sha256|digests"):
        verify_effective(lookup_bundle, input_files={role: bad})
    with pytest.raises(VerificationError, match="without a lookup"):
        verify_effective(lookup_bundle, input_files={"lookup:nope": LOOKUP_FILE})
    with pytest.raises(VerificationError, match="exceeds"):
        verify_effective(lookup_bundle, input_files={role: LOOKUP_FILE}, max_input_file_bytes=10)


def test_mutated_effective_pointer_is_a_tier2_failure_not_a_tier1_one(tmp_path: Path) -> None:
    b = save(make_result(tally=True), tmp_path / "b")

    def bump(d: dict[str, Any]) -> None:  # coherently: the N-channel bound is max_steps 2 pieces
        d["effective"]["max_steps"] += 1
        for c in d["effective"]["tallies"]["channels"]:
            if c["kind"] == "N":
                c["bound_per_history"] = (
                    float(d["effective"]["max_steps"]) * 2.0 * d["effective"]["scoring_pieces"]
                )

    _mutate(b, bump)
    # the persisted files are rewritten by _mutate: result.json is not hashed in `files`
    assert load(b).binding_limits  # tier 1 accepts (a declared limit)
    rep = verify_effective(b)
    assert not rep.ok and rep.status == "differs" and "/effective/max_steps" in rep.differences


def test_missing_cache_object_is_an_explicit_tier2_failure(
    alpha_bundle: Path, tmp_path: Path, monkeypatch: Any
) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    b = tmp_path / "copy"
    shutil.copytree(alpha_bundle, b)
    monkeypatch.setenv("IONMC_CACHE_DIR", str(empty))
    assert load(b).effective_binding == "bundle-local"  # tier 1 never needs the cache
    with pytest.raises(VerificationError, match="not re-derivable"):
        verify_effective(b)
    with pytest.raises(VerificationError, match="not re-derivable"):
        verify_effective(b, cache_dir=empty)


def test_input_file_size_must_equal_the_manifest_record() -> None:
    from ionmc.results import _input_lookup

    entry = {"sha256": {"status": "unavailable"}, "content_sha256": "0" * 64,
             "item_sha256": "0" * 64, "bytes": LOOKUP_FILE.stat().st_size + 1}  # fmt: skip
    with pytest.raises(VerificationError, match="manifest records"):
        _input_lookup("lookup:x", LOOKUP_FILE, entry, 2**30)


def test_user_stopping_object_is_reported_not_rederivable(tmp_path: Path) -> None:
    import dataclasses

    from ionmc.physics.stopping import BetheStoppingSource

    cfg = dataclasses.replace(
        _ion_config(), physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource())
    )
    rep = verify_effective(save(Simulation(cfg).run(), tmp_path / "b"))
    assert rep.status == "not-rederivable" and not rep.ok
    assert rep.not_rederivable and "user stopping object" in rep.not_rederivable[0]
