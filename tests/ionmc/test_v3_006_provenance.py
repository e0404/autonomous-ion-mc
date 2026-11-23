"""Rows PROV-1 .. PROV-6 and the PROV-7 expectation helper of plan v3-006 (decision 0045 sections
2-3): package tree oid, code identity (stamp, git, unknown), environment/data/models blocks, the
privacy rule of the path set P and the build stamp of a wheel (PROV-5, local validation)."""

from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import zipfile
from pathlib import Path
from typing import Any

import pytest

from ionmc import provenance as prov
from ionmc import results
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.results import ResultFormatError, load, save
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "src" / "ionmc"
HEX40 = "0123456789abcdef" * 2 + "01234567"
GIT = shutil.which("git")


def _git(args: list[str], cwd: Path, env_extra: dict[str, str] | None = None) -> str:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update({"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_SYSTEM": os.devnull})
    env.update(env_extra or {})
    cmd = ["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", *args]
    return subprocess.run(
        cmd, cwd=cwd, env=env, check=True, capture_output=True, text=True
    ).stdout.strip()


def _copy_pkg(dest: Path) -> Path:
    shutil.copytree(PKG, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    return dest


def _config(model: str) -> SimulationConfig:
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0),
        BoxPhantom((-3.0, -3.0, 0.0), (6.0, 6.0, 6.0), WATER),
        (ScoringGrid((-1.2, -1.2, 0.1), (0.3, 0.3, 0.7), (8, 8, 6)),),
        PhysicsOptions(nuclear=False, stopping_model=model),
        RunOptions("python", "float64", 5, 12, n_batches=3),
    )


def _require_datasets(*dataset_ids: str) -> None:
    """Skip a cache-backed case without the dataset (CI tier); fail when the dataset is absent
    under IONMC_REQUIRE_DATA=1 or IONMC_LV_TIER=1. A corrupt object, malformed manifest or pin
    mismatch (``cache.IntegrityError``) is never a skip."""
    from ionmc.data import cache

    cdir = cache.resolve_cache_dir(None)
    required = "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER"))
    for dataset_id in dataset_ids:
        try:
            cache.verify(dataset_id, cdir)
        except FileNotFoundError as exc:
            if required:
                pytest.fail(f"{dataset_id} is required here but not cached: {exc}")
            pytest.skip(f"{dataset_id} is not in the cache (IONMC_CACHE_DIR)")


_MODEL_DATASETS = {
    "nist-pstar-astar-water-2005": ("nist-pstar-water-2005",),
    "icru90-water-v1": ("geant4-icru90-stopping-11.4.2",),
}


def _bundle(tmp: Path, model: str, name: str = "b") -> Path:
    _require_datasets(*_MODEL_DATASETS.get(model, ()))
    return save(Simulation(_config(model)).run(), tmp / name)  # failures propagate


def _doc(bundle: Path) -> dict[str, Any]:
    return json.loads((bundle / "result.json").read_text(encoding="utf-8"))  # type: ignore[no-any-return]


# ------------------------------------------------------------------------------------- PROV-1


def test_prov1_tree_oid_equals_git(tmp_path: Path) -> None:
    if GIT is None:
        _git_required()
    # the committed package tree is reproduced byte for byte (git archive of HEAD) ...
    arch = tmp_path / "arch"
    arch.mkdir()
    archive = subprocess.run(
        ["git", "archive", "HEAD", "src/ionmc"], cwd=ROOT, check=True, capture_output=True
    ).stdout
    subprocess.run(["tar", "-x", "-C", str(arch)], input=archive, check=True)
    head_tree = _git(["rev-parse", "HEAD:src/ionmc"], ROOT)
    assert prov.package_tree_oid(arch / "src" / "ionmc") == head_tree
    # ... and on a clean checkout the live directory equals it too
    if not _git(["status", "--porcelain", "--", "src/ionmc"], ROOT):
        assert prov.package_tree_oid() == head_tree
    # build stamp, caches and bytecode never enter the oid
    clean = prov.package_tree_oid(arch / "src" / "ionmc")
    (arch / "src" / "ionmc" / "_build_info.json").write_text("{}", encoding="utf-8")
    (arch / "src" / "ionmc" / "__pycache__").mkdir()
    (arch / "src" / "ionmc" / "__pycache__" / "x.pyc").write_bytes(b"x")
    assert prov.package_tree_oid(arch / "src" / "ionmc") == clean


def test_prov1_git_sorting_and_modes(tmp_path: Path) -> None:
    """Directories sort as ``name/``; the executable bit selects mode 100755."""
    if GIT is None:
        _git_required()
    pkg = tmp_path / "repo" / "src" / "ionmc"
    (pkg / "a").mkdir(parents=True)
    (pkg / "a" / "z.py").write_text("z\n")
    (pkg / "a.py").write_text("a\n")  # "a.py" < "a/" in git's order, ".py" sorts before "/"
    (pkg / "a-b").write_text("ab\n")  # "a-b" < "a.py" < "a/"
    (pkg / "run.sh").write_text("#!/bin/sh\n")
    (pkg / "run.sh").chmod(0o755)
    (pkg / "empty").mkdir()
    _git(["init", "-q"], tmp_path / "repo")
    _git(["add", "-A"], tmp_path / "repo")
    _git(["commit", "-qm", "x"], tmp_path / "repo")
    assert prov.package_tree_oid(pkg) == _git(["rev-parse", "HEAD:src/ionmc"], tmp_path / "repo")


# ------------------------------------------------------------------------------------- PROV-2


def test_prov2_changed_byte_stamp_and_unknown(tmp_path: Path) -> None:
    pkg = _copy_pkg(tmp_path / "pkg" / "ionmc")
    oid = prov.package_tree_oid(pkg)
    # no stamp, no git (the copy lies in no work tree): nulls, "unknown", a reason
    code = prov.code_identity(pkg)
    if code["source"] == "git":
        pytest.skip("the temporary directory lies inside a Git work tree of src/ionmc")
    assert code["source"] == "unavailable" and code["status"] == "unknown"
    assert code["git_commit"] is None and code["git_dirty"] is None and code["reason"]
    assert code["package_tree_oid"] == oid
    # stamp present: dirty iff the runtime oid differs from the stamped one
    stamp = {"git_commit": "a" * 40, "package_tree_oid": oid, "status": "ok"}
    (pkg / prov.BUILD_STAMP_NAME).write_text(json.dumps(stamp), encoding="utf-8")
    clean = prov.code_identity(pkg)
    assert (clean["source"], clean["git_commit"], clean["git_dirty"]) == ("stamp", "a" * 40, False)
    assert clean["package_tree_oid"] == oid and clean["status"] == "ok"
    with open(pkg / "results.py", "ab") as fh:
        fh.write(b"#")  # one changed byte
    dirty = prov.code_identity(pkg)
    assert dirty["package_tree_oid"] != oid and dirty["git_dirty"] is True
    # a stamp without a commit (build from an sdist): unknown, never an error
    stamp = {"git_commit": None, "package_tree_oid": oid, "status": "unknown"}
    (pkg / prov.BUILD_STAMP_NAME).write_text(json.dumps(stamp), encoding="utf-8")
    unk = prov.code_identity(pkg)
    assert unk["status"] == "unknown" and unk["git_commit"] is None and unk["source"] == "stamp"


def test_prov2_git_fallback_ignores_foreign_git_vars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    if GIT is None:
        _git_required()
    repo, foreign = tmp_path / "a", tmp_path / "foreign"
    pkg = _copy_pkg(repo / "src" / "ionmc")
    for r in (repo, foreign):
        r.mkdir(exist_ok=True)
        _git(["init", "-q"], r)
    (foreign / "f.txt").write_text("f\n")
    _git(["add", "-A"], foreign)
    _git(["commit", "-qm", "f"], foreign)
    _git(["add", "-A"], repo)
    _git(["commit", "-qm", "a"], repo)
    commit = _git(["rev-parse", "HEAD"], repo)
    monkeypatch.setenv("GIT_DIR", str(foreign / ".git"))
    monkeypatch.setenv("GIT_WORK_TREE", str(foreign))
    monkeypatch.setenv("GIT_INDEX_FILE", str(foreign / "nowhere"))
    code = prov.code_identity(pkg)
    assert (code["source"], code["git_commit"], code["git_dirty"]) == ("git", commit, False)
    assert code["package_tree_oid"] == _git(["rev-parse", "HEAD:src/ionmc"], repo)
    with open(pkg / "results.py", "ab") as fh:
        fh.write(b"#")
    assert prov.code_identity(pkg)["git_dirty"] is True
    # a package that is not src/ionmc of the work tree it lies in is not that repository's package
    odd = _copy_pkg(repo / "vendor" / "ionmc")
    assert prov.code_identity(odd)["source"] == "unavailable"


# ------------------------------------------------------------------------------------- PROV-3/4


def test_prov3_default_model_and_registry_append(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    b = _bundle(tmp_path, "bethe-ionmc-v1")
    doc = _doc(b)
    st, m = doc["effective"]["stopping"], doc["provenance"]["models"]
    assert (m["stopping_model_id_effective"], m["stopping_version"]) == (
        st["model_id_effective"], st["version"],
    )  # fmt: skip
    assert m["stopping_spec_sha256"] == st["spec_sha256"]
    assert m["evidence_at_run"] == st["evidence_at_run"] and m["evidence_at_run"]
    assert doc["run"]["cache_dir_source"] == "not-used"
    # a registry append after the save changes nothing in what load returns
    from ionmc.physics import stopping_evidence as ev

    extra = dict(ev.STOPPING_EVIDENCE)
    for k, v in list(extra.items()):
        extra[k] = (*v, v[0])
    monkeypatch.setattr(ev, "STOPPING_EVIDENCE", extra)
    again = load(b).doc["provenance"]["models"]
    assert again["evidence_at_run"] == m["evidence_at_run"]


@pytest.mark.parametrize("model", ["nist-pstar-astar-water-2005", "icru90-water-v1"])
def test_prov3_prov4_tabulated_models(tmp_path: Path, model: str) -> None:
    b = _bundle(tmp_path, model)
    doc = _doc(b)
    st, p = doc["effective"]["stopping"], doc["provenance"]
    assert p["models"]["stopping_model_id_effective"] == st["model_id_effective"]
    assert p["models"]["stopping_spec_sha256"] == st["spec_sha256"]
    assert p["models"]["evidence_at_run"] == st["evidence_at_run"]
    if model == "icru90-water-v1":
        statuses = [r["status"] for r in st["evidence_at_run"]]
        assert len(statuses) == 2 and statuses[0] == "consistency-only"
        assert statuses[1] == "topas-agreement-shared-lineage:V-C1:PASS"
    pins = [x["pin"] for x in st["data_used"]]
    assert pins and [x["pin"] for x in p["data"]["stopping_data_used"]] == pins
    assert all("cache_dir" not in json.dumps(x) for x in p["data"]["stopping_data_used"])
    want = "IONMC_CACHE_DIR" if os.environ.get("IONMC_CACHE_DIR") else "default"
    assert doc["run"]["cache_dir_source"] in (want, "explicit")
    load(b)
    # PROV-4: a removed pin is refused by load
    doc["provenance"]["data"]["stopping_data_used"].pop()
    (b / "result.json").write_bytes(results.canonical_bytes(doc))
    with pytest.raises(ResultFormatError):
        load(b)


def test_prov4_data_block_nuclear_and_lookups() -> None:
    summary = {
        "stopping": {"data_used": [{"pin": {"id": "d"}, "cache_dir": "/c", "x": [1]}]},
        "tables": {"sha256": "s", "stopping_source": "t", "materials": []},
        "nuclear": {
            "table_id": "tab",
            "cache_dir": "/c",
            "pins": [{"cache_dir": "/c", "s": "q"}],
        },
    }
    lk = [{"name": "f", "content_sha256": "c"}]
    d = prov.data_block(summary, lk)
    assert d["nuclear"] == {"table_id": "tab", "pins": [{"s": "q"}]}
    assert d["stopping_data_used"] == [{"pin": {"id": "d"}, "x": [1]}] and d["lookups"] == lk
    assert prov.data_block(summary, [], {"table_id": "o"})["nuclear"] == {"table_id": "o"}
    assert prov.data_block({**summary, "nuclear": None}, [])["nuclear"] is None


# ------------------------------------------------------------------------------------- PROV-6


def test_prov6_environment_block_has_no_identity() -> None:
    env = prov.environment_block()
    text = json.dumps(env)
    assert set(env) == {
        "python_version", "platform", "numpy_version", "warp_version", "warp_devices", "warp_error",
    }  # fmt: skip
    for secret in (os.uname().nodename, str(Path.home())):
        assert secret not in text


def test_prov6_privacy_rule() -> None:
    names, home = ["root", "mybox"], "/home/ann"
    doc: dict[str, Any] = {
        "run": {"output_path": "/home/ann/out/root", "cache_dir": "/home/ann/.cache/ionmc"},
        "input_manifest": {"files": [{"path": "/home/ann/root/f.toml"}]},
        "effective": {"stopping": {"cache_dir": "/home/ann/.cache"}, "mcs_mapping": {"e": "root"}},
        "provenance": {"environment": {"platform": "Linux-6.1-x86_64", "n": "channel/groot"}},
    }
    assert prov.privacy_violations(doc, names=names, home=home) == []  # P holds paths; "root" ok
    bad = json.loads(json.dumps(doc))
    bad["effective"]["x"] = "/home/ann/y"  # home path outside P
    bad["effective"]["y"] = "a/root/b"  # a path component equal to the user name
    bad["effective"]["z"] = "C:\\Users\\mybox\\x"  # backslash separator
    bad["provenance"]["environment"]["host"] = "mybox"  # bare identity string in an identity block
    found = prov.privacy_violations(bad, names=names, home=home)
    assert {f.split(":")[0] for f in found} == {
        "/effective/x", "/effective/y", "/effective/z", "/provenance/environment/host",
    }  # fmt: skip


def test_prov6_real_bundle_and_save_guard(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    b = _bundle(tmp_path, "bethe-ionmc-v1")
    doc = _doc(b)
    assert prov.privacy_violations(doc) == []
    report = (b / "report.md").read_text(encoding="utf-8")
    assert str(tmp_path) not in report and doc["run"]["output_path"] not in report
    assert doc["provenance"]["code"]["package_tree_oid"] and doc["run"]["cache_dir_source"]
    result = Simulation(_config("bethe-ionmc-v1")).run()
    result.environment["provenance_snapshot"]["environment"]["platform"] = (
        f"/srv/{socket.gethostname()}/x"
    )
    with pytest.raises(ResultFormatError):
        save(result, tmp_path / "leak")
    assert not (tmp_path / "leak").exists()


# ------------------------------------------------------------------------------------- PROV-7


def test_prov7_expectation_helper() -> None:
    sha = "b" * 40
    exp = prov.host_runner_identity_expectation(sha, "c" * 40)
    assert exp == {
        "source": "git", "status": "ok", "git_dirty": False, "git_commit": sha,
        "package_tree_oid": "c" * 40,
    }  # fmt: skip
    code = {**exp, "ionmc_version": "x", "reason": ""}
    assert prov.host_runner_identity_mismatches(code, exp) == []
    assert len(prov.host_runner_identity_mismatches({**code, "git_dirty": True}, exp)) == 1
    assert "package_tree_oid" not in prov.host_runner_identity_expectation(sha)


# ------------------------------------------------------------------------------------- PROV-5


# stderr texts of an offline build environment that cannot obtain hatchling (never a hook error)
_ENV_BUILD_ERRORS = (
    "failed to resolve", "no solution found", "network", "offline", "failed to fetch",
    "failed to download", "dns error", "connection",
)  # fmt: skip


def _lv_tier() -> bool:
    """The declared local-validation selector: every PROV-5 failure fails there."""
    return os.environ.get("IONMC_LV_TIER") == "1"


def _git_required() -> None:
    """Git is a required tool of the local-validation tier (IONMC_LV_TIER=1); CI may skip."""
    if _lv_tier():
        pytest.fail("IONMC_LV_TIER=1 but git is not available")
    pytest.skip("git is not available")


def _uv_build(src: Path, out: Path) -> Path:
    uv = shutil.which("uv")
    if uv is None:
        if _lv_tier():
            pytest.fail("IONMC_LV_TIER=1 but uv is not on PATH")
        pytest.skip("uv is not available")
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    cache = ROOT / ".uv-cache"
    if cache.is_dir():
        env.setdefault("UV_CACHE_DIR", str(cache))
    proc = subprocess.run(
        [uv, "build", "--wheel", "--offline", "--out-dir", str(out), str(src)],
        cwd=src, env=env, capture_output=True, text=True, timeout=300, check=False,
    )  # fmt: skip
    if proc.returncode != 0:
        msg = f"uv build --wheel --offline failed: {proc.stderr.strip()[-600:]}"
        environmental = any(k in proc.stderr.lower() for k in _ENV_BUILD_ERRORS)
        if environmental and not _lv_tier():
            pytest.skip(msg)
        pytest.fail(msg)
    (wheel,) = out.glob("ionmc-*.whl")
    return wheel


def test_prov5_wheel_build_stamp(tmp_path: Path) -> None:
    if GIT is None:
        _git_required()
    before = _git(["status", "--porcelain"], ROOT)
    wheel = _uv_build(ROOT, tmp_path / "out")
    assert _git(["status", "--porcelain"], ROOT) == before  # the checkout is not touched
    assert not (PKG / prov.BUILD_STAMP_NAME).exists()
    unpacked = tmp_path / "unpacked"
    with zipfile.ZipFile(wheel) as zf:
        zf.extractall(unpacked)
    stamp = json.loads((unpacked / "ionmc" / prov.BUILD_STAMP_NAME).read_text(encoding="utf-8"))
    assert set(stamp) == {"git_commit", "package_tree_oid", "status"}
    assert stamp["git_commit"] == _git(["rev-parse", "HEAD"], ROOT)
    clean = prov.package_tree_oid(PKG) == _git(["rev-parse", "HEAD:src/ionmc"], ROOT)
    assert stamp["status"] == ("ok" if clean else "dirty")
    assert stamp["package_tree_oid"] == prov.package_tree_oid(unpacked / "ionmc")
    code = prov.code_identity(unpacked / "ionmc")
    assert (code["source"], code["git_dirty"], code["git_commit"]) == (
        "stamp",
        not clean,
        stamp["git_commit"],
    )


def test_prov5_build_without_repository_stamps_unknown(tmp_path: Path) -> None:
    src = tmp_path / "src_copy"
    src.mkdir()
    for name in ("pyproject.toml", "hatch_build.py", "README.md"):
        shutil.copy2(ROOT / name, src / name)
    _copy_pkg(src / "src" / "ionmc")
    wheel = _uv_build(src, tmp_path / "out")
    with zipfile.ZipFile(wheel) as zf:
        stamp = json.loads(zf.read(f"ionmc/{prov.BUILD_STAMP_NAME}"))
    assert stamp["git_commit"] is None and stamp["status"] == "unknown"
    assert not (src / "src" / "ionmc" / prov.BUILD_STAMP_NAME).exists()


def test_prov5_pyproject_and_lock_unchanged() -> None:
    if GIT is None:
        _git_required()
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'requires = ["hatchling==1.32.4"]' in text
    assert '[tool.hatch.build.hooks.custom]\npath = "hatch_build.py"' in text
    assert _git(["diff", "--stat", "HEAD", "--", "uv.lock"], ROOT) == ""
    assert "src/ionmc/_build_info.json" in (ROOT / ".gitignore").read_text(encoding="utf-8")


def test_prov5_dirty_build_is_stamped_dirty(tmp_path: Path) -> None:
    if GIT is None:
        _git_required()
    repo = tmp_path / "repo"
    repo.mkdir()
    for name in ("pyproject.toml", "hatch_build.py", "README.md", ".gitignore"):
        shutil.copy2(ROOT / name, repo / name)
    _copy_pkg(repo / "src" / "ionmc")
    _git(["init", "-q"], repo)
    _git(["add", "-A"], repo)
    _git(["commit", "-qm", "x"], repo)
    commit = _git(["rev-parse", "HEAD"], repo)
    with open(repo / "src" / "ionmc" / "results.py", "ab") as fh:
        fh.write(b"#")  # the build tree differs from the commit's package tree
    wheel = _uv_build(repo, tmp_path / "out")
    unpacked = tmp_path / "unpacked"
    with zipfile.ZipFile(wheel) as zf:
        zf.extractall(unpacked)
    stamp = json.loads((unpacked / "ionmc" / prov.BUILD_STAMP_NAME).read_text(encoding="utf-8"))
    assert stamp["status"] == "dirty" and stamp["git_commit"] == commit
    assert stamp["package_tree_oid"] == prov.package_tree_oid(unpacked / "ionmc")
    assert stamp["package_tree_oid"] != _git(["rev-parse", "HEAD:src/ionmc"], repo)
    code = prov.code_identity(unpacked / "ionmc")
    assert (code["source"], code["git_commit"], code["git_dirty"]) == ("stamp", commit, True)


def test_prov_run_time_identity_survives_save(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = Simulation(_config("bethe-ionmc-v1")).run()
    snap = result.environment["provenance_snapshot"]
    run_oid = snap["code"]["package_tree_oid"]
    monkeypatch.setattr(prov, "package_tree_oid", lambda *a, **k: "f" * 40)
    code = _doc(save(result, tmp_path / "b"))["provenance"]["code"]
    assert code["package_tree_oid"] == run_oid and code["saved_with"] == "f" * 40
    monkeypatch.undo()
    same = _doc(save(result, tmp_path / "c"))["provenance"]
    assert "saved_with" not in same["code"] and same["environment"] == snap["environment"]
    load(tmp_path / "b")


def test_prov_run_time_effective_snapshot_survives_registry_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from ionmc.physics import stopping_evidence as ev

    result = Simulation(_config("bethe-ionmc-v1")).run()
    run_summary = result.environment["provenance_snapshot"]["effective_summary"]
    run_evidence = json.loads(json.dumps(run_summary["stopping"]["evidence_at_run"]))
    extra = dict(ev.STOPPING_EVIDENCE)
    for k, v in list(extra.items()):
        extra[k] = (*v, v[0])
    monkeypatch.setattr(ev, "STOPPING_EVIDENCE", extra)
    live = result.effective_config.summary()["stopping"]["evidence_at_run"]
    assert len(live) == len(run_evidence) + 1  # the live registry differs from the run's state
    doc = _doc(save(result, tmp_path / "b"))
    assert doc["effective"]["stopping"]["evidence_at_run"] == run_evidence
    assert doc["provenance"]["models"]["evidence_at_run"] == run_evidence
    assert doc["effective"] == json.loads(results.canonical_bytes(run_summary))
    load(tmp_path / "b")
