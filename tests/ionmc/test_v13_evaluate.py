# ruff: noqa: E501
"""Tests of the sealed V3-013 evaluator on synthetic inputs (temporary git repositories).

Every check (a)-(e) is shown to pass on consistent data and to fail closed on a tampered input.
Nothing here reads a real run, TOPAS output or host record.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
STOP = ROOT / "validation/scripts/stopping"
for p in (str(SRC), str(STOP)):
    if p not in sys.path:
        sys.path.insert(0, p)
import v13_evaluate as ev  # noqa: E402


def _gen():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location(
        "_t_gen_v13_manifest", STOP / "gen_v13_manifest.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GEN = _gen()
AT = ev.rio.attest_module()
KEY_CW, KEY_HOST, KEY_REF = b"\x01" * 32, b"\x02" * 32, b"\x03" * 32
HZ, NZ, R_PLAN = 0.05, 200, 22.5
J_MEV = 1.602176634e-13
ENV = {"schema": "s", "digests": {"numpy": "a" * 64}, "versions": {"numpy": "2"},
       "file_counts": {"numpy": 1}, "sys_version": "3.12", "platform": "p"}  # fmt: skip
T_SEAL_END, T_SAMPLES_START, T_SAMPLES_END = (
    "2026-11-20T10:00:00+00:00", "2026-11-20T10:40:00+00:00", "2026-11-20T11:00:00+00:00",
)  # fmt: skip
T_ION, T_TOPAS = "2026-11-20T10:30:00+00:00", "2026-11-20T12:00:00+00:00"
ENERGY = {1: 75, 2: 125, 3: 200}
SHIFT = {1: 0.006, 2: 0.10, 3: 0.20}  # synthetic dR_pred [mm]
MODEL_ICRU, MODEL_DEF = ev.MODEL_ICRU90, ev.MODEL_DEFAULT


def git(repo: Path, *args: str) -> str:
    env = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@x", "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@x", "PATH": os.environ["PATH"], "HOME": str(repo)}  # fmt: skip
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True,
                          check=True, env=env).stdout.strip()  # fmt: skip


def commit(repo: Path, msg: str) -> str:
    git(repo, "add", "-A")
    git(repo, "-c", "commit.gpgsign=false", "commit", "-q", "-m", msg)
    return git(repo, "rev-parse", "HEAD")


def put(repo: Path, rel: str, data: bytes | str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data.encode() if isinstance(data, str) else data)


def bragg(zc: np.ndarray, shift: float = 0.0) -> np.ndarray:
    z = zc - shift
    tail = 0.5 * np.array([math.erfc((v - 7.0) / (math.sqrt(2.0) * 0.15)) for v in z])
    return (1.0 + 0.05 * z) * tail + 2.0 * np.exp(-((z - 6.9) ** 2) / (2 * 0.25**2)) * tail


ZC = (np.arange(NZ) + 0.5) * HZ


def sign(payload: dict[str, Any], seed: bytes, signer: str) -> dict[str, Any]:
    return dict(AT.sign_record(payload, AT.keypair_from_seed(seed), signer))


def review(
    rid: str, head: str, started: str, finished: str, *, seed: bytes | None = KEY_CW, **over: Any
) -> dict:
    report = {"head_sha": head, "verdict": "pass", "summary": "synthetic", "findings": []}
    payload = {
        "schema_version": 1, "review_id": rid, "task_id": "V3-013", "head_sha": head, "base_sha": "0" * 40,
        "status": "passed", "started_at": started, "finished_at": finished, "model": "m", "provider": "p",
        "tool": "t", "experiment_id": ev.review_records_module().current_experiment_id(),
        "report_sha256": hashlib.sha256(json.dumps(report).encode()).hexdigest(), "report": report,
    }  # fmt: skip
    payload.update(over)
    return payload if seed is None else sign(payload, seed, "codex-worker")


def topas_stdout(n: int) -> str:
    return (
        "Registered Physics Processes:\nmsc, eIoni, ionIoni\n\n"
        "Use ICRU90 data                                    1\n\n"
        "msc:  for e-  SubType= 10\n  GoudsmitSaunderson : Emin=    0 eV  Emax=  100 MeV Nbins=120\n\n"
        "msc:  for alpha  SubType= 10\n  UrbanMsc : Emin=    0 eV  Emax=  100 TeV Nbins=120\n\n"
        "msc:  for GenericIon  SubType= 10\n  UrbanMsc : Emin=    0 eV  Emax=  100 TeV Nbins=120\n\n"
        f"Begin processing of event 0\nParticle source Beam: Total number of histories: {n}\n"
    )


def idd_csv(energy_total: np.ndarray) -> str:
    dose = energy_total * J_MEV / (120.0 * 120.0 * HZ * 1e-6)
    head = (
        "# TOPAS Version: 4.3\n# Parameter File: input.txt\n# Results for scorer: IDD\n"
        "# Scored in component: Water\n# X in 1 bin  of 12 cm\n# Y in 1 bin  of 12 cm\n"
        f"# Z in {dose.size} bins of {HZ / 10} cm\n# DoseToMedium ( Gy ) : Sum   Standard_Deviation\n"
    )
    return head + "".join(f"0, 0, {i}, {float(v)!r}, 0\n" for i, v in enumerate(dose))


def predictions(tree: str, steps_sha: str) -> dict[str, Any]:
    r_def = 100.0
    energies = []
    for c in (1, 2, 3):
        d_d = SHIFT[c]
        energies.append({
            "case": c, "e_mev_per_u": ENERGY[c], "r_plan_mm": R_PLAN, "delta_r_pred_mm": SHIFT[c],
            "delta_r_pred_d_mm": d_d, "r_csda_default_mm": r_def,
            "vc2_precondition_holds": c != 1, "vc2d_informative": c != 1,
            "vc2d_gate": "substitute: CC-T",
            "vc2d_round_trip_tolerance_mm": {MODEL_DEF: 1e-4 * r_def, MODEL_ICRU: 1e-4 * (r_def + d_d)},
            "vc2d_difference_tolerance_mm": 1e-4 * r_def,
        })  # fmt: skip
    models, blocks = {}, {}
    for m in (MODEL_ICRU, MODEL_DEF):
        models[m] = {"model_id_effective": m, "version": 1, "spec_sha256": "s" + m, "content_sha256": "c" + m,
                     "tables_sha256": {str(c): f"t{m}{c}" for c in (1, 2, 3)},
                     "tables_materials": {str(c): [{"I_eV_effective": 78.0}] for c in (1, 2, 3)}}  # fmt: skip
        blocks[m] = {
            str(c): {
                "model_id_effective": m,
                "version": 1,
                "spec_sha256": "s" + m,
                "content_sha256": "c" + m,
                "evidence": [c],
            }
            for c in (1, 2, 3)
        }
    return {"schema": ev.PREDICTIONS_SCHEMA, "energies": energies, "models": models, "identity_blocks": blocks,
            "src_ionmc_tree": tree, "script": {"path": ev.STEPS_REL, "sha256": steps_sha}}  # fmt: skip


class World:
    """A synthetic repository, run records and the helpers to change them."""

    def __init__(self, tmp: Path, **o: Any) -> None:
        self.o = o
        self.repo = tmp / "repo"
        self.repo.mkdir()
        self.records = tmp / "records"
        self.out = tmp / "out"
        r = self.repo
        git(r, "init", "-q", "-b", "main")
        put(r, "src/ionmc/__init__.py", "x = 1\n")
        git(r, "-c", "commit.gpgsign=false", "commit", "-q", "--allow-empty", "-m", "empty")
        commit(r, "src")
        tree = git(r, "rev-parse", "HEAD:src/ionmc")
        for rel in (ev.OWN_REL, ev.LAUNCHER_REL, ev.ENV_DIGEST_REL, "pyproject.toml", "uv.lock"):
            put(r, rel, f"# {rel}\n")
        put(r, ev.STEPS_REL, "# steps\n")
        steps_sha = hashlib.sha256((r / ev.STEPS_REL).read_bytes()).hexdigest()
        kp = {
            s: AT.keypair_from_seed(sd)
            for s, sd in (
                ("codex-worker", KEY_CW),
                ("host-runner", KEY_HOST),
                ("reference-runner", KEY_REF),
            )
        }
        put(
            r,
            ev.TRUST_PATH,
            json.dumps(
                {s: [{"key_id": k.key_id, "public_key": k.public_key.hex()}] for s, k in kp.items()}
            ),
        )
        put(r, ev.PREDICTIONS, json.dumps(predictions(tree, steps_sha), indent=1))
        put(r, ev.ENV_SEAL, json.dumps(ENV))
        put(r, ev.VC3, json.dumps({"row": "V-C3"}))
        put(r, ev.PLAN_PATH, "# plan\n")
        for c in (1, 2, 3):
            for part in (*range(8), 19):
                name = f"he-water-{ENERGY[c]}mevu-v13-idd-em-part{part}" + (
                    "-smoke" if part == 19 else ""
                )
                src = ROOT / ev.BUNDLE_DIR / name
                for f in ("case.json", "input.txt"):
                    put(r, f"{ev.BUNDLE_DIR}/{name}/{f}", (src / f).read_bytes())
        self.manifest_doc = GEN.build(r)
        put(r, ev.MANIFEST, GEN.text(self.manifest_doc))
        self.seal = commit(r, "seal")
        self.manifest_sha = hashlib.sha256((r / ev.MANIFEST).read_bytes()).hexdigest()
        put(r, "samples-marker.txt", "records\n")
        self.samples = commit(r, "samples")
        self.pred = json.loads((r / ev.PREDICTIONS).read_text())
        self.n_run = 0
        self.ion_index: list[dict] = []
        self.top_index: list[dict] = []
        self.tomb: dict[str, Any] = {}
        self.build_ionmc(self.seal)
        self.build_topas(self.samples, parts=range(4))
        self.finish()

    # -- review pins and the final evaluated commit --
    def finish(
        self,
        *,
        seal_note_digest: str | None = None,
        drop_seal_review: bool = False,
        bad_signature: bool = False,
    ) -> None:
        r = self.repo
        rv = review(
            "REVIEW-" + "a" * 32, self.seal, "2026-11-20T09:00:00+00:00", T_SEAL_END, seed=KEY_CW
        )
        if bad_signature:
            rv["finished_at"] = "2026-11-19T00:00:00+00:00"  # changed after signing
        if drop_seal_review:
            (r / ev.SEAL_REVIEW).unlink(missing_ok=True)
        else:
            put(r, ev.SEAL_REVIEW, json.dumps(rv))
        put(
            r,
            ev.SAMPLES_REVIEW,
            json.dumps(review("REVIEW-" + "b" * 32, self.samples, T_SAMPLES_START, T_SAMPLES_END)),
        )
        digest = seal_note_digest or self.manifest_sha
        put(r, ev.PLAN_PATH, f"# plan\n\n- seal review REVIEW-{'a' * 32} at {self.seal}; manifest sha256 {digest}\n"
            f"- samples review REVIEW-{'b' * 32} at {self.samples}\n")  # fmt: skip
        self.head = commit(r, "pins")

    # -- ionmc runs --
    def curves(self, case: int, model: str) -> np.ndarray:
        shift = SHIFT[case] if model == MODEL_ICRU else 0.0
        rng = np.random.default_rng(case)
        noise = 1.0 + 1e-3 * rng.standard_normal((ev.N_BATCHES, 1))  # shared by both models (CRN)
        return noise * (bragg(ZC, shift) * HZ)[None, :]

    def build_ionmc(self, sha: str) -> None:
        o = self.o
        for row, case, model, backend in ev.expected_ionmc_runs():
            self.n_run += 1
            rid = f"RUN-20261120T10{self.n_run:04d}Z-{self.n_run:08x}"
            ddir = self.records / "published/validation/generated" / rid
            ddir.mkdir(parents=True)
            doc: dict[str, Any] = {**ev.step_expectation(row, case, backend), "model": model,
                                   "counters": {"a": 0, "b": 0}, "environment_digests": ENV,
                                   "warp_modules": None if backend == "python" else {"m/k.so": "1" * 64}}  # fmt: skip
            blk = self.pred["identity_blocks"][model][str(case)]
            pins = self.pred["models"][model]
            doc["effective_config"] = {"stopping": dict(blk), "tables": {"sha256": pins["tables_sha256"][str(case)],
                                                                          "materials": pins["tables_materials"][str(case)]}}  # fmt: skip
            doc["stopping_identity"] = doc["effective_config"]["stopping"]
            if row == "vc-sample":
                arr = self.curves(case, model)[:, None, None, :]
                npz = ddir / "sample.npz"
                np.savez_compressed(
                    npz, idd_batch_energy_mev=arr, idd_edges2_mm=np.arange(NZ + 1) * HZ
                )
                doc |= {
                    "npz": "sample.npz",
                    "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                }
            else:
                r_e0 = 100.0 + (SHIFT[case] if model == MODEL_ICRU else 0.0)
                doc |= {"r_e0_mm": r_e0, "z_m_mm": r_e0 + 1e-6}
            if (row, case, model, backend) == o.get("bad_identity"):
                doc["effective_config"]["stopping"]["spec_sha256"] = "forged"
            if (row, case, model, backend) == o.get("bad_warp"):
                doc["warp_modules"] = {"m/k.so": "2" * 64}
            (ddir / "step.stdout.txt").write_text(
                "#JSON-BEGIN\n" + json.dumps(doc) + "\n#JSON-END\n"
            )
            digest, _ = ev.rio.published_manifest(self.records / "published", rid)
            rec = {"run_id": rid, "sha": sha, "task_id": "V3-013", "succeeded": True, "exit_code": 0,
                   "timed_out": False, "dirty_before": False, "dirty_after": False,
                   "committed_state_unchanged": True, "output_publication_error": None,
                   "started_at": o.get("ion_start", T_ION), "finished_at": "2026-11-20T10:50:00+00:00",
                   "stdout_sha256": "0" * 64, "output_manifest_sha256": digest}  # fmt: skip
            signed = sign(
                rec,
                KEY_HOST if (row, case, model, backend) != o.get("bad_host_sig") else KEY_REF,
                "host-runner",
            )
            (self.records / "host-records").mkdir(exist_ok=True)
            (self.records / "host-records" / f"{rid}.json").write_text(json.dumps(signed))
            self.ion_index.append({"run_id": rid, "row": row, "case": case, "model": model, "backend": backend,
                                   "k": 0, "step_stdout": "step.stdout.txt"})  # fmt: skip

    # -- TOPAS runs --
    def build_topas(self, sha: str, parts: Any, case_list: Any = (1, 2, 3)) -> None:
        for case in case_list:
            for part in parts:
                self.topas_run(sha, case, part)

    def topas_run(
        self,
        sha: str,
        case: int,
        part: int,
        *,
        started: str | None = None,
        tamper_input: bool = False,
        stderr: str = "",
    ) -> None:
        name = f"he-water-{ENERGY[case]}mevu-v13-idd-em-part{part}"
        bdir = self.repo / ev.BUNDLE_DIR / name
        files = {}
        for f in ("case.json", "input.txt"):
            b = (bdir / f).read_bytes()
            files[f] = {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}
        if tamper_input:
            files["input.txt"] = {**files["input.txt"], "sha256": "f" * 64}
        case_doc = json.loads((bdir / "case.json").read_text())
        request = {"schema_version": 2, "code_sha": sha, "input_files": files, "case": case_doc,
                   "dirty": False, "engine": {"engine": "topas"}}  # fmt: skip
        cfg = ev.rio.canonical_digest(request)
        rid = f"REF-{cfg[:20]}-{hashlib.sha256(name.encode()).hexdigest()[:8]}"
        d = self.records / "topas" / rid
        (d / "work").mkdir(parents=True)
        rng = np.random.default_rng(100 * case + part)
        energy = (
            bragg(ZC, SHIFT[case])
            * HZ
            * ev.TOPAS_HISTORIES
            * (1.0 + 2e-3 * rng.standard_normal(NZ))
        )
        (d / "work/idd_dose.csv").write_text(idd_csv(energy))
        (d / "request.json").write_text(json.dumps(request))
        (d / "stdout.txt").write_text(topas_stdout(ev.TOPAS_HISTORIES))
        (d / "stderr.txt").write_text(stderr)
        listing = [{"path": p.relative_to(d).as_posix(), "bytes": p.stat().st_size,
                    "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(d.rglob("*")) if p.is_file()]  # fmt: skip
        by = {e["path"]: e["sha256"] for e in listing}
        work = {
            e["path"][5:]: {"sha256": e["sha256"], "bytes": e["bytes"]}
            for e in listing
            if e["path"].startswith("work/")
        }
        rec = sign({"schema_version": 3, "run_id": rid, "configuration_id": cfg, "code_sha": sha, "task_id": "V3-013",
                    "engine": "topas", "request_sha256": by["request.json"], "started_at": started or self.o.get("topas_start", T_TOPAS),
                    "finished_at": "2026-11-20T13:00:00+00:00", "seconds": 60.0, "exit_code": 0, "timed_out": False,
                    "missing_outputs": [], "succeeded": True, "output_count": len(work),
                    "output_files_sha256": ev.rio.canonical_digest(work), "stdout_sha256": by["stdout.txt"],
                    "stderr_sha256": by["stderr.txt"]}, KEY_REF, "reference-runner")  # fmt: skip
        (d / "record.json").write_text(json.dumps(rec))
        rb = (d / "record.json").read_bytes()
        listing.append(
            {"path": "record.json", "sha256": hashlib.sha256(rb).hexdigest(), "bytes": len(rb)}
        )
        (d / "artifacts.json").write_text(json.dumps({"run_id": rid, "sha": sha, "files": listing}))
        self.top_index.append({"bundle": name, "run_id": rid})

    def drop_topas(self, case: int, part: int) -> None:
        name = f"he-water-{ENERGY[case]}mevu-v13-idd-em-part{part}"
        for e in [e for e in self.top_index if e["bundle"] == name]:
            shutil.rmtree(self.records / "topas" / e["run_id"])
            self.top_index.remove(e)

    def write_index(self) -> None:
        (self.records / "index.json").write_text(
            json.dumps({"ionmc": self.ion_index, "topas": self.top_index})
        )

    def run(self, **kw: Any) -> dict[str, Any]:
        self.write_index()
        return ev.evaluate(self.repo, self.seal, self.head, self.records, self.out, launch_check=False,
                           environment=lambda: ENV, **kw)  # fmt: skip


@pytest.fixture()
def world(tmp_path: Path) -> World:
    return World(tmp_path)


# ------------------------------------------------------------------ consistent data passes
def test_consistent_synthetic_evidence_passes_every_check(world: World) -> None:
    res = world.run()
    assert all(r["valid"] for r in res["runs"]["ionmc"]), [
        r["problems"] for r in res["runs"]["ionmc"] if r["problems"]
    ]
    assert all(t["valid"] for t in res["runs"]["topas"])
    assert [v["status"] for v in res["vc1"].values()] == ["PASS"] * 3
    assert res["evidence_record"] == ev.EVIDENCE_PASS
    assert (
        res["vc2"]["1"]["gating"] is False
        and res["vc2"]["2"]["verdict"] == "PASS"
        and res["vc2"]["3"]["verdict"] == "PASS"
    )
    assert abs(res["vc2"]["2"]["delta_r80_mm"] - SHIFT[2]) < 0.01
    assert (
        res["vc2d"]["1/python"]["gating"] is False
        and res["vc2d"]["2/warp-cpu"]["verdict"] == "PASS"
    )
    assert res["integration"]["science_gates"] == "PASS" and res["integration"][
        "substitute_gate_energies"
    ] == ["1"]
    assert (world.out / ev.RESULT_NAME).is_file()
    assert set(res["vc1a"]) == {"1", "2", "3"} and res["vc1a"]["1"]["gating"] is False


# ------------------------------------------------------------------ (a) review records and pins
def test_missing_seal_review_record_stops_the_evaluation(world: World) -> None:
    world.finish(drop_seal_review=True)
    with pytest.raises(ev.RecordError, match="not a review record"):
        world.run()


def test_bad_signature_of_the_seal_review_stops_the_evaluation(world: World) -> None:
    world.finish(bad_signature=True)
    with pytest.raises(ev.RecordError, match="review record"):
        world.run()


def test_review_of_another_head_sha_stops_the_evaluation(world: World) -> None:
    rv = review("REVIEW-" + "a" * 32, world.samples, "2026-11-20T09:00:00+00:00", T_SEAL_END)
    put(world.repo, ev.SEAL_REVIEW, json.dumps(rv))
    world.head = commit(world.repo, "wrong head")
    with pytest.raises(ev.RecordError, match="head_sha"):
        world.run()


def test_plan_pin_with_a_wrong_manifest_digest_stops_the_evaluation(world: World) -> None:
    world.finish(seal_note_digest="0" * 64)
    with pytest.raises(ev.RecordError, match="manifest sha256"):
        world.run()


def test_pin_notes_must_occur_exactly_once() -> None:
    with pytest.raises(ev.RecordError, match="seal review pin note"):
        ev.parse_pins("nothing\n")
    note = f"seal review REVIEW-{'a' * 32} at {'1' * 40}; manifest sha256 {'2' * 64}\n"
    with pytest.raises(ev.RecordError, match="samples review pin note"):
        ev.parse_pins(note)
    with pytest.raises(ev.RecordError, match="found 2"):
        ev.parse_pins(note * 2)


# ------------------------------------------------------------------ (b) the sealed closure
def test_a_changed_sealed_file_breaks_the_manifest_check(world: World) -> None:
    put(world.repo, "uv.lock", "tampered\n")
    world.finish()
    with pytest.raises(ev.RecordError, match="uv.lock: sha256 differs"):
        world.run()


def test_a_changed_manifest_breaks_the_seal_blob_check(world: World) -> None:
    doc = json.loads((world.repo / ev.MANIFEST).read_text())
    doc["note"] = "edited after the seal"
    put(world.repo, ev.MANIFEST, json.dumps(doc))
    world.finish()
    with pytest.raises(ev.RecordError, match="differs between the seal"):
        world.run()


def test_a_manifest_that_lists_itself_or_misses_a_required_path_is_refused(tmp_path: Path) -> None:
    for mutate, msg in ((lambda f: f.__setitem__(ev.MANIFEST, "0" * 64), "list itself"),
                        (lambda f: f.pop(ev.OWN_REL), "required sealed paths")):  # fmt: skip
        sub = tmp_path / msg.split()[0]
        sub.mkdir()
        w = World(sub)
        doc = json.loads((w.repo / ev.MANIFEST).read_text())
        mutate(doc["files"])
        put(w.repo, ev.MANIFEST, json.dumps(doc))
        w.seal = commit(w.repo, "new seal")
        w.manifest_sha = hashlib.sha256((w.repo / ev.MANIFEST).read_bytes()).hexdigest()
        w.samples = w.seal
        w.finish()
        with pytest.raises(ev.RecordError, match=msg):
            w.run()


def test_a_changed_src_ionmc_tree_breaks_the_pin(world: World) -> None:
    doc = json.loads((world.repo / ev.MANIFEST).read_text())
    doc["src_ionmc_tree"] = "0" * 40
    put(world.repo, ev.MANIFEST, json.dumps(doc))
    world.seal = commit(world.repo, "bad tree pin")
    world.manifest_sha = hashlib.sha256((world.repo / ev.MANIFEST).read_bytes()).hexdigest()
    world.samples = world.seal
    world.finish()
    with pytest.raises(ev.RecordError, match="src/ionmc tree at the seal"):
        world.run()


# ------------------------------------------------------------------ (c) launch and environment
def test_a_different_environment_stops_the_evaluation(world: World) -> None:
    world.write_index()
    other = {**ENV, "digests": {"numpy": "b" * 64}}
    with pytest.raises(ev.RecordError, match="environment differs"):
        ev.evaluate(
            world.repo,
            world.seal,
            world.head,
            world.records,
            world.out,
            launch_check=False,
            environment=lambda: other,
        )


def test_the_evaluator_refuses_a_direct_launch(tmp_path: Path) -> None:
    env = {**os.environ, "PYTHONPATH": str(SRC)}
    p = subprocess.run([sys.executable, "-B", str(STOP / "v13_evaluate.py"), "--seal-sha", "0" * 40, "--code-sha",
                        "0" * 40, "--records", str(tmp_path), "--out", str(tmp_path / "o")],
                       capture_output=True, text=True, env=env, check=False)  # fmt: skip
    assert p.returncode == 2 and "not launched through sealed_launch.py" in p.stderr
    assert not (tmp_path / "o" / ev.RESULT_NAME).exists()


def test_the_sealed_launcher_allows_exactly_the_new_targets() -> None:
    text = (STOP / "sealed_launch.py").read_text()
    assert (
        '"validation/scripts/stopping/v13_evaluate.py"' in text
        and '"validation/scripts/stopping/env_digest.py"' in text
    )


def test_evaluation_commit_must_be_head(world: World) -> None:
    world.write_index()
    with pytest.raises(ev.RecordError, match="not the evaluated HEAD"):
        ev.evaluate(
            world.repo,
            world.seal,
            world.seal,
            world.records,
            world.out,
            launch_check=False,
            environment=lambda: ENV,
        )


# ------------------------------------------------------------------ (d)/(e) run checks
def _one(res: dict, key: str, **where: Any) -> dict:
    return next(r for r in res["runs"][key] if all(r.get(k) == v for k, v in where.items()))


def test_topas_run_before_the_samples_review_finished_is_invalid_and_fails_the_case(
    tmp_path: Path,
) -> None:
    w = World(
        tmp_path, topas_start="2026-11-20T10:50:00+00:00"
    )  # after the seal review, before the samples one
    res = w.run()
    t = _one(res, "topas", case=1, part=0)
    assert not t["valid"] and any("samples review" in p for p in t["problems"])
    assert res["vc1"]["1"]["status"] == "FAIL" and res["evidence_record"] == ev.EVIDENCE_FAIL
    assert (
        "missing" in res["vc1"]["1"]["reasons"][0]
    )  # the repeat part was never run: FAIL, not omitted


def test_changed_input_hash_invalidates_the_topas_run(tmp_path: Path) -> None:
    w = World(tmp_path)
    w.drop_topas(1, 0)
    w.topas_run(w.samples, 1, 0, tamper_input=True)
    res = w.run()
    t = _one(res, "topas", case=1, part=0)
    assert not t["valid"] and any("differ from the sealed bundle" in p for p in t["problems"])
    assert res["vc1"]["1"]["status"] == "FAIL"


def test_topas_stderr_invalidates_a_part_and_the_repeat_replaces_it_once(tmp_path: Path) -> None:
    w = World(tmp_path)
    w.drop_topas(2, 2)
    w.topas_run(w.samples, 2, 2, stderr="killed track\n")
    w.topas_run(w.samples, 2, 6)
    res = w.run()
    assert res["vc1"]["2"]["status"] == "PASS" and 6 in res["vc1"]["2"]["topas_parts"]
    (tmp_path / "b").mkdir()
    w2 = World(tmp_path / "b")
    w2.drop_topas(2, 2)
    w2.topas_run(w2.samples, 2, 2, stderr="killed track\n")
    w2.topas_run(w2.samples, 2, 6, stderr="killed track\n")
    res2 = w2.run()
    assert res2["vc1"]["2"]["status"] == "FAIL" and "both invalid" in res2["vc1"]["2"]["reasons"][0]


def test_repeat_without_a_trigger_is_a_fail(tmp_path: Path) -> None:
    w = World(tmp_path)
    w.topas_run(w.samples, 3, 5)
    res = w.run()
    assert (
        res["vc1"]["3"]["status"] == "FAIL"
        and "without a validity failure" in res["vc1"]["3"]["reasons"][0]
    )


def test_ionmc_identity_block_mismatch_invalidates_the_run(tmp_path: Path) -> None:
    key = ("vc-sample", 2, MODEL_ICRU, "warp-cpu")
    res = World(tmp_path, bad_identity=key).run()
    r = _one(res, "ionmc", row="vc-sample", case=2, model=MODEL_ICRU)
    assert not r["valid"] and any("identity block" in p for p in r["problems"])
    assert res["vc1"]["2"]["status"] == "FAIL" and res["vc2"]["2"]["verdict"] == "FAIL"
    assert res["integration"]["science_gates"] == "BLOCKED"


def test_differing_warp_modules_invalidate_every_run_of_the_backend(tmp_path: Path) -> None:
    res = World(tmp_path, bad_warp=("vc-sample", 1, MODEL_DEF, "warp-cpu")).run()
    warp = [r for r in res["runs"]["ionmc"] if r["backend"] == "warp-cpu"]
    assert warp and all(not r["valid"] for r in warp)
    assert all(r["valid"] for r in res["runs"]["ionmc"] if r["backend"] == "python")
    assert res["evidence_record"] == ev.EVIDENCE_FAIL


def test_ionmc_run_started_before_the_seal_review_finished_is_invalid(tmp_path: Path) -> None:
    res = World(tmp_path, ion_start="2026-11-20T09:30:00+00:00").run()
    assert all(
        not r["valid"] and any("seal review" in p for p in r["problems"])
        for r in res["runs"]["ionmc"]
    )


def test_unsigned_host_record_invalidates_the_run(tmp_path: Path) -> None:
    res = World(tmp_path, bad_host_sig=("vc2d-sample", 3, MODEL_DEF, "python")).run()
    r = _one(res, "ionmc", row="vc2d-sample", case=3, model=MODEL_DEF, backend="python")
    assert not r["valid"] and any("host-run record" in p for p in r["problems"])
    assert res["vc2d"]["3/python"]["verdict"] == "FAIL"


def test_a_run_at_a_commit_with_a_changed_sealed_blob_is_invalid(tmp_path: Path) -> None:
    w = World(tmp_path)
    put(w.repo, "uv.lock", "changed after the seal\n")
    bad = commit(w.repo, "execution commit with a changed sealed file")
    w.drop_topas(1, 0)
    w.topas_run(bad, 1, 0)
    w.head = bad  # the evaluator commit itself differs too: stops
    w.finish()
    with pytest.raises(ev.RecordError, match="sealed"):
        w.run()


def test_index_must_list_every_ionmc_run(tmp_path: Path) -> None:
    w = World(tmp_path)
    w.ion_index.pop()
    with pytest.raises(ev.RecordError, match="lacks ionmc runs"):
        w.run()


# ------------------------------------------------------------------ arithmetic
def test_td_em_equals_the_r2_td_statistics(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(ev.rio, "case_info", lambda c: {"r_plan_mm": R_PLAN, "ion": "he"})
    rng = np.random.default_rng(1)
    ion = np.array([bragg(ZC) * HZ * (1 + 1e-3 * rng.standard_normal(NZ)) for _ in range(40)])
    parts = [
        (bragg(ZC) * HZ * 25000 * (1 + 2e-3 * rng.standard_normal(NZ)), 25000) for _ in range(4)
    ]
    mine = ev.td_em(ion, parts, HZ, R_PLAN)
    ref = ev.r2e.td_statistics("TD-EM", 1, ion, parts, HZ)
    assert mine["status"] == ref["status"] == "PASS"
    assert [c["verdict"] for c in mine["criteria"]] == [c["verdict"] for c in ref["criteria"]]
    assert mine["pooled"] == ref["pooled"]
    assert mine["criteria"][0]["statistic"] == ref["criteria"][0]["statistic"]
    shifted = [(bragg(ZC, 0.5) * HZ * 25000, 25000)] * 4 + []
    shifted[1:] = [(bragg(ZC, 0.5) * HZ * 25000 * (1 + 1e-3 * i), 25000) for i in range(1, 4)]
    bad = ev.td_em(ion, shifted, HZ, R_PLAN)
    assert bad["status"] == ev.r2e.td_statistics("TD-EM", 1, ion, shifted, HZ)["status"] == "FAIL"


def test_vc2_tost_arithmetic() -> None:
    v = ev.vc2_verdict(0.10, 0.005, 0.10)  # tol = max(0.02, 0.01) = 0.02; 1.645 * 0.005 = 0.0082
    assert v["verdict"] == "PASS" and v["tol_mm"] == 0.02
    assert ev.vc2_verdict(0.10, 0.02, 0.10)["verdict"] == "FAIL"  # 1.645 * 0.02 = 0.0329 > 0.02
    assert ev.vc2_verdict(0.14, 0.0, 0.10)["verdict"] == "FAIL"  # |0.04| > 0.02
    assert ev.vc2_verdict(0.50, 0.001, 0.50)["tol_mm"] == pytest.approx(
        0.05
    )  # 10 % of the prediction
    assert ev.vc2_verdict(None, None, 0.1)["verdict"] == "FAIL"
    assert ev.vc2_verdict(math.nan, 0.0, 0.1)["verdict"] == "FAIL"


def test_vc2d_gate_arithmetic() -> None:
    ok = ev.vc2d_verdict(110.0 + 1e-3, 110.0, 100.0 + 1e-3, 100.0, 10.0, 0.011, 0.01, 0.01)
    assert ok["verdict"] == "PASS" and ok["gate_i"] == ok["gate_ii"] == "PASS"
    assert ev.vc2d_verdict(110.02, 110.0, 100.0, 100.0, 10.0, 0.011, 0.01, 0.01)["gate_i"] == "FAIL"
    assert (
        ev.vc2d_verdict(110.0, 110.0, 100.0, 100.0, 10.05, 0.011, 0.01, 0.01)["gate_ii"] == "FAIL"
    )


def test_jackknife_recovers_a_known_shift_and_sigma_scales_with_noise() -> None:
    rng = np.random.default_rng(3)
    a = np.array([bragg(ZC, 0.3) * HZ * (1 + 5e-3 * rng.standard_normal()) for _ in range(40)])
    b = np.array([bragg(ZC) * HZ * (1 + 5e-3 * rng.standard_normal()) for _ in range(40)])
    j = ev.jackknife_delta_r80(a, b, HZ, R_PLAN)
    assert (
        j["delta_mm"] == pytest.approx(0.3, abs=0.02)
        and 0.0 <= j["sigma_mm"] < 0.01
        and j["n_batches"] == 40
    )
    flat = ev.jackknife_delta_r80(np.zeros((4, NZ)), np.zeros((4, NZ)), HZ, R_PLAN)
    assert flat["delta_mm"] is None
    with pytest.raises(ev.RecordError):
        ev.jackknife_delta_r80(a, b[:10], HZ, R_PLAN)


def test_warp_modules_rule() -> None:
    same = [{"run_id": "a", "backend": "warp-cpu", "warp_modules": {"x": "1"}},
            {"run_id": "b", "backend": "warp-cpu", "warp_modules": {"x": "1"}},
            {"run_id": "c", "backend": "python", "warp_modules": None}]  # fmt: skip
    assert ev.warp_modules_problems(same) == {}
    diff = [*same[:1], {"run_id": "b", "backend": "warp-cpu", "warp_modules": {"x": "2"}}, same[2]]
    assert set(ev.warp_modules_problems(diff)) == {"a", "b"}
    assert set(
        ev.warp_modules_problems([{"run_id": "a", "backend": "warp-cpu", "warp_modules": {}}])
    ) == {"a"}
    assert set(
        ev.warp_modules_problems([{"run_id": "c", "backend": "python", "warp_modules": {"x": "1"}}])
    ) == {"c"}


# ------------------------------------------------------------------ manifest generator
def test_manifest_generator_check_round_trip_and_tamper(world: World) -> None:
    out = world.repo.parent / "m.json"
    assert GEN.main(["--repo", str(world.repo), "--write", "--out", str(out)]) == 0
    assert GEN.main(["--repo", str(world.repo), "--check", "--manifest", str(out)]) == 0
    put(world.repo, "uv.lock", "edited\n")
    assert GEN.main(["--repo", str(world.repo), "--check", "--manifest", str(out)]) == 1
    doc = json.loads(out.read_text())
    assert GEN.MANIFEST_REL not in doc["files"]  # does not list itself
    assert doc["file_count"] == len(doc["files"]) and ev.PREDICTIONS in doc["files"]


def test_generator_and_evaluator_constants_agree() -> None:
    assert GEN.MANIFEST_REL == ev.MANIFEST and GEN.SCHEMA == ev.MANIFEST_SCHEMA
    for name in (
        "EVALUATOR_REL:OWN_REL",
        "LAUNCHER_REL:LAUNCHER_REL",
        "ENV_DIGEST_REL:ENV_DIGEST_REL",
        "STEPS_REL:STEPS_REL",
        "TRUST_REL:TRUST_PATH",
        "ENV_SEAL_REL:ENV_SEAL",
        "PREDICTIONS_REL:PREDICTIONS",
    ):
        g, e = name.split(":")
        assert getattr(GEN, g) == getattr(ev, e)
    assert set(ev.REQUIRED_SEALED) <= set(GEN.FIXED) | {ev.ENV_SEAL}


def test_the_real_closure_contains_the_decision_code() -> None:
    closure = GEN.import_closure(ROOT, [GEN.EVALUATOR_REL, GEN.ENV_DIGEST_REL, GEN.STEPS_REL])
    for rel in ("validation/scripts/stopping/env_seal.py", "validation/scripts/ion_em/r2_io.py",
                "validation/scripts/ion_em/r2_statistics.py", "validation/scripts/ion_em/r2_evaluate.py",
                "infrastructure/validation/review_records.py", "infrastructure/validation/host_run_records.py"):  # fmt: skip
        assert rel in closure


def _cache_dir_pred(blocks: dict[str, dict[str, Any]]) -> dict[str, Any]:
    models = {
        m: {"model_id_effective": m, "version": 1, "spec_sha256": "s" + m, "content_sha256": "c" + m,
            "tables_sha256": {"1": "t"}, "tables_materials": {"1": [{"I_eV_effective": 78.0}]}}
        for m in blocks
    }  # fmt: skip
    return {"models": models, "identity_blocks": {m: {"1": b} for m, b in blocks.items()}}


def _cache_dir_doc(model: str, block: dict[str, Any]) -> dict[str, Any]:
    stop = {"model_id_effective": model, "version": 1, "spec_sha256": "s" + model,
            "content_sha256": "c" + model, **block}  # fmt: skip
    return {"effective_config": {"stopping": stop, "tables": {"sha256": "t", "materials": [{"I_eV_effective": 78.0}]}},
            "stopping_identity": stop}  # fmt: skip


def test_identity_cache_dir_rules_follow_the_sealed_pin() -> None:
    """Codex REVIEW-f5d6f67e finding 1 (operator-authorised fix, IR-20261118-182151-FD0238): the
    analytic default records null and must stay null; a data-backed model must record a real
    directory at the top level and in every data_used entry; the placeholder itself, an empty
    string or null are refused there."""
    ph = ev.CACHE_DIR_PLACEHOLDER
    base = {"model_id_effective": None, "version": 1}
    pred = _cache_dir_pred({
        MODEL_DEF: {**base, "model_id_effective": MODEL_DEF, "spec_sha256": "s" + MODEL_DEF,
                    "content_sha256": "c" + MODEL_DEF, "cache_dir": None, "data_used": []},
        MODEL_ICRU: {**base, "model_id_effective": MODEL_ICRU, "spec_sha256": "s" + MODEL_ICRU,
                     "content_sha256": "c" + MODEL_ICRU, "cache_dir": ph,
                     "data_used": [{"key": "icru90", "cache_dir": ph}]},
    })  # fmt: skip
    real = "/abs/run/cache/ionmc-data"

    def problems(model: str, block: dict[str, Any]) -> list[str]:
        return ev.identity_problems(_cache_dir_doc(model, block), pred, model, 1)

    # analytic default: null passes, any directory fails
    assert problems(MODEL_DEF, {"cache_dir": None, "data_used": []}) == []
    bad = problems(MODEL_DEF, {"cache_dir": real, "data_used": []})
    assert any("analytic model" in p for p in bad) and any("differs" in p for p in bad)
    # data-backed model: a real directory everywhere passes
    assert (
        problems(
            MODEL_ICRU, {"cache_dir": real, "data_used": [{"key": "icru90", "cache_dir": real}]}
        )
        == []
    )
    # ... null, the empty string or the placeholder itself fail at the top level ...
    for v in (None, "", ph):
        bad = problems(
            MODEL_ICRU, {"cache_dir": v, "data_used": [{"key": "icru90", "cache_dir": real}]}
        )
        assert any("the identity block records no cache directory" in p for p in bad), (v, bad)
    # ... and in a nested data_used entry
    for v in (None, "", ph):
        bad = problems(
            MODEL_ICRU, {"cache_dir": real, "data_used": [{"key": "icru90", "cache_dir": v}]}
        )
        assert any("data_used[0] records no cache directory" in p for p in bad), (v, bad)
    # a missing data_used entry is a mismatch, not silently accepted
    bad = problems(MODEL_ICRU, {"cache_dir": real, "data_used": []})
    assert any("differ in number" in p for p in bad)


def test_evaluator_normaliser_matches_the_steps_normaliser_on_null_and_strings() -> None:
    import importlib.util as _iu

    spec = _iu.spec_from_file_location("_steps_v13_for_eval_test", ev.REPO_ROOT / ev.STEPS_REL)
    assert spec is not None and spec.loader is not None
    steps = _iu.module_from_spec(spec)
    spec.loader.exec_module(steps)
    for block in (
        {"cache_dir": None, "data_used": []},
        {
            "cache_dir": "/a/b",
            "data_used": [{"k": 1, "cache_dir": "/a/b"}, {"k": 2, "cache_dir": None}],
        },
        {"model_id_effective": "m"},
    ):
        assert ev._normalise_cache_dir(block) == steps.normalise_cache_dir(block)
