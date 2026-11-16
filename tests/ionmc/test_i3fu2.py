"""Amendment I3-FU-1 step 2: bundle generator, frozen manifest and evaluator (synthetic runs).

Bundles ``validation/reference_cases/topas/i3fu2-*`` and their frozen manifest come from
``validation/scripts/stopping/gen_i3fu2_bundles.py``; the evaluator
``validation/scripts/stopping/i3fu2_evaluate.py`` is exercised on synthetic TOPAS outputs (phase
space with the ``ss_evaluate`` header layout, a TOPAS CSV ``edep.csv``) whose energy losses are
the ionmc prediction times planted factors, and on a synthetic libdEdx archive whose ICRU73_NEW
carbon/oxygen water tables are the ionmc model times a known factor. No real libdEdx, ICRU or
TOPAS data are used.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import re
import shutil
import sys
import tarfile
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from ionmc.materials import water
from ionmc.physics.projectiles import CARBON12, OXYGEN16
from ionmc.physics.stopping import bethe_mass_stopping

ROOT = Path(__file__).resolve().parents[2]
TOPAS = ROOT / "validation" / "reference_cases" / "topas"
MANIFEST = TOPAS / "i3fu2-frozen-manifest.json"
TD_EM = TOPAS / "carbon-water-150mevu-r2-idd-em-part0" / "input.txt"
SCRIPTS = ROOT / "validation" / "scripts"
GRID = np.geomspace(0.025, 1000.0, 53)
SYN_HISTORIES = 2000
TABLE_FACTOR = 1.05  # synthetic ICRU73_NEW = ionmc model x this factor (s ~ +5 %)


def _load(name: str, path: Path) -> ModuleType:
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


gen = _load("gen_i3fu2_bundles", SCRIPTS / "stopping" / "gen_i3fu2_bundles.py")
ev = _load("i3fu2_evaluate", SCRIPTS / "stopping" / "i3fu2_evaluate.py")
gen_ss = _load("gen_ss_bundles", SCRIPTS / "ion_em" / "gen_ss_bundles.py")
seeds = _load("_v3_008_seed_guard", ROOT / "tests" / "ionmc" / "test_v3_008_seed_families.py")

BUNDLES = sorted(d for d in TOPAS.iterdir() if d.is_dir() and d.name.startswith("i3fu2-"))
EVIDENCE = [d for d in BUNDLES if d.name != gen.SMOKE_NAME]


def _native(d: Path, key: str) -> str:
    m = re.findall(rf"^{re.escape(key)} = (.+)$", (d / "input.txt").read_text(), re.M)
    assert len(m) == 1, (d.name, key)
    return str(m[0]).strip()


def _case(d: Path) -> dict[str, Any]:
    return dict(json.loads((d / "case.json").read_text()))


# =============================================================================================
# bundles
# =============================================================================================


def test_bundle_set() -> None:
    assert len(BUNDLES) == 29 and len(EVIDENCE) == 28
    assert [d.name for d in BUNDLES] == gen.bundle_names()
    assert sorted(gen.evidence_keys()) == sorted(
        [(c, 0) for c in range(13, 22)]
        + [(c, 1) for c in (13, 17, 21)]
        + [(c, 2) for c in range(13, 22)]
        + [(c, 3) for c in (13, 17, 21)]
        + [(c, 4) for c in (17, 21)]
        + [(c, 5) for c in (13, 21)]
    )
    assert not set(gen.PARTS) & set(gen.RESERVED_PARTS)
    energies = {c: float(e) for c, e in gen.CASES.items()}
    assert energies == {13: 2.4, 14: 3, 15: 5, 16: 7, 17: 10, 18: 15, 19: 20, 20: 30, 21: 50}
    for d in EVIDENCE:
        m = re.fullmatch(
            r"i3fu2-(carbon|oxygen|he)-water-([0-9p]+)mevu-case(\d\d)-part(\d)", d.name
        )
        assert m, d.name
        info = _case(d)["i3fu2"]
        assert (info["ion"], info["case"], info["part"]) == (
            m.group(1),
            int(m.group(3)),
            int(m.group(4)),
        )
        assert info["kind"] == "evidence"


def test_regeneration_is_byte_identical_and_manifest(tmp_path: Path) -> None:
    made = gen.generate(tmp_path)
    assert sorted(nm for nm, _ in made) == [d.name for d in BUNDLES]
    for nm, _ in made:
        for f in ("case.json", "input.txt"):
            assert (TOPAS / nm / f).read_bytes() == (tmp_path / nm / f).read_bytes(), nm
    assert gen.manifest_text(tmp_path) == MANIFEST.read_text()
    assert gen.check(TOPAS, MANIFEST)
    m = json.loads(MANIFEST.read_text())
    assert m["schema"] == "ionmc-i3fu2-frozen-bundles-1"
    assert m["bundle_count"] == 29 and m["file_count"] == 58
    for d in BUNDLES:
        for f in ("case.json", "input.txt"):
            assert m["bundles"][d.name][f] == hashlib.sha256((d / f).read_bytes()).hexdigest()
    for p in [MANIFEST] + [d / f for d in BUNDLES for f in ("case.json", "input.txt")]:
        b = p.read_bytes()
        assert b.endswith(b"\n") and not b.endswith(b"\n\n"), p


def test_seeds_formula_and_reservation() -> None:
    topas = seeds.topas_seeds()
    used = []
    for d in EVIDENCE:
        c = _case(d)
        case_no, part = c["i3fu2"]["case"], c["i3fu2"]["part"]
        s = seeds.TOPAS_8_BASE + 100 * part + case_no
        assert c["seeds"] == [s] == [gen.topas_seed(case_no, part)]
        assert int(_native(d, "i:Ts/Seed")) == s
        assert s in topas and s in seeds.V3_008_RESERVED
        assert case_no in seeds.TOPAS_CASES and part in seeds.TOPAS_PARTS
        used.append(s)
    assert len(set(used)) == 28 and min(used) == 20619013 and max(used) == 20619521
    smoke = _case(TOPAS / gen.SMOKE_NAME)
    assert smoke["seeds"] == [20600440] and smoke["histories"] == 2000
    assert 20600440 in seeds.REHEARSAL_8 and 20600440 in seeds.V3_008_RESERVED
    assert smoke["i3fu2"]["kind"] == "smoke" and "not evidence" in smoke["rationale"]


def test_thickness_rule_and_exit_energy() -> None:
    proj = {"carbon": CARBON12, "oxygen": OXYGEN16}
    for d in BUNDLES:
        info = _case(d)["i3fu2"]
        ion, a, e0, t = info["ion"], info["a"], info["e0_per_u_mev"], info["slab_mm"]
        p = gen.IONS[ion][0]
        if ion in proj:
            assert p is proj[ion]
        s = float(bethe_mass_stopping(e0, p, water(78.0))[0])
        t_full = float(f"{10 * 0.03 * e0 * a / s:.3g}")
        half = info["part"] in (1, 3)
        assert t == (float(f"{t_full / 2:.3g}") if half else t_full), d.name
        frac = 10 * 0.03 * e0 * a / s
        assert abs(frac / t_full - 1.0) <= 0.005, d.name  # 3 % to within 3-digit rounding
        de = ev.slab_energy_loss(a * e0, a, t, ev.model_stopping(p))  # raises below 2 MeV/u
        assert (a * e0 - de) / a > 2.0, d.name
        assert abs(de / (a * e0) - (0.015 if half else 0.03)) < 0.003, d.name
        hlz = float(_native(d, "d:Ge/Water/HLZ").removesuffix(" mm"))
        assert _native(d, "d:Ge/Water/TransZ") == _native(d, "d:Ge/Water/HLZ")
        assert abs(2 * hlz - t) < 1e-12, d.name
        assert re.fullmatch(r"\d+\.\d{1,6} mm", _native(d, "d:Ge/Water/HLZ"))


def test_physics_block_identical_to_td_em() -> None:
    def ph(text: str) -> list[str]:
        return [ln for ln in text.splitlines() if re.match(r"^\w+:Ph/", ln)]

    td = ph(TD_EM.read_text())
    assert td == gen.PHYSICS_LINES
    for d in BUNDLES:
        lines = ph((d / "input.txt").read_text())
        part = _case(d)["i3fu2"]["part"]
        extra = [gen.PH_VERBOSITY_LINE]
        if part == 5 and d.name != gen.SMOKE_NAME:
            extra.append(gen.FLUCTUATIONS_OFF_LINE)
            assert _case(d)["i3fu2"]["fluctuations"] is False
        else:
            assert _case(d)["i3fu2"]["fluctuations"] is True
        assert sorted(ln for ln in lines if ln not in td) == sorted(extra), d.name
        assert [ln for ln in lines if ln in td] == td, d.name


def test_native_lines_scorers_and_case_keys() -> None:
    r2_keys = list(_case(TOPAS / "he-water-150mevu-r2-lat-em-part4"))
    assert gen.SSEXIT_LINES == gen_ss.SCORER_LINES
    for d in BUNDLES:
        c = _case(d)
        assert list(c) == r2_keys + ["i3fu2"]
        assert c["expected_outputs"] == ["ss_exit.header", "ss_exit.phsp", "edep.csv"]
        assert c["histories"] == (2000 if d.name == gen.SMOKE_NAME else 50000)
        assert int(_native(d, "i:So/Beam/NumberOfHistoriesInRun")) == c["histories"]
        info = c["i3fu2"]
        _p, z, a, particle = gen.IONS[info["ion"]]
        assert (info["z"], info["a"]) == (z, a)
        assert _native(d, "s:So/Beam/BeamParticle") == f'"{particle}"'
        energy = _native(d, "d:So/Beam/BeamEnergy")
        assert abs(float(energy.removesuffix(" MeV")) - a * info["e0_per_u_mev"]) < 1e-9
        assert re.fullmatch(r"\d+(\.\d+)? MeV", energy) and not energy.startswith("0")
        for key, want in (
            ("i:Ts/MaxStepNumber", "10000000"),
            ("i:Ts/RunVerbosity", "2"),
            ("i:Ph/Verbosity", "1"),
            ("i:Ts/TrackingVerbosity", "0"),
            ("i:Ts/NumberOfThreads", "1"),
            ("b:Ts/PauseBeforeQuit", '"False"'),
            ("b:Ts/ShowCPUTime", '"True"'),
            ("d:Ge/Source/TransZ", "-0.001 mm"),
            ("u:So/Beam/BeamEnergySpread", "0"),
            ("d:Ge/World/HLZ", "40 cm"),
            ("d:Ge/Water/HLX", "60 mm"),
        ):
            assert _native(d, key) == want, (d.name, key)
        text = (d / "input.txt").read_text()
        assert [
            ln for ln in text.splitlines() if ln.startswith(("s:Sc", "b:Sc", "i:Sc", "sv:Sc"))
        ] == (gen.SSEXIT_LINES + gen.EDEP_LINES)
        head = text.splitlines()[0]
        assert "I3-FU-1 step 2" in head and "no validation claim" in head
        if d.name != gen.SMOKE_NAME:
            assert "20619000 + 100*part + case" in head
            assert "diagnostic, no validation claim" in c["rationale"]


# =============================================================================================
# evaluator on synthetic runs
# =============================================================================================

HEADER_TEMPLATE = """TOPAS ASCII Phase Space

Number of Original Histories: {n}
Number of Original Histories that Reached Phase Space: {reached}
Number of Scored Particles: {scored}

Columns of data are as follows:
 1: Position X [cm]
 2: Position Y [cm]
 3: Position Z [cm]
 4: Direction Cosine X
 5: Direction Cosine Y
 6: Energy [MeV]
 7: Weight
 8: Particle Type (in PDG Format)
 9: Flag to tell if Third Direction Cosine is Negative (1 means true)
10: Flag to tell if this is the First Scored Particle from this History (1 means true)

"""
EDEP_TEMPLATE = """# TOPAS Version: 4.3
# Parameter File: input.txt
# Results for scorer: EDep
# Scored in component: Water
# X in 1 bin  of 12 cm
# Y in 1 bin  of 12 cm
# Z in 1 bin  of {z_cm} cm
# EnergyDeposit ( MeV ) : Sum   Standard_Deviation
0, 0, 0, {total:.10g}, {sd:.6g}
"""


def _fixture_header_layout() -> None:
    fx = (ROOT / "tests" / "ionmc" / "data" / "v5ss" / "ss_exit.header").read_text()
    syn = HEADER_TEMPLATE.format(n=3, reached=3, scored=3)
    assert ev.ss.parse_header(fx)["columns"] == ev.ss.parse_header(syn)["columns"]


def test_synthetic_header_matches_fixture_layout() -> None:
    _fixture_header_layout()


def _archive(tmp: Path, factor: float) -> tuple[Path, str, dict[str, str]]:
    """Synthetic libdEdx archive: ICRU73_NEW water C/O = ionmc model x ``factor`` (E >= 1)."""
    energy_text = f"{GRID.size}\n" + "".join(f"{e:.10g}\n" for e in GRID)
    blocks = []
    for target in (104, 276):
        for z in range(1, 10):
            if target == 276 and z in (6, 8):
                proj = {6: CARBON12, 8: OXYGEN16}[z]
                s = bethe_mass_stopping(np.maximum(GRID, 1.0), proj, water(78.0)) * factor
            else:
                s = np.full(GRID.size, 1.0e6 + z)  # decoys
            blocks.append(f"#{target}:{z:02d}:{GRID.size}\n" + "".join(f"{v:.12E}\n" for v in s))
    files = {
        "libdedx-main/data/raw/icru73_newEng.dat": energy_text.encode(),
        "libdedx-main/data/raw/ICRU73_NEW.dat": "\n\n".join(blocks).encode(),
    }
    path = tmp / "archive.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    members = {
        m: hashlib.sha256(files[f"libdedx-main/{m}"]).hexdigest()
        for m in ("data/raw/ICRU73_NEW.dat", "data/raw/icru73_newEng.dat")
    }
    return path, hashlib.sha256(path.read_bytes()).hexdigest(), members


@pytest.fixture(scope="module")
def syn(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Synthetic bundle set (2000 histories), its manifest, archive and the model predictions."""
    root = tmp_path_factory.mktemp("i3fu2")
    bundles = root / "bundles"
    gen.generate(bundles, histories=SYN_HISTORIES)
    manifest = root / "manifest.json"
    manifest.write_text(gen.manifest_text(bundles))
    archive, sha, members = _archive(root, TABLE_FACTOR)
    de_m = {}
    for c, p in gen.evidence_keys():
        info = json.loads((bundles / gen.bundle_name(c, p) / "case.json").read_text())["i3fu2"]
        a, e0 = info["a"], info["e0_per_u_mev"]
        stop = ev.model_stopping(gen.IONS[info["ion"]][0])
        de_m[(c, p)] = ev.slab_energy_loss(a * e0, a, info["slab_mm"], stop)
    return {
        "root": root,
        "bundles": bundles,
        "manifest": manifest,
        "archive": archive,
        "sha": sha,
        "members": members,
        "de_m": de_m,
    }


def _write_run(
    run_dir: Path,
    case_doc: dict[str, Any],
    de_m: float,
    r: float,
    rng: np.random.Generator,
    *,
    deficit: int = 0,
    edep_scale: float = 1.0,
) -> None:
    info = case_doc["i3fu2"]
    n = int(case_doc["histories"])
    a, e0, t = info["a"], info["e0_per_u_mev"], info["slab_mm"]
    e_in = a * e0
    de_rows = de_m * (1.0 + r) * (1.0 + 0.02 * rng.standard_normal(n - deficit))
    e_exit = e_in - de_rows
    pdg = ev.PDG[info["ion"]]
    rows = np.zeros((n - deficit, 10))
    rows[:, 2] = t / 10.0
    rows[:, 5] = e_exit
    rows[:, 6] = 1.0
    rows[:, 7] = pdg
    rows[:, 9] = 1.0
    buf = io.StringIO()
    np.savetxt(buf, rows, fmt=["%.6g"] * 5 + ["%.10g", "%g", "%d", "%d", "%d"])
    run_dir.mkdir(parents=True)
    (run_dir / "ss_exit.header").write_text(
        HEADER_TEMPLATE.format(n=n, reached=n - deficit, scored=n - deficit)
    )
    (run_dir / "ss_exit.phsp").write_text(buf.getvalue())
    total = float(np.sum(de_rows)) * edep_scale
    (run_dir / "edep.csv").write_text(
        EDEP_TEMPLATE.format(z_cm=t / 10.0, total=total, sd=0.01 * total)
    )
    (run_dir / "record.json").write_text(json.dumps({"run_id": f"REF-syn-{run_dir.name}"}))


def _scenario(
    syn: dict[str, Any],
    name: str,
    r_of: Any,
    *,
    skip: tuple[tuple[int, int], ...] = (),
) -> Path:
    runs = Path(syn["root"]) / f"runs-{name}"
    for c, p in gen.evidence_keys():
        if (c, p) in skip:
            continue
        bname = gen.bundle_name(c, p)
        doc = json.loads((syn["bundles"] / bname / "case.json").read_text())
        rng = np.random.default_rng([20261116, c, p])  # per run: independent of skipped runs
        _write_run(runs / bname, doc, syn["de_m"][(c, p)], r_of(c, p), rng)
    return runs


def _evaluate(syn: dict[str, Any], runs: Path) -> dict[str, Any]:
    doc = ev.evaluate(
        runs,
        syn["archive"],
        syn["manifest"],
        "0" * 40,
        syn["bundles"],
        archive_sha256=syn["sha"],
        member_sha256=syn["members"],
    )
    return dict(doc)


def _r_icru(c: int, p: int) -> float:
    return 0.0 if p == 4 else TABLE_FACTOR - 1.0


def _r_model(c: int, p: int) -> float:
    return 0.0


def _r_mixed(c: int, p: int) -> float:
    """C sides with the table, O with the model; carbon thickness, He and fluctuation controls
    planted to fail; the oxygen thickness control passes."""
    return {0: 0.05, 1: 0.06, 2: 0.0, 3: 0.0, 4: 0.02, 5: 0.056}[p]


@pytest.fixture(scope="module")
def results(syn: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out = {}
    for name, fn in (("icru", _r_icru), ("model", _r_model), ("mixed", _r_mixed)):
        out[name] = _evaluate(syn, _scenario(syn, name, fn))
    out["missing"] = _evaluate(syn, _scenario(syn, "missing", _r_icru, skip=((21, 4), (13, 1))))
    return out


def _bands(doc: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return dict(doc["bands"])


def test_icru_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["icru"]
    assert doc["pattern"] == "ICRU"
    for ion in ("carbon", "oxygen"):
        for b, v in _bands(doc)[ion].items():
            assert v["class"] == "I", (ion, b)
            assert v["mean_r_M"] == pytest.approx(0.05, abs=1.5e-3)
            assert v["mean_r_I"] == pytest.approx(0.0, abs=4e-3)
            assert v["mean_s"] == pytest.approx(0.05, abs=4e-3)
            assert v["x_b"] == 0.01
            assert v["n"] == {"[2,5)": 2, "[5,10)": 2, "[10,20)": 2, "[20,50]": 3}[b]
            assert v["sigma_b"] == pytest.approx(np.hypot(v["se_b"], 0.0025))
            assert 1e-4 < v["se_b"] < 1e-3  # 2 % row noise, 2000 rows, n runs
    c = doc["controls"]
    assert c["controls_pass"] is True and doc["conclusive"] is True
    for ion in ("carbon", "oxygen"):
        assert len(c["thickness"][ion]["differences"]) == 3 and c["thickness"][ion]["pass"]
    assert len(c["helium"]["r_M"]) == 2 and c["helium"]["pass"]
    assert all(abs(v) < 2e-3 for v in c["helium"]["r_M"])
    assert len(c["fluctuations_off"]["differences"]) == 2 and c["fluctuations_off"]["pass"]
    assert c["missing_runs"] == []
    runs = doc["provenance"]["runs"]
    assert len(runs) == 28
    assert all(
        v["scored"] == SYN_HISTORIES and v["run_id"].startswith("REF-syn-") for v in runs.values()
    )


def test_model_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["model"]
    assert doc["pattern"] == "MODEL"
    for ion in ("carbon", "oxygen"):
        for b, v in _bands(doc)[ion].items():
            assert v["class"] == "M", (ion, b)
            assert v["mean_r_M"] == pytest.approx(0.0, abs=1.5e-3)
            assert v["mean_r_I"] == pytest.approx(-0.05 / 1.05, abs=4e-3)
    assert doc["controls"]["controls_pass"] is True and doc["conclusive"] is True


def test_mixed_scenario_and_failing_controls(results: dict[str, dict[str, Any]]) -> None:
    doc = results["mixed"]
    assert doc["pattern"] == "MIXED"
    assert all(v["class"] == "I" for v in _bands(doc)["carbon"].values())
    assert all(v["class"] == "M" for v in _bands(doc)["oxygen"].values())
    c = doc["controls"]
    assert c["thickness"]["carbon"]["pass"] is False
    assert all(d == pytest.approx(-0.01, abs=2e-3) for d in c["thickness"]["carbon"]["differences"])
    assert c["thickness"]["oxygen"]["pass"] is True
    assert c["helium"]["pass"] is False
    assert all(v == pytest.approx(0.02, abs=2e-3) for v in c["helium"]["r_M"])
    assert c["fluctuations_off"]["pass"] is False
    assert all(d == pytest.approx(0.006, abs=2e-3) for d in c["fluctuations_off"]["differences"])
    assert c["controls_pass"] is False and doc["conclusive"] is False


def test_missing_control_run_reports_bands_but_not_conclusive(
    results: dict[str, dict[str, Any]],
) -> None:
    doc = results["missing"]
    assert doc["bands"] == results["icru"]["bands"]  # same per-run noise as the ICRU scenario
    assert doc["pattern"] == "ICRU"
    c = doc["controls"]
    assert sorted(c["missing_runs"]) == sorted([gen.bundle_name(21, 4), gen.bundle_name(13, 1)])
    assert c["helium"]["complete"] is False and c["helium"]["pass"] is False
    assert c["thickness"]["carbon"]["complete"] is False
    assert len(c["thickness"]["carbon"]["differences"]) == 2
    assert c["controls_pass"] is False and doc["conclusive"] is False


def _walk(x: Any, path: str = "") -> list[tuple[str, Any]]:
    out = [(path, x)]
    if isinstance(x, dict):
        for k, v in x.items():
            out += _walk(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            out += _walk(v, f"{path}[{i}]")
    return out


def test_no_per_energy_leak(results: dict[str, dict[str, Any]]) -> None:
    band_keys = {"n", "mean_r_M", "mean_r_I", "mean_s", "se_b", "sigma_b", "x_b", "class"}
    for doc in results.values():
        for path, v in _walk(doc):
            if isinstance(v, list):
                assert len(v) <= 3, path
        assert set(doc["bands"]) == {"carbon", "oxygen"}
        for ion in doc["bands"].values():
            assert set(ion) == {"[2,5)", "[5,10)", "[10,20)", "[20,50]"}
            for band in ion.values():
                assert set(band) == band_keys
        for run in doc["provenance"]["runs"].values():
            assert set(run) == {"seed", "histories", "scored", "run_id"}
        text = json.dumps(doc)
        for bad in ("r_I_per", "delta_e", "de_m", "de_i", "exit_energy"):
            assert bad not in text.lower()
        assert set(doc["controls"]["helium"]) == {
            "e0_per_u_mev",
            "r_M",
            "tolerance",
            "complete",
            "pass",
        }


def test_classify_rule() -> None:
    sig = 0.0025
    k = 1.645 * sig
    assert ev.classify(0.0, 0.0, 0.003, sig) == "U"  # X = 0.001 <= k
    assert ev.classify(0.05, 0.0, 0.05, sig) == "I"
    assert ev.classify(0.0, -0.05, 0.05, sig) == "M"
    assert ev.classify(-0.03, -0.076, 0.05, sig) == "L"
    assert ev.classify(0.02, -0.03, 0.05, sig) == "N"
    assert ev.classify(0.0, 0.0, 3 * k, sig) == "U"  # X == k: U checked first
    assert ev.band_of(2.4) == "[2,5)" and ev.band_of(5.0) == "[5,10)"
    assert ev.band_of(20.0) == "[20,50]" and ev.band_of(50.0) == "[20,50]"
    with pytest.raises(ValueError):
        ev.band_of(50.1)


def test_rk4_against_constant_stopping_and_domain_guard() -> None:
    assert ev.slab_energy_loss(120.0, 12, 0.1, lambda e: 500.0) == pytest.approx(5.0, rel=1e-12)
    with pytest.raises(ValueError, match="200"):
        ev.slab_energy_loss(120.0, 12, 0.1, lambda e: 500.0, n_steps=100)
    with pytest.raises(ValueError, match="domain"):
        ev.slab_energy_loss(28.8, 12, 1.0, lambda e: 500.0)  # 28.8 - 50 MeV < 2 MeV/u
    s = ev.table_stopping(GRID, 100.0 * GRID**-0.8)
    assert s(10.0) == pytest.approx(100.0 * 10.0**-0.8, rel=1e-12)


# --- fail-closed paths ---------------------------------------------------------------------


def _one(syn: dict[str, Any], tmp: Path, **kw: Any) -> tuple[Path, dict[str, Any]]:
    bname = gen.bundle_name(17, 0)
    doc = json.loads((syn["bundles"] / bname / "case.json").read_text())
    _write_run(tmp / bname, doc, syn["de_m"][(17, 0)], 0.0, np.random.default_rng(1), **kw)
    return tmp / bname, doc


def test_measure_run_validity(syn: dict[str, Any], tmp_path: Path) -> None:
    run, doc = _one(syn, tmp_path / "ok")
    m = ev.measure_run(run, doc)
    assert m["de"] == pytest.approx(syn["de_m"][(17, 0)], rel=3e-3)
    run, doc = _one(syn, tmp_path / "deficit", deficit=1)
    with pytest.raises(ev.RecordError, match="deficit"):
        ev.measure_run(run, doc)
    run, doc = _one(syn, tmp_path / "edep", edep_scale=1.002)
    with pytest.raises(ev.RecordError, match="cross check"):
        ev.measure_run(run, doc)
    run, doc = _one(syn, tmp_path / "pdg")
    other = dict(doc, i3fu2=dict(doc["i3fu2"], ion="oxygen"))
    with pytest.raises(ev.RecordError, match="PDG"):
        ev.measure_run(run, other)
    run, doc = _one(syn, tmp_path / "z")
    thick = dict(doc, i3fu2=dict(doc["i3fu2"], slab_mm=2 * doc["i3fu2"]["slab_mm"]))
    with pytest.raises(ev.RecordError, match="far face"):
        ev.measure_run(run, thick)
    run, doc = _one(syn, tmp_path / "missing")
    (run / "edep.csv").unlink()
    with pytest.raises(ev.RecordError, match="missing output"):
        ev.measure_run(run, doc)


def test_missing_evidence_run_and_manifest_mismatch(syn: dict[str, Any], tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    runs.mkdir()
    with pytest.raises(ev.RecordError, match="missing evidence run"):
        _evaluate(syn, runs)
    bad = tmp_path / "bundles"
    shutil.copytree(syn["bundles"], bad)
    p = bad / gen.bundle_name(13, 0) / "input.txt"
    p.write_text(p.read_text().replace("i:Ts/Seed = 20619013", "i:Ts/Seed = 20619014"))
    with pytest.raises(ev.RecordError, match="manifest"):
        ev.evaluate(
            runs,
            syn["archive"],
            syn["manifest"],
            "0" * 40,
            bad,
            archive_sha256=syn["sha"],
            member_sha256=syn["members"],
        )
    with pytest.raises(ev.RecordError, match="schema"):
        ev.evaluate(runs, syn["archive"], MANIFEST.parent / "v5ss-frozen-manifest.json", "x")


def test_run_input_copy_must_match_manifest(syn: dict[str, Any], tmp_path: Path) -> None:
    run, _doc = _one(syn, tmp_path)
    entries = json.loads(syn["manifest"].read_text())["bundles"]
    shutil.copy(syn["bundles"] / run.name / "input.txt", run / "input.txt")
    ev.check_run_inputs(run, entries[run.name])
    (run / "input.txt").write_text("changed\n")
    with pytest.raises(ev.RecordError, match="frozen manifest"):
        ev.check_run_inputs(run, entries[run.name])
