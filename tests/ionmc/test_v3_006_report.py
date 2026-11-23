"""Rows REP-1..REP-4 of plan v3-006 (step S5): ``ionmc.svg``, ``ionmc.report`` and
``regenerate_report``.

Two real tiny transport runs (8 histories, python backend, single process) are written with
``ionmc.results.save``; the golden fixture is the committed synthetic bundle under
``data/v3_006/golden`` (``result.json``, ``report.md``, ``plots/*.svg`` and the arrays as ``.npy``
files; REG-4 forbids a committed ``.npz``, so the test assembles ``arrays.npz`` with the writer of
``ionmc.results``). No physics table value and no stopping data appear in the fixture.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc import report, results, svg
from ionmc.config import SimulationConfig
from ionmc.configio import load_toml
from ionmc.physics.stopping_evidence import STOPPING_EVIDENCE
from ionmc.report import regenerate_report, render_plots, render_report
from ionmc.results import load, save
from ionmc.scoring import TallyRequest
from ionmc.simulation import Result, Simulation

DATA = Path(__file__).resolve().parent / "data" / "v3_006"
GOLDEN = DATA / "golden"
SRC = Path(__file__).resolve().parents[2] / "src"
BOX_TOML = DATA / "proton_box.toml"
RESERVED = re.compile(r"\bvalidated\b")
FENCE = re.compile(r"```evidence\n(.*?)\n```\n", re.DOTALL)


def _run(cfg: SimulationConfig) -> Result:
    assert cfg.run.cpu_workers == 1
    return Simulation(cfg).run()


@pytest.fixture(scope="module")
def toml_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A TOML-configured run with its manifest (a config file path is a value of the set P)."""
    cfg, man = load_toml(BOX_TOML)
    return save(_run(cfg), tmp_path_factory.mktemp("toml") / "bundle", manifest=man)


@pytest.fixture(scope="module")
def tilted_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """An API run with a tilted source (dominant axis x) and a LET tally."""
    cfg, _ = load_toml(BOX_TOML)
    src = dataclasses.replace(cfg.source, position_mm=(-14.0, 0.0, 10.0), direction=(1.0, 0.2, 0.1))
    cfg = dataclasses.replace(cfg, source=src, tallies=(TallyRequest("let", "dose", "let_d"),))
    return save(_run(cfg), tmp_path_factory.mktemp("tilt") / "bundle")


def _assemble_golden(dst: Path) -> Path:
    """The committed golden bundle with ``arrays.npz`` assembled from its ``.npy`` files by the
    deterministic writer (stored, ``zlib_level = 0``, as the generator wrote it)."""
    shutil.copytree(GOLDEN, dst, ignore=shutil.ignore_patterns("arrays"))
    doc = json.loads((dst / "result.json").read_text(encoding="utf-8"))
    arrs = []
    for e in doc["arrays"]:
        a = np.load(GOLDEN / "arrays" / f"{e['name']}.npy", allow_pickle=False)
        arrs.append(
            results._Arr(
                e["name"], a, tuple(e["axes"]), e["grid"], e["units"], e["definition"], e["origin"]
            )
        )
    results.write_npz(dst / "arrays.npz", arrs, 0)
    sha, size = results._sha256_path(dst / "arrays.npz")
    if doc["files"]["arrays.npz"] != {"sha256": sha, "size": size}:
        # the file digest of a zip depends on the zip library; the content digests of every
        # array (checked by load) are the authoritative identity (decision 0045 section 1)
        doc["files"]["arrays.npz"] = {"sha256": sha, "size": size}
        (dst / "result.json").write_bytes(results.canonical_bytes(doc))
    return dst


@pytest.fixture(scope="module")
def golden(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _assemble_golden(tmp_path_factory.mktemp("golden") / "bundle")


def _derived(bundle: Path) -> dict[str, bytes]:
    files = json.loads((bundle / "result.json").read_text(encoding="utf-8"))["files"]
    return {n: (bundle / n).read_bytes() for n in files if n != "arrays.npz"}


def _generated(text: str) -> str:
    """The report without its verbatim evidence records (decision 0045 section 5)."""
    return FENCE.sub("", text)


# ---- REP-1 ------------------------------------------------------------------------------------

_CHILD = """
import hashlib, json, sys
sys.path.insert(0, sys.argv[1])
from ionmc.report import regenerate_report
out = regenerate_report(sys.argv[2])
print("REGEN:" + json.dumps({n: hashlib.sha256(b).hexdigest() for n, b in sorted(out.items())}))
"""


def _clean_process_digests(bundle: Path) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    env.update(IONMC_SINGLE_PROCESS="1", OMP_NUM_THREADS="1")
    p = subprocess.run(
        [sys.executable, "-I", "-c", _CHILD, str(SRC), str(bundle)],
        capture_output=True, text=True, timeout=300, env=env, check=False, cwd=bundle,
    )  # fmt: skip
    assert p.returncode == 0, p.stderr[-2000:]
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("REGEN:")]
    assert len(lines) == 1
    return json.loads(lines[0][6:])  # type: ignore[no-any-return]


def test_rep1_regeneration_equals_written_files(
    toml_bundle: Path, tilted_bundle: Path, golden: Path
) -> None:
    for b in (toml_bundle, tilted_bundle, golden):
        written = _derived(b)
        assert "report.md" in written
        # in this process
        assert regenerate_report(b) == written
        # in a clean interpreter (python -I)
        digests = _clean_process_digests(b)
        assert digests == {n: hashlib.sha256(v).hexdigest() for n, v in sorted(written.items())}
    assert any(n.startswith("plots/") for n in _derived(tilted_bundle))


def test_rep1_save_lists_derived_files_and_tree_is_strict(tilted_bundle: Path) -> None:
    doc = json.loads((tilted_bundle / "result.json").read_text(encoding="utf-8"))
    names = sorted(doc["files"])
    assert names[0] == "arrays.npz" and "report.md" in names
    assert all(n.startswith("plots/") and n.endswith(".svg") for n in names if "/" in n)
    assert sorted(p.name for p in tilted_bundle.iterdir()) == [
        "arrays.npz", "plots", "report.md", "result.json",
    ]  # fmt: skip
    assert sorted(f"plots/{p.name}" for p in (tilted_bundle / "plots").iterdir()) == [
        n for n in names if n.startswith("plots/")
    ]
    assert load(tilted_bundle).doc["files"] == doc["files"]


# ---- REP-2 ------------------------------------------------------------------------------------


def test_rep2_golden_bytes(golden: Path) -> None:
    loaded = load(golden)  # a complete valid 1.0 bundle: tier 1 passes
    out = regenerate_report(golden)
    committed = {
        str(p.relative_to(GOLDEN)): p.read_bytes()
        for p in sorted(GOLDEN.rglob("*"))
        if p.is_file() and p.relative_to(GOLDEN).parts[0] in ("report.md", "plots")
    }
    assert out == committed
    assert sorted(out) == sorted(n for n in loaded.doc["files"] if n != "arrays.npz")
    for name, rec in loaded.doc["files"].items():
        if name in committed:
            assert hashlib.sha256(committed[name]).hexdigest() == rec["sha256"], name
    assert loaded.doc["effective_extra"]["profile_axis"] == 1


def test_rep2_golden_is_kb_scale_synthetic_and_has_no_npz() -> None:
    files = [p for p in GOLDEN.rglob("*") if p.is_file()]
    assert files
    for p in files:
        assert p.stat().st_size < 64 * 1024, p
        assert p.suffix != ".npz", p
        if p.suffix == ".npy":
            assert np.load(p, allow_pickle=False).size < 2000, p
    doc = json.loads((GOLDEN / "result.json").read_text(encoding="utf-8"))
    assert doc["run"]["output_path"] == "golden"  # no machine path
    assert doc["provenance"]["code"]["source"] == "unavailable"  # a synthetic code identity
    assert doc["provenance"]["data"]["stopping_data_used"] == []  # no stopping data
    assert not (GOLDEN / "arrays.npz").exists()


def test_rep2_generator_reproduces_the_committed_fixture(tmp_path: Path) -> None:
    spec = importlib.util.spec_from_file_location("gen_golden", DATA / "gen_golden.py")
    assert spec is not None and spec.loader is not None
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    out = gen.build(tmp_path / "golden")
    mine = sorted(p.relative_to(out) for p in out.rglob("*") if p.is_file())
    theirs = sorted(p.relative_to(GOLDEN) for p in GOLDEN.rglob("*") if p.is_file())
    assert mine == theirs
    # Byte identity holds in the local-validation tier (IONMC_LV_TIER=1). Elsewhere (the hosted CI
    # runner) the stopping tables built by `validate` differ in their last bits (plan v3-006, REG-2
    # clarification), so only the table-derived pointers of result.json and the report lines that
    # print them may differ; every other pointer, line, array and SVG is compared exactly.
    lv = os.environ.get("IONMC_LV_TIER") == "1"
    for rel in theirs:
        a, b = (out / rel).read_bytes(), (GOLDEN / rel).read_bytes()
        if lv or rel.suffix not in (".json", ".md"):
            assert a == b, rel
        elif rel.suffix == ".json":
            diffs = _diff_pointers(json.loads(a), json.loads(b))
            foreign = [d for d in diffs if not d.startswith(TABLE_DERIVED_POINTERS)]
            assert not foreign, f"{rel}: platform-independent pointers differ: {foreign}"
        else:
            keep_a = [ln for ln in a.decode("utf-8").splitlines() if not _table_derived_line(ln)]
            keep_b = [ln for ln in b.decode("utf-8").splitlines() if not _table_derived_line(ln)]
            assert keep_a == keep_b, rel


# result.json pointers whose values derive from the stopping tables (digests, ranges, step and
# accumulator bounds computed from them) and the report lines that print such values
TABLE_DERIVED_POINTERS = (
    "/effective/tables", "/effective/max_steps", "/effective/tallies/channels",
    "/provenance/data/table_sha256", "/provenance/data/table_materials",
    "/files/report.md/sha256",
)  # fmt: skip
_TABLE_DERIVED_LINE = re.compile(
    r"table sha256:|^\| max-steps \||^\| tallies \||^\| stopping-identity \|"
)


def _table_derived_line(line: str) -> bool:
    return _TABLE_DERIVED_LINE.search(line) is not None


def _diff_pointers(a: Any, b: Any, path: str = "") -> list[str]:
    """Every JSON pointer at which ``a`` and ``b`` differ (a differing container is reported once
    per differing leaf; a length or key-set difference as the container's pointer)."""
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return [path or "/"]
        return [d for k in a for d in _diff_pointers(a[k], b[k], f"{path}/{k}")]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [path or "/"]
        return [d for i, (x, y) in enumerate(zip(a, b, strict=True))
                for d in _diff_pointers(x, y, f"{path}/{i}")]  # fmt: skip
    return [] if a == b else [path or "/"]


# ---- REP-3 ------------------------------------------------------------------------------------


def test_rep3_content_of_a_real_run(toml_bundle: Path, tmp_path: Path) -> None:
    loaded = load(toml_bundle)
    doc = loaded.doc
    text = (toml_bundle / "report.md").read_text(encoding="utf-8")
    code = doc["provenance"]["code"]
    for needle in (
        "## Identity", "## Validity", "## Requested versus effective configuration",
        "## Source, geometry and materials", "## Physics and stopping model",
        "## Data identities and roles", "## Run and random numbers", "## Statistics definition",
        "## Per-grid summary", "## Energy balance and counters", "## Missing items",
        f"package tree object id: {code['package_tree_oid']}",
        f"commit: {code['git_commit'] or 'unknown'}",
        "**Result status: valid**", "| row | class | rule | requested | effective | reason |",
        "normalisation: per primary history", "uncertainty: batch standard error",
        "closure_residual_mev", "counter step_truncation", "proton_box.toml",
    ):  # fmt: skip
        assert needle in text, needle
    assert text.splitlines()[-1] == report.LAST_LINE == "Not a clinical claim."
    assert text.endswith("Not a clinical claim.\n") and not text.endswith("\n\n")
    # every evidence_at_run record verbatim, in order
    blocks = [json.loads(b) for b in FENCE.findall(text)]
    assert blocks == doc["provenance"]["models"]["evidence_at_run"] and blocks
    # input files by basename and sha256; no timing, no path of the set P, no absolute path
    assert doc["input_manifest"]["config_file"]["sha256"] in text
    p_values = [
        doc["run"]["output_path"],
        doc["input_manifest"]["base_dir"],
        doc["input_manifest"]["config_file"]["path"],
    ]
    for v in p_values:
        assert v not in text, v
    assert str(tmp_path) not in text and str(DATA) not in text and str(toml_bundle) not in text
    assert not re.search(r"(?<![\w.])/(home|tmp|Users|var|mnt|root)/", text)
    assert "timings" in text and "not printed" in text
    for v in doc["run"]["timings"]["values_s"].values():
        s = format(v, ".6g")
        assert len(s) < 5 or s not in text, s
    assert not RESERVED.search(_generated(text))
    for p in sorted((toml_bundle / "result.json").parent.rglob("*.svg")):
        t = p.read_text(encoding="utf-8")
        assert doc["run"]["output_path"] not in t and str(DATA) not in t


def test_rep3_banners_statuses_and_notes() -> None:
    doc = json.loads((GOLDEN / "result.json").read_text(encoding="utf-8"))
    arrays = {
        e["name"]: np.load(GOLDEN / "arrays" / f"{e['name']}.npy", allow_pickle=False)
        for e in doc["arrays"]
    }
    base = render_report(doc, arrays)
    assert render_report(doc, arrays) == base  # pure
    assert "INVALID RESULT" not in base and report.SHARED_LINEAGE_NOTE not in base
    bad = copy.deepcopy(doc)
    bad["valid"] = False
    bad["accounting"]["valid"] = False
    bad["accounting"]["reasons"] = ["stall"]
    text = render_report(bad, arrays)
    assert "**INVALID RESULT: transport-limit counters are nonzero" in text
    assert "- invalidating counter: stall" in text
    # every registry record: reproduced verbatim, never paraphrased; the lineage note follows
    # a topas-agreement-shared-lineage record and only that
    seen = set()
    for recs in STOPPING_EVIDENCE.values():
        for rec in recs:
            d = copy.deepcopy(doc)
            d["provenance"]["models"]["evidence_at_run"] = [dataclasses.asdict(rec)]
            t = render_report(d, arrays)
            assert [json.loads(b) for b in FENCE.findall(t)] == [dataclasses.asdict(rec)]
            assert f'"status": "{rec.status}"' in t
            shared = rec.status.startswith("topas-agreement-shared-lineage:")
            assert (report.SHARED_LINEAGE_NOTE in _generated(t)) is shared
            if shared:
                assert f"```\n\n{report.SHARED_LINEAGE_NOTE}\n" in t
            assert not RESERVED.search(_generated(t))
            seen.add(rec.status.split(":")[0])
    assert "topas-agreement-shared-lineage" in seen and "unvalidated" in seen


def test_rep3_registry_guard_and_reserved_word_in_user_text() -> None:
    for recs in STOPPING_EVIDENCE.values():
        for rec in recs:  # no registry text carries the reserved word (decision 0044)
            assert not RESERVED.search(json.dumps(dataclasses.asdict(rec)))
    doc = json.loads((GOLDEN / "result.json").read_text(encoding="utf-8"))
    doc["effective_extra"]["materials"][0]["name"] = "validated"  # user-controlled free text
    doc["statistics"]["grids"][0]["name"] = "a validated grid"
    text = render_report(doc, {})
    assert not RESERVED.search(_generated(text))
    assert "(word withheld)" in text


# ---- REP-4 ------------------------------------------------------------------------------------


def _svg_ok(text: str) -> ET.Element:
    root = ET.fromstring(text)  # parses
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert text.endswith("</svg>\n") and not text.endswith("\n\n")
    assert not re.search(r"(?i)\b(nan|inf|infinity)\b", text)
    for el in root.iter():
        for key in ("x", "y", "x1", "y1", "x2", "y2", "cx", "cy", "points", "width", "height"):
            if key in el.attrib:
                for tok in re.split(r"[ ,]", el.attrib[key]):
                    assert re.fullmatch(r"-?\d+(\.\d+)?(e[+-]?\d+)?", tok), (key, tok)
    assert "http" not in re.sub(r'xmlns="http://www.w3.org/2000/svg"', "", text)
    assert "<script" not in text and "href" not in text and "@import" not in text
    return root


def test_rep4_svgs_parse_for_real_runs_and_golden(
    toml_bundle: Path, tilted_bundle: Path, golden: Path
) -> None:
    for b in (toml_bundle, tilted_bundle, golden):
        plots = sorted((b / "plots").iterdir())
        assert len(plots) >= 3  # dose, energy, lateral (and the LET tally)
        for p in plots:
            _svg_ok(p.read_text(encoding="utf-8"))
    names = {p.name for p in (golden / "plots").iterdir()}
    assert {
        "dose__dose.svg",
        "energy__dose.svg",
        "lateral__dose.svg",
        "let__dose__let.svg",
    } <= names
    assert any(p.name.startswith("let__") for p in (tilted_bundle / "plots").iterdir())


def _golden_inputs() -> tuple[dict[str, Any], dict[str, np.ndarray[Any, Any]]]:
    doc = json.loads((GOLDEN / "result.json").read_text(encoding="utf-8"))
    arrays = {
        e["name"]: np.load(GOLDEN / "arrays" / f"{e['name']}.npy", allow_pickle=False)
        for e in doc["arrays"]
    }
    return doc, arrays


def _count(text: str, tag: str) -> int:
    return len(ET.fromstring(text).findall(f".//{{http://www.w3.org/2000/svg}}{tag}"))


def test_rep4_undefined_voxels_are_gaps() -> None:
    doc, arrays = _golden_inputs()
    plots = render_plots(doc, arrays)
    let_name = "plots/let__dose__let.svg"
    assert _count(plots[let_name], "polyline") == 1  # the first y layer is undefined: a leading gap
    mask_name, mean_name = "grid__dose__q__let__defined_mask", "grid__dose__q__let__mean"
    arrays = {k: v.copy() for k, v in arrays.items()}
    arrays[mask_name][:, 4, :] = False  # an interior undefined layer
    arrays[mean_name][:, 4, :] = np.nan
    gapped = render_plots(doc, arrays)[let_name]
    assert _count(gapped, "polyline") == 2  # the path is broken, not bridged
    _svg_ok(gapped)
    # the same for the dose band: an undefined interior dose voxel splits line and polygon
    arrays["grid__dose__defined_mask"][0, 3, 0] = False
    split = render_plots(doc, arrays)["plots/dose__dose.svg"]
    _svg_ok(split)
    assert _count(split, "polyline") == 2 and _count(split, "polygon") == 2
    # a grid without any defined voxel still yields a parseable plot saying so
    none = {k: v.copy() for k, v in arrays.items()}
    none["grid__dose__defined_mask"][:] = False
    empty = render_plots(doc, none)
    assert "no defined values" in empty["plots/dose__dose.svg"]
    for t in empty.values():
        _svg_ok(t)


def test_rep4_profile_axis_is_the_dominant_source_axis(
    tilted_bundle: Path, toml_bundle: Path, golden: Path
) -> None:
    for b, axis in ((tilted_bundle, 0), (toml_bundle, 2), (golden, 1)):
        doc = load(b).doc
        d = [abs(c) for c in doc["effective"]["unit_direction"]]
        assert doc["effective_extra"]["profile_axis"] == axis == d.index(max(d))
        text = (b / "report.md").read_text(encoding="utf-8")
        assert f"profile axis: {axis} ({'xyz'[axis]})" in text
        dose = (b / "plots" / "dose__dose.svg").read_text(encoding="utf-8")
        assert f"{'xyz'[axis]} (mm), profile axis" in dose
    # ties go to the lowest axis index
    assert results._profile_axis((1.0, 1.0, 0.0)) == 0
    assert results._profile_axis((0.0, -1.0, 1.0)) == 1
    assert results._profile_axis((0.3, -0.9, 0.9)) == 1
    assert results._profile_axis((1.0, 1.0, 1.0)) == 0


def test_rep4_writer_is_deterministic_and_escapes() -> None:
    fig = svg.Figure('a <b> & "c"', "x", "y", "note")
    fig.line("s<1>", [0.0, 1.0, 2.0, 3.0], [1.0, float("nan"), 2.0, 3.0])
    first = fig.render()
    assert first == fig.render()
    root = _svg_ok(first)
    assert "&lt;b&gt;" in first and "&amp;" in first and "s&lt;1&gt;" in first
    assert len(root.findall(".//{http://www.w3.org/2000/svg}circle")) == 1  # an isolated point
    assert len(root.findall(".//{http://www.w3.org/2000/svg}polyline")) == 1
    assert svg.num(0.1 + 0.2) == "0.3" and svg.num(-0.0) == "0"
    with pytest.raises(ValueError, match="non-finite"):
        svg.num(float("nan"))
    assert svg.nice_ticks(0.0, 1.0) == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def test_rep4_derived_files_never_print_timings(tilted_bundle: Path) -> None:
    loaded = load(tilted_bundle)
    assert isinstance(loaded.doc["run"]["timings"], dict)
    for n, b in _derived(tilted_bundle).items():
        t = b.decode("utf-8")
        for v in loaded.doc["run"]["timings"]["values_s"].values():
            s = format(v, ".6g")
            assert len(s) < 5 or s not in t, (n, s)


# ---- REP-4: the +-1 standard-error band is never clipped (REVIEW-92b2e8cf finding 1) -----------


def test_rep4_band_crosses_zero_where_the_standard_error_exceeds_the_dose() -> None:
    """The band is drawn as mean - std to mean + std; with std > mean its lower edge is negative
    and the plot's y axis extends below zero (a negative tick label appears)."""
    import xml.etree.ElementTree as ET

    from ionmc import report, svg

    fig = svg.Figure("t", "x", "y")
    fig.band([0.0, 1.0, 2.0], [-0.5, -0.2, -0.1], [1.0, 1.2, 1.1])
    assert fig._ranges()[2] < 0.0
    doc, arrays = _golden_inputs()
    dose = arrays["grid__dose__dose_gy"]
    big = dict(arrays)
    big["grid__dose__dose_std_gy"] = np.where(np.isnan(dose), np.nan, 10.0 * dose + 1e-9)
    plots = report.render_plots(doc, big)

    def y_ticks(text: str) -> list[str]:  # the y-axis tick labels are the end-anchored texts
        root = ET.fromstring(text)
        return [t.text or "" for t in root.iter("{http://www.w3.org/2000/svg}text")
                if t.get("text-anchor") == "end"]  # fmt: skip

    for name in ("plots/dose__dose.svg", "plots/lateral__dose.svg"):
        assert any(label.startswith("-") for label in y_ticks(plots[name])), name
        assert "band = +-1 standard error" in plots[name]
    # the control: the committed arrays (std < mean everywhere defined) have no negative y tick
    control = report.render_plots(doc, arrays)["plots/dose__dose.svg"]
    assert y_ticks(control) and not any(t.startswith("-") for t in y_ticks(control))
