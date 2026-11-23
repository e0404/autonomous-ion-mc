"""Row DOC-1 (V3-006 step S7): the new documentation pages are checked against the code.

CI, no data cache. Every expectation is derived from the modules; nothing here copies a constant.
Not a clinical claim.
"""

from __future__ import annotations

import dataclasses
import inspect
import re
import tomllib
from pathlib import Path

import pytest

from ionmc import cli, configio, contract, results
from ionmc import config as ionconfig

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
NEW_PAGES = (
    "results/index.md",
    "results/format.md",
    "configuration/index.md",
    "cli.md",
    "workflow.md",
)


def _page(name: str) -> str:
    return (DOCS / name).read_text(encoding="utf-8")


def _sections(text: str) -> dict[str, str]:
    """Heading line -> the text up to the next heading."""
    out: dict[str, str] = {}
    head: str | None = None
    buf: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            if head is not None:
                out[head] = "\n".join(buf)
            head, buf = line.lstrip("#").strip(), []
        elif head is not None:
            buf.append(line)
    if head is not None:
        out[head] = "\n".join(buf)
    return out


def test_doc1_results_format_matches_schema_constants() -> None:
    text = _page("results/format.md")
    for token in (results.SCHEMA, results.SCHEMA_VERSION, results.KIND):
        assert token in text
    roles = [
        *results.GRID_ARRAYS,
        *results.QUANTITY_ARRAYS,
        *results.QUANTITY_BATCH_ARRAYS,
        *results.GEOMETRY_ARRAYS,
        *results.DIAGNOSTIC_ARRAYS,
    ]
    for name in roles:
        assert f"`{name}`" in text, name
    for table in (results.GRID_ARRAYS, results.GEOMETRY_ARRAYS, results.DIAGNOSTIC_ARRAYS):
        for name, (units, _definition) in table.items():
            row = next(line for line in text.splitlines() if line.startswith(f"| `{name}` |"))
            assert f"| {units} |" in row, name
    for unit in results.UNIT_VOCABULARY:
        assert f"`{unit}`" in text, unit
    for cls, unit in results.UNITS_TABLE.items():
        assert f"| {cls} | `{unit}` |" in text, cls
    budgets = {
        n: v for n, v in vars(results).items() if n.startswith("MAX_") and isinstance(v, int)
    }
    assert {"MAX_TOTAL_BYTES", "MAX_JSON_TOKENS", "MAX_NDIM_SPECTRUM_BATCH"} <= set(budgets)
    for name, value in budgets.items():
        assert name in text, name
        assert re.search(rf"\b{value}\b", text), (name, value)
    for fixed in results._FIXED_FILES:
        assert f"`{fixed}`" in text
    for name in results.LOOKUP_PLAN_EXCLUSIONS:
        assert name in text, name


def test_doc1_budget_options_documented_with_defaults() -> None:
    fmt, cli_page = _page("results/format.md"), _page("cli.md")
    for flag, dest, _what in cli.BUDGET_OPTIONS:
        assert f"`{flag}`" in cli_page and f"`{flag}`" in fmt, flag
        const = "MAX_" + dest[len("max_") :].upper()
        assert hasattr(results, const), const
        assert f"`{const}`" in fmt


SCHEMA_SECTIONS = {
    "PHYSICS_FIELDS": "[physics]",
    "RUN_FIELDS": "[run]",
    "DIAGNOSTICS_FIELDS": "[diagnostics]",
    "SOURCE_FIELDS": "[source]",
    "PROJECTILE_FIELDS": "Inline projectile definition",
    "SCORING_FIELDS": "[[scoring]]",
    "TALLY_FIELDS": "[[tallies]]",
    "MATERIAL_FIELDS": "[[materials]]",
    "STERNHEIMER_FIELDS": "Sternheimer parameters",
    "BOX_FIELDS": "Geometry type box",
    "VOXEL_FIELDS": "Geometry type voxel",
    "OUTPUT_FIELDS": "[output]",
}


@pytest.mark.parametrize("table", sorted(SCHEMA_SECTIONS))
def test_doc1_configuration_names_every_key(table: str) -> None:
    sections = _sections(_page("configuration/index.md"))
    body = sections[SCHEMA_SECTIONS[table]]
    for key in getattr(configio, table):
        assert f"`{key}`" in body, (table, key)


def test_doc1_configuration_top_level_lookups_and_defaults() -> None:
    text = _page("configuration/index.md")
    sections = _sections(text)
    for key in configio.TOP_KEYS:
        assert f"| `{key}` |" in sections["Top level"], key
    for key in configio.REQUIRED_TOP:
        row = next(x for x in sections["Top level"].splitlines() if x.startswith(f"| `{key}` |"))
        assert "| yes |" in row, key
    lookups = sections["[[lookups]]"]
    for key in ("file", "sha256"):  # the lookup record keys of ``configio._lookups_from``
        assert f"`{key}`" in lookups
    # defaults of the dataclasses appear in their tables (numbers and strings)
    for cls, table in (
        (ionconfig.PhysicsOptions, "[physics]"),
        (ionconfig.RunOptions, "[run]"),
        (ionconfig.DiagnosticsOptions, "[diagnostics]"),
    ):
        for f in dataclasses.fields(cls):
            if f.name == "stopping_model" or f.name not in getattr(
                configio,
                {"[physics]": "PHYSICS_FIELDS", "[run]": "RUN_FIELDS"}.get(
                    table, "DIAGNOSTICS_FIELDS"
                ),
            ):
                continue
            d = f.default
            if d is dataclasses.MISSING or d is None or isinstance(d, bool):
                continue
            row = next(x for x in sections[table].splitlines() if x.startswith(f"| `{f.name}` |"))
            assert f"| {d} |" in row, (f.name, d)
    for key in ("path", "persist_batches"):
        assert f"`{key}`" in sections["[output]"]
    assert configio.SCHEMA in text and configio.SCHEMA_VERSION in text


def test_doc1_configuration_lists_every_reconciliation_row() -> None:
    body = _sections(_page("configuration/index.md"))["Requested and effective"]
    for row in contract.RECONCILIATION_TABLE:
        line = next((x for x in body.splitlines() if x.startswith(f"| `{row.id}` |")), None)
        assert line is not None, row.id
        for cls, rule in row.outcomes:
            assert f"{cls} `{rule}`" in line, (row.id, cls, rule)


def _cli_table() -> list[tuple[str, int]]:
    """(exception name, code) of every row of ``cli._exit_code``, from its source text."""
    src = inspect.getsource(cli._exit_code)
    rows: list[tuple[str, int]] = []
    for m in re.finditer(r"isinstance\(e, ([\w| ]+)\).*?,\s*(\d+)\),", src):
        for name in (n.strip() for n in m.group(1).split("|")):
            rows.append((name, int(m.group(2))))
    fallthrough = re.search(r"return (\d+)\s+#", src)
    assert fallthrough is not None
    rows.append(("any other exception", int(fallthrough.group(1))))
    return rows


def test_doc1_cli_lists_every_exit_code_with_its_exception() -> None:
    rows = _cli_table()
    assert len(rows) >= 9
    text = _page("cli.md")
    lines = [x for x in text.splitlines() if x.startswith("| ")]
    cells = [[c.strip() for c in x.split("|")[1:3]] for x in lines]
    for name, code in rows:
        label = name if name == "any other exception" else f"`{name}`"
        assert [label, str(code)] in cells, (name, code)
    for sub in ("run", "check", "capabilities", "verify", "report"):
        assert f"ionmc {sub}" in text
    for flag in ("--rederive", "--check", "--overwrite", "--out", "--inputs", "--cache-dir"):
        assert flag in text, flag


_P_PATHS = re.compile(r"(?<![\w.])(?:/home/|/Users/|/root/|/tmp/|/mnt/)\S*|[A-Za-z]:\\\S*")
_CLINICAL = re.compile(r"clinical use|treatment planning", re.IGNORECASE)


@pytest.mark.parametrize("name", NEW_PAGES)
def test_doc1_page_hygiene(name: str) -> None:
    text = _page(name)
    assert text.endswith("\n") and not text.endswith("\n\n")
    assert not re.search(r"\bvalidated\b", text)
    assert not _P_PATHS.search(text)
    for line in text.splitlines():
        if _CLINICAL.search(line):
            assert re.search(r"\bnot\b", line, re.IGNORECASE), line
    assert text.rstrip("\n").endswith("Not a clinical claim.")


def test_doc1_pages_are_in_the_navigation() -> None:
    nav = tomllib.loads((ROOT / "zensical.toml").read_text(encoding="utf-8"))["project"]["nav"]
    found: set[str] = set()

    def walk(node: object) -> None:
        if isinstance(node, str):
            found.add(node)
        elif isinstance(node, dict):
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(nav)
    for name in NEW_PAGES:
        assert name in found, name
        assert (DOCS / name).is_file()


# ---- the reproduction procedure of docs/workflow.md section 8 is executed, not only described


def test_workflow_reproduction_procedure_is_executable(tmp_path: Path) -> None:
    """``from_dict(result.json.requested)`` + the same seed, backend and precision reproduces
    every array content digest and the effective block (docs/workflow.md section 8; row OUT-3,
    small). A ``not-reconstructible`` record is refused by ``from_dict``, as the page states."""
    from ionmc.simulation import Simulation

    toml = Path(__file__).parent / "data" / "v3_006" / "proton_box.toml"
    cfg = configio.load_toml_full(toml).config
    cfg = dataclasses.replace(cfg, run=dataclasses.replace(cfg.run, n_histories=8, n_batches=4))
    first = results.save(Simulation(cfg).run(), tmp_path / "first")
    loaded = results.load(first)
    again = configio.from_dict(loaded.doc["requested"])
    second = results.save(Simulation(again).run(), tmp_path / "second")
    doc1, doc2 = loaded.doc, results.load(second).doc
    digests = lambda d: {e["name"]: e["content_sha256"] for e in d["arrays"]}  # noqa: E731
    assert digests(doc1) == digests(doc2) and digests(doc1)
    assert doc1["effective"] == doc2["effective"]
    assert doc1["requested"] == doc2["requested"]
    # the limit the page states: a not-reconstructible record is refused, never guessed
    broken = dict(loaded.doc["requested"])
    broken["geometry"] = {**broken["geometry"], "material_index": {"status": "not-reconstructible"}}
    with pytest.raises(configio.ConfigError):
        configio.from_dict(broken)
    page = (Path(__file__).parents[2] / "docs" / "workflow.md").read_text(encoding="utf-8")
    assert "configio.from_dict" in page and "not-reconstructible" in page
