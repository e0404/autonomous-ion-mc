"""V3-005D Amendment 21 (f): the committed M2' record scripts and the projected table contain no
hard-coded measured values (the test itself holds none). Pure file reads relative to the
repository root; skipped when the record directory is absent."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "validation" / "plans" / "v3-005d-m2prime"
SCRIPTS = REC / "scripts"
PROJECTED = REC / "M2PRIME-projected.md"

pytestmark = pytest.mark.skipif(not REC.is_dir(), reason="M2' record directory absent")


def _files() -> list[Path]:
    return sorted(SCRIPTS.glob("*.py")) + [PROJECTED]


def test_scripts_read_measured_values_from_the_loader() -> None:
    for name in ("m2prime.py", "m2prime2.py"):
        text = (SCRIPTS / name).read_text(encoding="utf-8")
        assert re.search(r"^from exfor_tot import .*\bmeasured_sigma_tot\b", text, re.M), name
        # no `TOT = {...}` dict literal with numeric values
        for m in re.finditer(r"^\s*TOT\s*=\s*\{(.*?)\}", text, re.M | re.S):
            assert not re.search(r"\d", m.group(1)), f"{name}: hard-coded TOT dict"
        assert not re.search(r"^\s*TOT\s*=\s*\{", text, re.M), f"{name}: TOT dict literal"


def test_no_signed_percentages_or_numeric_result_tokens() -> None:
    pct = re.compile(r"[+-]\d+\.\d+ ?%")
    tok = re.compile(r"\b(?:sTot|sR|sEl)\b\s*=\s*[-+]?\d+\.\d+(?![\d.]|\s*[*/(\w])")
    for f in _files():
        assert f.is_file(), f
        text = f.read_text(encoding="utf-8")
        assert not pct.search(text), f"{f.name}: signed percentage"
        assert not tok.search(text), f"{f.name}: sTot/sR/sEl assigned a number"


def test_projected_table_withholds_the_m2_reference_twelve_times() -> None:
    text = PROJECTED.read_text(encoding="utf-8")
    assert text.count("withheld (M2 reference)") == 12
