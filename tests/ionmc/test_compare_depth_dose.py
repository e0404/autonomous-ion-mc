"""Fail-closed behaviour of the reference depth-dose comparison script."""

from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/reference/compare_depth_dose.py"


def _load():
    spec = importlib.util.spec_from_file_location("compare_depth_dose", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_batches_option_removed_for_any_metadata(tmp_path: Path) -> None:
    mod = _load()
    cases = {
        "duplicate_seeds": {"n_batches": 2, "seeds": [1, 1]},
        "mismatched_counts": {"n_batches": 5, "seeds": [1, 2]},
        "unrelated": {"foo": "bar"},
    }
    for name, meta in cases.items():
        f = tmp_path / f"{name}.json"
        f.write_text(json.dumps(meta))
        with pytest.raises(SystemExit) as exc:
            mod.main(
                ["--runs", str(tmp_path), "--output", str(tmp_path / "o.json"), "--batches", str(f)]
            )
        assert exc.value.code == 2
        assert not (tmp_path / "o.json").exists()


def test_status_is_exploratory() -> None:
    assert _load().EXPLORATORY.startswith("exploratory")


def test_code_sha_must_equal_head() -> None:
    mod = _load()
    head, dirty = mod.code_state(None)
    assert len(head) == 40 and isinstance(dirty, bool)
    assert mod.code_state(head) == (head, dirty)
    with pytest.raises(SystemExit) as exc:
        mod.code_state("0" * 40)
    assert exc.value.code != 0


def test_head_matches_git() -> None:
    head = subprocess.run(
        ["git", "-C", str(SCRIPT.parent), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    assert _load().code_state(None)[0] == head
