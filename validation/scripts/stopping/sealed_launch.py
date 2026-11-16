"""Sealed launcher of the amendment I3-FU-1 diagnostics (stdlib only; Codex REVIEW-cf6fa397).

Usage (the only accepted launch of ``i3_band_diagnostic.py`` and ``i3fu2_evaluate.py``)::

    .venv/bin/python -I -S -B -X pycache_prefix="$(mktemp -d)" \
        validation/scripts/stopping/sealed_launch.py <target script> <args...>

``-I`` (isolated: no PYTHON* environment variables, no user site, no script directory on the
path), ``-S`` (no ``site``: no ``.pth`` files, no ``sitecustomize`` / ``usercustomize``), ``-B``
(no bytecode written) and a fresh empty ``pycache_prefix`` (no bytecode read for any module
imported after start-up). The launcher

(a) checks ``sys.flags.isolated == 1``, ``sys.flags.no_site == 1``,
    ``sys.flags.dont_write_bytecode == 1`` and that ``sys.pycache_prefix`` is an existing
    directory holding no file;
(b) checks the import machinery by identity: ``sys.meta_path`` is exactly
    [``importlib.machinery.BuiltinImporter``, ``FrozenImporter``, ``PathFinder``] (the class
    objects themselves); ``sys.path_hooks`` is exactly [``zipimport.zipimporter``, the
    ``FileFinder`` path hook] (the second by module ``_frozen_importlib_external``, qualified
    name ``FileFinder.path_hook.<locals>.path_hook_for_FileFinder`` and code file
    ``<frozen importlib._bootstrap_external>``); ``sys.path_importer_cache`` holds only
    ``FileFinder`` / ``zipimporter`` entries or ``None``;
(c) sets ``sys.path`` to exactly [<repo>/src, <repo>/validation/scripts/stopping, <venv
    site-packages>] followed by the stdlib entries already present. Under ``-S`` the interpreter
    does not process ``pyvenv.cfg`` (``sys.prefix`` is the base installation), so the venv's
    site-packages is resolved from the launcher's own interpreter path: the directory two
    levels above ``sys.executable`` must hold ``pyvenv.cfg``, and site-packages is
    ``<venv>/lib/python<major>.<minor>/site-packages``;
(d) records the sha256 of every stdlib module file loaded so far (modules whose ``__file__``
    lies under ``sysconfig.get_paths()["stdlib"]``), its own sha256 and the above, and passes
    the record to the target as ``builtins.__ionmc_sealed_launch__`` (a dict; set immediately
    before the target runs, read by ``env_seal.capture_launch``);
(e) runs the target as ``__main__`` with ``runpy.run_path``.

Any failed check: a message on stderr and exit status 2; the target is not run.
"""

from __future__ import annotations

import builtins
import hashlib
import importlib.machinery
import os
import runpy
import sys
import sysconfig
import zipimport
from pathlib import Path
from typing import Any

RECORD_NAME = "__ionmc_sealed_launch__"
ORIGIN = "validation/scripts/stopping/sealed_launch.py"
FILEFINDER_HOOK = (
    "_frozen_importlib_external",
    "FileFinder.path_hook.<locals>.path_hook_for_FileFinder",
    "<frozen importlib._bootstrap_external>",
)
REPO = Path(__file__).resolve().parents[3]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def flag_problems() -> list[str]:
    """(a): isolation flags and the fresh empty pycache_prefix."""
    out = []
    for name in ("isolated", "no_site", "dont_write_bytecode"):
        if getattr(sys.flags, name) != 1:
            out.append(f"sys.flags.{name} is not 1 (launch with python -I -S -B)")
    prefix = getattr(sys, "pycache_prefix", None)
    if not prefix:
        out.append("sys.pycache_prefix is not set (-X pycache_prefix=<fresh empty directory>)")
    elif not Path(prefix).is_dir():
        out.append("the pycache_prefix directory does not exist")
    elif any(p.is_file() for p in Path(prefix).rglob("*")):
        out.append("the pycache_prefix directory is not empty")
    return out


def machinery_problems(
    meta_path: list[Any] | None = None, path_hooks: list[Any] | None = None
) -> list[str]:
    """(b): ``sys.meta_path`` / ``sys.path_hooks`` are the standard ones, by identity."""
    meta = list(sys.meta_path if meta_path is None else meta_path)
    hooks = list(sys.path_hooks if path_hooks is None else path_hooks)
    want = [
        importlib.machinery.BuiltinImporter,
        importlib.machinery.FrozenImporter,
        importlib.machinery.PathFinder,
    ]
    out = []
    if len(meta) != len(want) or any(m is not w for m, w in zip(meta, want, strict=False)):
        out.append(f"sys.meta_path is not the standard finder list: {meta!r}")
    ok_hooks = len(hooks) == 2 and hooks[0] is zipimport.zipimporter
    if ok_hooks:
        hook = hooks[1]
        code = getattr(hook, "__code__", None)
        ok_hooks = (
            getattr(hook, "__module__", None),
            getattr(hook, "__qualname__", None),
            getattr(code, "co_filename", None),
        ) == FILEFINDER_HOOK
    if not ok_hooks:
        out.append(f"sys.path_hooks is not the standard hook list: {hooks!r}")
    for key, finder in sys.path_importer_cache.items():
        if finder is not None and not isinstance(
            finder, importlib.machinery.FileFinder | zipimport.zipimporter
        ):
            out.append(f"sys.path_importer_cache[{key!r}] holds a foreign finder {finder!r}")
    return out


def venv_site_packages() -> Path:
    """The venv's site-packages of the launcher's own interpreter (see (c))."""
    venv = Path(sys.executable).parent.parent
    if not (venv / "pyvenv.cfg").is_file():
        raise RuntimeError(f"{sys.executable} is not a venv interpreter (no pyvenv.cfg)")
    major, minor = sys.version_info[:2]
    site = venv / "lib" / f"python{major}.{minor}" / "site-packages"
    if not site.is_dir():
        raise RuntimeError(f"no site-packages directory {site}")
    return site.resolve()


def stdlib_dir() -> Path:
    return Path(sysconfig.get_paths()["stdlib"]).resolve()


def loaded_stdlib_hashes() -> dict[str, str]:
    """(d): {path relative to the stdlib dir: sha256} of every loaded stdlib module file."""
    base = stdlib_dir()
    out = {}
    for mod in list(sys.modules.values()):
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        path = Path(f).resolve()
        if path.is_relative_to(base) and path.is_file():
            out[path.relative_to(base).as_posix()] = _sha256(path)
    return dict(sorted(out.items()))


def launch_record(target: Path) -> dict[str, Any]:
    """Run checks (a)-(d); raise RuntimeError on any problem; return the record."""
    problems = flag_problems() + machinery_problems()
    if problems:
        raise RuntimeError("; ".join(problems))
    base = stdlib_dir()
    stdlib_entries = [p for p in sys.path if Path(p).resolve().is_relative_to(base.parent)]
    sys.path[:] = [
        str(REPO / "src"),
        str(REPO / "validation" / "scripts" / "stopping"),
        str(venv_site_packages()),
        *stdlib_entries,
    ]
    return {
        "origin": ORIGIN,
        "launcher_sha256": _sha256(Path(__file__).resolve()),
        "target": str(target),
        "python_executable": str(Path(sys.executable).resolve()),
        "flags": {
            "isolated": sys.flags.isolated,
            "no_site": sys.flags.no_site,
            "dont_write_bytecode": sys.flags.dont_write_bytecode,
        },
        "pycache_prefix": sys.pycache_prefix,
        "pycache_prefix_empty": True,
        "machinery_ok": True,
        "sys_path": list(sys.path),
        "sys_prefix": sys.prefix,
        "sys_base_prefix": sys.base_prefix,
        "stdlib_dir": str(base),
        "stdlib_modules_at_launch": loaded_stdlib_hashes(),
        "pid": os.getpid(),
    }


def main(argv: list[str]) -> int:
    if not argv:
        print("usage: sealed_launch.py <target script> <args...>", file=sys.stderr)
        return 2
    target = Path(argv[0]).resolve()
    try:
        record = launch_record(target)
    except RuntimeError as exc:
        print(f"sealed_launch: refused: {exc}", file=sys.stderr)
        return 2
    setattr(builtins, RECORD_NAME, record)
    sys.argv = [str(target), *argv[1:]]
    runpy.run_path(str(target), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
