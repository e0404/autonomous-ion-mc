"""Environment seal shared by the amendment I3-FU-1 diagnostics (step 1 and step 2).

Codex REVIEW-c4f3fdd1 finding 3: the lock file pins versions, not bytes. This module seals the
numerical runtime of the diagnostics by content:

* the installed numpy distribution is verified against its own wheel manifest
  (``numpy-<version>.dist-info/RECORD``: csv rows ``path, sha256=<urlsafe base64>, size``):
  every listed file with a hash is hashed and compared (size too); the RECORD file itself and
  entries without a hash are skipped; a missing file or a mismatch fails closed;
* the running environment (:func:`current_environment`): numpy version, sha256 of the RECORD
  file, sha256 of numpy's ``_core/_multiarray_umath*`` extension, sha256 of the resolved Python
  executable (the real binary, not the venv symlink), ``sys.version``, ``platform.platform()``,
  ``sys.prefix``, ``sys.base_prefix`` and the sha256 of the stdlib directory listing;
* the stdlib (Codex REVIEW-cf6fa397): the seal holds ``stdlib_modules`` = {path relative to the
  stdlib directory: sha256} of every stdlib module file (``.py`` and compiled extensions; the
  directories ``site-packages``, ``test`` and ``__pycache__`` excluded), a superset of any loaded
  set so that modules imported lazily during a computation are covered; every loaded stdlib
  module must be in the seal with a matching hash (:func:`stdlib_modules_check`), at the start
  and again after the computation;
* the committed environment seal (``validation/results/stopping/v3-008-i3fu1/
  environment-seal.json``, :data:`SEAL_FIELDS` and ``stdlib_modules``), written once from the
  reviewed environment by ``i3_band_diagnostic.py --write-environment-seal`` launched through
  ``sealed_launch.py``; both diagnostics require every field to equal the running environment
  (:func:`verify_environment`);
* launch control (Codex REVIEW-8081e0e7 finding 3, REVIEW-cf6fa397): both diagnostics must be
  launched through the stdlib-only launcher ``sealed_launch.py`` as ``python -I -S -B -X
  pycache_prefix=<fresh empty directory> sealed_launch.py <script> ...``: no PYTHON*
  environment variables, no ``site``, no ``.pth`` files, no ``sitecustomize`` /
  ``usercustomize``, no bytecode read or written; the launcher checks the flags, the empty
  prefix and the import machinery by identity (standard ``sys.meta_path`` finders and
  ``sys.path_hooks``), sets an explicit ``sys.path`` (repo ``src``, the scripts directory, the
  venv's site-packages, the stdlib entries) and records the stdlib modules loaded so far.
  :func:`capture_launch` runs as the first statement of each diagnostic (only stdlib imported
  before it) and :func:`check_launch` ENFORCES before any computation: the launcher record is
  present (origin ``sealed_launch.py``), the flags are -I -S -B, the prefix held no file, the
  machinery check passes again at the capture, and ``sys.path`` equals the launcher's list;
  without the launcher record the diagnostics refuse to run. Each diagnostic also requires the
  launcher's recorded sha256 to equal its committed (sealed) bytes;
* runtime controls (Codex REVIEW-183057ce finding 3): every loaded numpy module resolves under the
  verified distribution root and its file is listed with a hash in RECORD (hence verified;
  :func:`numpy_modules_check`), and, as a second line of defence only, no ``__cached__``
  bytecode file exists for any sealed or numpy module (:func:`bytecode_check`; the absence of a
  cache file at check time cannot prove that none was executed earlier, the prefix mechanism
  can). Both diagnostics repeat these checks and the stdlib check after the computation, just
  before writing, so lazily imported modules are covered.

What is sealed: the interpreter binary, the stdlib module files, the numpy distribution, the
ionmc and script closure (by the diagnostics), the launcher. Declared limitation (what remains
outside): the operating-system kernel, libc and the other system libraries the interpreter and
its extensions link against (the C runtime of the interpreter), bound only through
``platform.platform()``; numpy's own bundled libraries are covered where its RECORD lists them.
The interpreter's built-in and frozen modules are part of the sealed interpreter binary.
"""

from __future__ import annotations

import base64
import builtins
import csv
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
import platform
import sys
import sysconfig
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

SEAL_REL = "validation/results/stopping/v3-008-i3fu1/environment-seal.json"
SEAL_SCHEMA = "ionmc-i3fu1-environment-seal-2"
SEAL_FIELDS = (
    "numpy_version",
    "numpy_record_sha256",
    "numpy_multiarray_sha256",
    "python_executable_sha256",
    "sys_version",
    "platform",
    "sys_prefix",
    "sys_base_prefix",
    "stdlib_dir_sha256_of_listing",
)
#: seal key of the {path relative to the stdlib dir: sha256} map of the stdlib module files
STDLIB_KEY = "stdlib_modules"
#: stdlib subdirectories outside the seal (third-party packages, the test suite, caches)
STDLIB_EXCLUDED = ("site-packages", "dist-packages", "test", "__pycache__")
LAUNCHER_REL = "validation/scripts/stopping/sealed_launch.py"


class EnvironmentSealError(RuntimeError):
    """The running environment does not match its wheel manifest or the committed seal."""


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def numpy_dist_info() -> Path:
    """The ``numpy-<version>.dist-info`` directory of the imported numpy."""
    import numpy

    dist = importlib.metadata.distribution("numpy")
    if dist.version != numpy.__version__:
        raise EnvironmentSealError(
            f"numpy distribution {dist.version} is not the imported numpy {numpy.__version__}"
        )
    path = Path(str(dist.locate_file(""))) / f"numpy-{dist.version}.dist-info"
    if not (path / "RECORD").is_file():
        raise EnvironmentSealError(f"no RECORD in {path}")
    return path


def verify_record(dist_info: Path) -> dict[str, Any]:
    """Hash every file listed with a hash in ``<dist_info>/RECORD`` against its row; returns
    ``{"files_verified", "record_sha256"}``. Missing file, unknown algorithm, size or hash
    mismatch: :class:`EnvironmentSealError`."""
    record = dist_info / "RECORD"
    base = dist_info.parent
    text = record.read_text(encoding="utf-8")
    n = 0
    for row in csv.reader(io.StringIO(text)):
        if not row:
            continue
        if len(row) != 3:
            raise EnvironmentSealError(f"RECORD row {row!r} does not have 3 fields")
        rel, digest, size = row
        if not digest or (base / rel).resolve() == record.resolve():
            continue
        algo, _, b64 = digest.partition("=")
        if algo != "sha256" or not b64:
            raise EnvironmentSealError(f"RECORD row {rel}: unsupported hash {digest!r}")
        path = base / rel
        if not path.is_file():
            raise EnvironmentSealError(f"RECORD lists {rel}, which is missing")
        blob = path.read_bytes()
        got = base64.urlsafe_b64encode(hashlib.sha256(blob).digest()).rstrip(b"=").decode()
        if got != b64:
            raise EnvironmentSealError(f"RECORD hash mismatch for {rel}")
        if size and int(size) != len(blob):
            raise EnvironmentSealError(f"RECORD size mismatch for {rel}")
        n += 1
    if n == 0:
        raise EnvironmentSealError(f"RECORD of {dist_info} lists no hashed file")
    return {"files_verified": n, "record_sha256": _sha256_file(record)}


def multiarray_extension() -> Path:
    """numpy's ``_core/_multiarray_umath*`` compiled extension (exactly one)."""
    import numpy

    found = sorted(Path(numpy.__file__).resolve().parent.glob("_core/_multiarray_umath*.so"))
    found += sorted(Path(numpy.__file__).resolve().parent.glob("_core/_multiarray_umath*.pyd"))
    if len(found) != 1:
        raise EnvironmentSealError(f"expected one _multiarray_umath extension, found {found}")
    return found[0]


def stdlib_dir() -> Path:
    """The stdlib directory of the running interpreter (``sysconfig``)."""
    return Path(sysconfig.get_paths()["stdlib"]).resolve()


def _stdlib_walk() -> list[Path]:
    base = stdlib_dir()
    out = []
    for root, dirs, files in os.walk(base):
        rel_root = Path(root).relative_to(base)
        dirs[:] = sorted(
            d
            for d in dirs
            if d != "__pycache__" and not (rel_root == Path(".") and d in STDLIB_EXCLUDED)
        )
        out += [Path(root) / f for f in sorted(files)]
    return out


def stdlib_listing_sha256() -> str:
    """sha256 of the sorted listing (one relative path per line) of every file of the stdlib
    directory outside :data:`STDLIB_EXCLUDED`."""
    base = stdlib_dir()
    rels = sorted(p.relative_to(base).as_posix() for p in _stdlib_walk())
    return hashlib.sha256("".join(r + "\n" for r in rels).encode()).hexdigest()


def stdlib_files() -> dict[str, str]:
    """{relative path: sha256} of every stdlib module file (``.py`` and compiled extensions)
    outside :data:`STDLIB_EXCLUDED`: the seal holds the whole importable stdlib, a superset of
    any loaded set, so that modules imported lazily during a computation are covered."""
    base = stdlib_dir()
    return {
        p.relative_to(base).as_posix(): _sha256_file(p)
        for p in _stdlib_walk()
        if p.suffix in (".py", ".so", ".pyd")
    }


def loaded_stdlib_modules(modules: Mapping[str, ModuleType] | None = None) -> dict[str, Path]:
    """{relative path: file} of every loaded module whose file lies in the stdlib directory
    (outside :data:`STDLIB_EXCLUDED`)."""
    base = stdlib_dir()
    mods = sys.modules if modules is None else modules
    out = {}
    for mod in list(mods.values()):
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        path = Path(f).resolve()
        if not path.is_relative_to(base):
            continue
        rel = path.relative_to(base)
        if rel.parts and rel.parts[0] in STDLIB_EXCLUDED:
            continue
        out[rel.as_posix()] = path
    return out


def stdlib_modules_check(
    seal_map: Mapping[str, str], modules: Mapping[str, ModuleType] | None = None
) -> dict[str, Any]:
    """Every loaded stdlib module file is in the seal with a matching sha256; else
    :class:`EnvironmentSealError`."""
    loaded = loaded_stdlib_modules(modules)
    for rel, path in sorted(loaded.items()):
        if rel not in seal_map:
            raise EnvironmentSealError(f"stdlib module {rel} is not in the environment seal")
        if _sha256_file(path) != seal_map[rel]:
            raise EnvironmentSealError(f"stdlib module {rel} differs from the environment seal")
    return {"stdlib_modules_verified": len(loaded)}


def current_environment(full: bool = False) -> dict[str, Any]:
    """The seal fields of the running environment (the numpy RECORD is verified first);
    ``full``: also the stdlib module map (for writing a seal)."""
    import numpy

    record = verify_record(numpy_dist_info())
    executable = Path(sys.executable).resolve()
    extra = {STDLIB_KEY: stdlib_files()} if full else {}
    return {
        **extra,
        "numpy_version": numpy.__version__,
        "numpy_record_sha256": record["record_sha256"],
        "numpy_multiarray_sha256": _sha256_file(multiarray_extension()),
        "python_executable_sha256": _sha256_file(executable),
        "sys_version": sys.version,
        "platform": platform.platform(),
        "sys_prefix": sys.prefix,
        "sys_base_prefix": sys.base_prefix,
        "stdlib_dir_sha256_of_listing": stdlib_listing_sha256(),
        "_numpy_record_files_verified": record["files_verified"],
        "_python_executable": str(executable),
    }


def seal_document(env: dict[str, Any]) -> dict[str, Any]:
    """The committed seal document of an environment (:data:`SEAL_FIELDS` and the stdlib module
    map; computed here if ``env`` lacks it)."""
    stdlib = env.get(STDLIB_KEY) or stdlib_files()
    return {"schema": SEAL_SCHEMA, **{k: env[k] for k in SEAL_FIELDS}, STDLIB_KEY: stdlib}


def seal_text(env: dict[str, Any]) -> str:
    return json.dumps(seal_document(env), indent=1, sort_keys=True) + "\n"


def parse_seal(raw: bytes) -> dict[str, Any]:
    try:
        doc = json.loads(raw.decode("utf-8"))
    except ValueError as exc:
        raise EnvironmentSealError(f"environment seal is not JSON: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("schema") != SEAL_SCHEMA:
        raise EnvironmentSealError(f"environment seal schema is not {SEAL_SCHEMA}")
    if sorted(doc) != sorted(("schema", *SEAL_FIELDS, STDLIB_KEY)):
        raise EnvironmentSealError(
            f"environment seal must hold exactly {(*SEAL_FIELDS, STDLIB_KEY)}"
        )
    for key in SEAL_FIELDS:
        if not isinstance(doc[key], str) or not doc[key]:
            raise EnvironmentSealError(f"environment seal field {key} is empty")
    stdlib = doc[STDLIB_KEY]
    if (
        not isinstance(stdlib, dict)
        or not stdlib
        or not all(isinstance(k, str) and isinstance(v, str) for k, v in stdlib.items())
    ):
        raise EnvironmentSealError("environment seal: stdlib_modules is not a {path: sha256} map")
    return doc


def verify_environment(raw: bytes, env: dict[str, Any] | None = None) -> dict[str, Any]:
    """Every seal field of the running environment (``env``, default
    :func:`current_environment`) equals the committed seal ``raw``; returns the provenance
    record ``{**fields, "numpy_record_files_verified", "seal_sha256", "match": True}``."""
    seal = parse_seal(raw)
    env = current_environment() if env is None else env
    bad = [k for k in SEAL_FIELDS if env.get(k) != seal[k]]
    if bad:
        raise EnvironmentSealError(f"running environment differs from the seal in {bad}")
    stdlib = stdlib_modules_check(seal[STDLIB_KEY])
    return {
        **{k: seal[k] for k in SEAL_FIELDS},
        "stdlib_modules_in_seal": len(seal[STDLIB_KEY]),
        **stdlib,
        "numpy_record_files_verified": env.get("_numpy_record_files_verified"),
        "seal_sha256": hashlib.sha256(raw).hexdigest(),
        "match": True,
    }


def _hashed_paths(dist_info: Path) -> set[Path]:
    """Resolved paths of the RECORD rows that carry a hash."""
    base = dist_info.parent
    out = set()
    for row in csv.reader(io.StringIO((dist_info / "RECORD").read_text(encoding="utf-8"))):
        if len(row) == 3 and row[1]:
            out.add((base / row[0]).resolve())
    return out


def numpy_modules_check(
    dist_info: Path | None = None, modules: Mapping[str, ModuleType] | None = None
) -> dict[str, Any]:
    """Every loaded numpy module with a ``__file__`` resolves under the distribution root
    (``dist_info.parent``) and is listed with a hash in RECORD (``verify_record`` hashes it):
    else :class:`EnvironmentSealError`. ``modules`` defaults to ``sys.modules``."""
    dist_info = numpy_dist_info() if dist_info is None else dist_info
    root = dist_info.parent.resolve()
    hashed = _hashed_paths(dist_info)
    mods = sys.modules if modules is None else modules
    n = 0
    for name, mod in sorted(mods.items()):
        if name != "numpy" and not name.startswith("numpy."):
            continue
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        path = Path(f).resolve()
        if not path.is_relative_to(root):
            raise EnvironmentSealError(f"numpy module {name} is loaded from outside {root}")
        if path not in hashed:
            raise EnvironmentSealError(f"numpy module {name} ({path.name}) is not hashed in RECORD")
        n += 1
    if n == 0:
        raise EnvironmentSealError("no loaded numpy module found")
    return {"numpy_modules_verified": n}


def numpy_modules(modules: Mapping[str, ModuleType] | None = None) -> dict[str, ModuleType]:
    mods = sys.modules if modules is None else modules
    return {k: v for k, v in mods.items() if k == "numpy" or k.startswith("numpy.")}


def bytecode_check(modules: Mapping[str, ModuleType]) -> dict[str, Any]:
    """``sys.dont_write_bytecode`` must be true and no module of ``modules`` may have its
    ``__cached__`` bytecode file on disk (the verified sources are compiled in memory):
    else :class:`EnvironmentSealError`."""
    if not sys.dont_write_bytecode:
        raise EnvironmentSealError(
            "sys.dont_write_bytecode is false (set PYTHONDONTWRITEBYTECODE=1)"
        )
    stale = sorted(
        name
        for name, mod in modules.items()
        if getattr(mod, "__cached__", None) and Path(str(mod.__cached__)).exists()
    )
    if stale:
        raise EnvironmentSealError(
            f"bytecode cache files exist for {stale[:5]}{' ...' if len(stale) > 5 else ''} "
            "(remove the __pycache__ directories before the run)"
        )
    return {"dont_write_bytecode": True, "modules_without_bytecode_cache": len(modules)}


def _launcher() -> ModuleType:
    """``sealed_launch.py`` of this directory (loaded by path; its import has no side effect)."""
    name = "_ionmc_sealed_launch_checks"
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).resolve().parent / "sealed_launch.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def capture_launch() -> dict[str, Any]:
    """The launch state, captured as the first statement of a diagnostic (see the module
    docstring): the flags, ``sys.pycache_prefix`` and whether it holds no file, the import
    machinery check (by identity, :mod:`sealed_launch`), ``sys.path`` and the launcher record
    (``builtins.__ionmc_sealed_launch__``). Never raises."""
    prefix = getattr(sys, "pycache_prefix", None)
    path = Path(prefix) if prefix else None
    exists = bool(path is not None and path.is_dir())
    files = sum(1 for p in path.rglob("*") if p.is_file()) if exists and path is not None else None
    launcher = _launcher()
    record = getattr(builtins, launcher.RECORD_NAME, None)
    return {
        "dont_write_bytecode_flag": int(sys.flags.dont_write_bytecode),
        "pycache_prefix": prefix,
        "pycache_prefix_exists": exists,
        "pycache_prefix_files_at_start": files,
        "isolated_flag": int(sys.flags.isolated),
        "no_site_flag": int(sys.flags.no_site),
        "meta_path": [
            type(f).__name__ if not isinstance(f, type) else f.__name__ for f in sys.meta_path
        ],
        "machinery_problems": launcher.machinery_problems(),
        "sys_path_at_capture": list(sys.path),
        "launcher_record": dict(record) if isinstance(record, dict) else None,
    }


def check_launch(launch: Mapping[str, Any]) -> dict[str, Any]:
    """Enforce the captured launch state (Codex REVIEW-cf6fa397): launched through
    ``sealed_launch.py`` (its record present, origin and launcher sha256 recorded), ``-I -S -B``,
    ``pycache_prefix`` an existing directory that held no file at the capture, the standard import
    machinery by identity, and ``sys.path`` equal to the launcher's explicit list; else
    :class:`EnvironmentSealError`. Returns the provenance record (the launcher record without its
    stdlib hash map, which is summarised by its count)."""
    record = launch.get("launcher_record")
    if not isinstance(record, dict) or record.get("origin") != LAUNCHER_REL:
        raise EnvironmentSealError("not launched through sealed_launch.py (no launcher record)")
    if launch.get("dont_write_bytecode_flag") != 1:
        raise EnvironmentSealError("not launched with python -B (sys.flags.dont_write_bytecode)")
    if launch.get("isolated_flag") != 1 or launch.get("no_site_flag") != 1:
        raise EnvironmentSealError("not launched with python -I -S (isolated, no site)")
    if not launch.get("pycache_prefix"):
        raise EnvironmentSealError("not launched with -X pycache_prefix=<fresh empty directory>")
    if not launch.get("pycache_prefix_exists"):
        raise EnvironmentSealError("the pycache_prefix directory does not exist")
    if launch.get("pycache_prefix_files_at_start") != 0:
        raise EnvironmentSealError("the pycache_prefix directory was not empty at launch")
    if launch.get("machinery_problems"):
        raise EnvironmentSealError(f"import machinery: {launch['machinery_problems']}")
    if launch.get("sys_path_at_capture") != record.get("sys_path"):
        raise EnvironmentSealError("sys.path at the capture is not the launcher's explicit list")
    summary = {k: v for k, v in record.items() if k != "stdlib_modules_at_launch"}
    summary["stdlib_modules_at_launch"] = len(record.get("stdlib_modules_at_launch") or {})
    return {
        **{k: v for k, v in launch.items() if k != "launcher_record"},
        "launcher_record": summary,
        "launch_ok": True,
    }
