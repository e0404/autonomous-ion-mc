"""Environment seal shared by the amendment I3-FU-1 diagnostics (step 1 and step 2).

Codex REVIEW-c4f3fdd1 finding 3: the lock file pins versions, not bytes. This module seals the
numerical runtime of the diagnostics by content:

* the installed numpy distribution is verified against its own wheel manifest
  (``numpy-<version>.dist-info/RECORD``: csv rows ``path, sha256=<urlsafe base64>, size``):
  every listed file with a hash is hashed and compared (size too); the RECORD file itself and
  entries without a hash are skipped; a missing file or a mismatch fails closed;
* the running environment (:func:`current_environment`): numpy version, sha256 of the RECORD
  file, sha256 of numpy's ``_core/_multiarray_umath*`` extension, sha256 of the resolved Python
  executable (the real binary, not the venv symlink), ``sys.version`` and
  ``platform.platform()``;
* the committed environment seal (``validation/results/stopping/v3-008-i3fu1/
  environment-seal.json``, :data:`SEAL_FIELDS`), written once from the reviewed environment by
  ``i3_band_diagnostic.py --write-environment-seal``; both diagnostics require every field of the
  running environment to equal it (:func:`verify_environment`);
* launch control (Codex REVIEW-8081e0e7 finding 3): both diagnostics are launched as
  ``python -B -X pycache_prefix=<fresh empty directory> <script> ...``.
  :func:`capture_launch` runs as the first statement of each script (only stdlib imported
  before it, in particular neither numpy nor ionmc) and records ``sys.flags.dont_write_bytecode``,
  ``sys.pycache_prefix`` and whether that directory exists and holds no file;
  :func:`check_launch` enforces it before any computation. With an empty prefix no bytecode can
  have been read for any module imported after the capture: the interpreter compiled the
  verified sources;
* runtime controls (Codex REVIEW-183057ce finding 3): every loaded numpy module resolves under the
  verified distribution root and its file is listed with a hash in RECORD (hence verified;
  :func:`numpy_modules_check`), and, as a second line of defence only, no ``__cached__``
  bytecode file exists for any sealed or numpy module (:func:`bytecode_check`; the absence of a
  cache file at check time cannot prove that none was executed earlier, the prefix mechanism
  can). Both diagnostics repeat these checks after the computation, just before writing, so
  lazily imported modules are covered.

Declared limitation: shared libraries outside the numpy distribution (libc, the system
libraries the interpreter links against) and the operating system are not sealed beyond
``platform.platform()``; numpy's own bundled libraries are covered where its RECORD lists them.
A consistent rewrite of a numpy file together with its RECORD row is detected through the
sealed RECORD sha256.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import importlib.metadata
import io
import json
import platform
import sys
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

SEAL_REL = "validation/results/stopping/v3-008-i3fu1/environment-seal.json"
SEAL_SCHEMA = "ionmc-i3fu1-environment-seal-1"
SEAL_FIELDS = (
    "numpy_version",
    "numpy_record_sha256",
    "numpy_multiarray_sha256",
    "python_executable_sha256",
    "sys_version",
    "platform",
)


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


def current_environment() -> dict[str, Any]:
    """The seal fields of the running environment (the numpy RECORD is verified first)."""
    import numpy

    record = verify_record(numpy_dist_info())
    executable = Path(sys.executable).resolve()
    return {
        "numpy_version": numpy.__version__,
        "numpy_record_sha256": record["record_sha256"],
        "numpy_multiarray_sha256": _sha256_file(multiarray_extension()),
        "python_executable_sha256": _sha256_file(executable),
        "sys_version": sys.version,
        "platform": platform.platform(),
        "_numpy_record_files_verified": record["files_verified"],
        "_python_executable": str(executable),
    }


def seal_document(env: dict[str, Any]) -> dict[str, Any]:
    """The committed seal document of an environment (:data:`SEAL_FIELDS` only)."""
    return {"schema": SEAL_SCHEMA, **{k: env[k] for k in SEAL_FIELDS}}


def seal_text(env: dict[str, Any]) -> str:
    return json.dumps(seal_document(env), indent=1, sort_keys=True) + "\n"


def parse_seal(raw: bytes) -> dict[str, Any]:
    try:
        doc = json.loads(raw.decode("utf-8"))
    except ValueError as exc:
        raise EnvironmentSealError(f"environment seal is not JSON: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("schema") != SEAL_SCHEMA:
        raise EnvironmentSealError(f"environment seal schema is not {SEAL_SCHEMA}")
    if sorted(doc) != sorted(("schema", *SEAL_FIELDS)):
        raise EnvironmentSealError(f"environment seal must hold exactly {SEAL_FIELDS}")
    for key in SEAL_FIELDS:
        if not isinstance(doc[key], str) or not doc[key]:
            raise EnvironmentSealError(f"environment seal field {key} is empty")
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
    return {
        **{k: seal[k] for k in SEAL_FIELDS},
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


def capture_launch() -> dict[str, Any]:
    """The launch state, captured as the first statement of a diagnostic (see the module
    docstring): ``-B`` (``sys.flags.dont_write_bytecode``), ``sys.pycache_prefix`` and whether
    the prefix directory exists and holds no file at this moment. Never raises."""
    prefix = getattr(sys, "pycache_prefix", None)
    path = Path(prefix) if prefix else None
    exists = bool(path is not None and path.is_dir())
    files = sum(1 for p in path.rglob("*") if p.is_file()) if exists and path is not None else None
    return {
        "dont_write_bytecode_flag": int(sys.flags.dont_write_bytecode),
        "pycache_prefix": prefix,
        "pycache_prefix_exists": exists,
        "pycache_prefix_files_at_start": files,
    }


def check_launch(launch: Mapping[str, Any]) -> dict[str, Any]:
    """Enforce the captured launch state: ``-B`` given, ``pycache_prefix`` an existing directory
    that held no file at the capture; else :class:`EnvironmentSealError`."""
    if launch.get("dont_write_bytecode_flag") != 1:
        raise EnvironmentSealError("not launched with python -B (sys.flags.dont_write_bytecode)")
    if not launch.get("pycache_prefix"):
        raise EnvironmentSealError("not launched with -X pycache_prefix=<fresh empty directory>")
    if not launch.get("pycache_prefix_exists"):
        raise EnvironmentSealError("the pycache_prefix directory does not exist")
    if launch.get("pycache_prefix_files_at_start") != 0:
        raise EnvironmentSealError("the pycache_prefix directory was not empty at launch")
    return {**launch, "launch_ok": True}
