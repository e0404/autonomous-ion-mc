"""Code, environment, data and model provenance of a persisted result (decision 0045 sections 2-3).

* :func:`package_tree_oid` is the Git *tree object id* of the imported ``ionmc`` package directory
  (SHA-1 blob/tree encoding exactly as git writes it; ``__pycache__``, ``*.pyc`` and the build stamp
  are excluded). It needs no git binary and equals ``git rev-parse <commit>:src/ionmc`` exactly when
  the code that ran equals that commit's package tree.
* :func:`code_identity` follows decision 0045 section 3: (a) build stamp, (b) one read-only,
  timeout-bounded ``git rev-parse`` with every inherited ``GIT_*`` variable removed, (c) unknown.
  This query is the only subprocess the library ever runs; the library never uses the network.
* :func:`privacy_violations` enforces the frozen path set P (decision 0045 section 2).
"""

from __future__ import annotations

import getpass
import hashlib
import json
import os
import platform
import re
import socket
import stat
import subprocess
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy

from ionmc._version import __version__

BUILD_STAMP_NAME = "_build_info.json"
GIT_TIMEOUT_S = 5.0
_PACKAGE_DIR = Path(__file__).resolve().parent
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_EXCLUDED_DIRS = frozenset({"__pycache__"})


# ---------------------------------------------------------------------------------------------
# package tree object id


def _blob_oid(data: bytes) -> bytes:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).digest()


def _tree_oid(directory: Path, top: bool) -> bytes | None:
    """Raw oid of the git tree of ``directory`` (None when it would be an empty tree)."""
    entries: list[tuple[bytes, bytes, bytes, bytes]] = []  # (sort key, mode, name, raw oid)
    for child in directory.iterdir():
        name = child.name
        if name in _EXCLUDED_DIRS or name.endswith(".pyc") or (top and name == BUILD_STAMP_NAME):
            continue
        raw = os.fsencode(name)
        st = child.lstat()
        if stat.S_ISLNK(st.st_mode):
            entries.append((raw, b"120000", raw, _blob_oid(os.fsencode(os.readlink(child)))))
        elif stat.S_ISDIR(st.st_mode):
            sub = _tree_oid(child, False)
            if sub is not None:
                entries.append((raw + b"/", b"40000", raw, sub))  # git sorts a tree as "name/"
        elif stat.S_ISREG(st.st_mode):
            mode = b"100755" if st.st_mode & stat.S_IXUSR else b"100644"
            entries.append((raw, mode, raw, _blob_oid(child.read_bytes())))
    if not entries:
        return None
    entries.sort(key=lambda e: e[0])
    body = b"".join(mode + b" " + name + b"\0" + oid for _, mode, name, oid in entries)
    return hashlib.sha1(b"tree %d\0" % len(body) + body).digest()


def package_tree_oid(package_dir: str | os.PathLike[str] | None = None) -> str:
    """Git tree object id (40 hex) of the package directory (default: the imported ``ionmc``)."""
    root = Path(package_dir) if package_dir is not None else _PACKAGE_DIR
    oid = _tree_oid(root, True)
    return (oid if oid is not None else hashlib.sha1(b"tree 0\0").digest()).hex()


# ---------------------------------------------------------------------------------------------
# code identity


def _git_env() -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env["LC_ALL"] = "C"
    return env


def _read_stamp(package_dir: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((package_dir / BUILD_STAMP_NAME).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _git_query(package_dir: Path) -> tuple[str, str] | str:
    """One read-only ``git rev-parse HEAD HEAD:src/ionmc --show-prefix``; (commit, tree oid) or
    a reason string. The work-tree prefix must be ``src/ionmc/``: a package installed elsewhere
    (for example a virtual environment inside some repository) is not that repository's package."""
    try:
        proc = subprocess.run(  # noqa: S603 - fixed argument list, no shell
            ["git", "rev-parse", "HEAD", "HEAD:src/ionmc", "--show-prefix"],  # noqa: S607
            cwd=package_dir,
            env=_git_env(),
            capture_output=True,
            text=True,
            timeout=GIT_TIMEOUT_S,
            check=False,
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"git query failed ({type(exc).__name__})"
    if proc.returncode != 0:
        return "the package is not in a Git work tree at src/ionmc (git rev-parse failed)"
    lines = proc.stdout.split("\n")
    if len(lines) < 3 or lines[2] != "src/ionmc/":
        return "the package directory is not src/ionmc of a Git work tree"
    if not (_HEX40.match(lines[0]) and _HEX40.match(lines[1])):
        return "git rev-parse returned an unexpected answer"
    return lines[0], lines[1]


def _code(
    version: str,
    oid: str | None,
    commit: str | None,
    dirty: bool | None,
    source: str,
    reason: str,
) -> dict[str, Any]:
    ok = commit is not None and dirty is not None
    return {
        "package_tree_oid": oid,
        "ionmc_version": version,
        "git_commit": commit,
        "git_dirty": dirty,
        "source": source,
        "status": "ok" if ok else "unknown",
        "reason": "" if ok else reason,
    }


def code_identity(package_dir: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """The ``provenance.code`` record (decision 0045 section 3, (a) stamp, (b) git, (c) unknown)."""
    pkg = Path(package_dir) if package_dir is not None else _PACKAGE_DIR
    try:
        oid: str | None = package_tree_oid(pkg)
    except OSError as exc:
        oid = None
        why = f"package tree not readable ({type(exc).__name__})"
        return _code(__version__, None, None, None, "unavailable", why)
    stamp = _read_stamp(pkg)
    if stamp is not None:
        commit, stamped = stamp.get("git_commit"), stamp.get("package_tree_oid")
        commit = commit if isinstance(commit, str) and _HEX40.match(commit) else None
        dirty = (oid != stamped) if isinstance(stamped, str) and _HEX40.match(stamped) else None
        if stamp.get("status") == "dirty" and dirty is not None:
            dirty = True  # built from a tree that differs from the commit's package tree
        return _code(
            __version__, oid, commit, dirty, "stamp",
            "the build stamp carries no commit or no tree id",
        )  # fmt: skip
    found = _git_query(pkg)
    if isinstance(found, str):
        return _code(__version__, oid, None, None, "unavailable", f"no build stamp; {found}")
    commit, head_tree = found
    return _code(__version__, oid, commit, head_tree != oid, "git", "")


# ---------------------------------------------------------------------------------------------
# environment, data, models


def environment_block() -> dict[str, Any]:
    """Software environment; never a hostname, user name, home path or absolute path."""
    from ionmc.environment import _describe_warp  # local import: environment imports this module

    warp_info = _describe_warp()
    devices = [dict(d) for d in warp_info["devices"]]
    if any(d.get("is_cuda") for d in devices):
        try:
            import warp

            driver = warp.get_cuda_driver_version()
        except Exception:  # noqa: BLE001 - diagnostics must never raise
            driver = None
        for d in devices:
            if d.get("is_cuda"):
                d["cuda_driver_version"] = None if driver is None else str(driver)
    return {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "numpy_version": numpy.__version__,
        "warp_version": warp_info["version"],
        "warp_devices": devices,
        "warp_error": warp_info["error"],
    }


def _strip(obj: Any, keys: frozenset[str]) -> Any:
    if isinstance(obj, Mapping):
        return {k: _strip(v, keys) for k, v in obj.items() if k not in keys}
    if isinstance(obj, list | tuple):
        return [_strip(v, keys) for v in obj]
    return obj


def data_block(
    effective_summary: Mapping[str, Any],
    lookups: Sequence[Mapping[str, Any]],
    nuclear_setup: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Dataset identities (never values, never cache paths): every ``stopping.data_used`` pin with
    archive/member sha256 and per-run role, table identities, the nuclear setup, lookup digests."""
    stopping = effective_summary.get("stopping") or {}
    tables = effective_summary["tables"]
    nuc = nuclear_setup if nuclear_setup is not None else effective_summary.get("nuclear")
    drop = frozenset({"cache_dir"})
    return {
        "stopping_data_used": _strip(stopping.get("data_used", []), drop),
        "stopping_content_sha256": stopping.get("content_sha256"),
        "table_sha256": tables["sha256"],
        "table_stopping_source": tables["stopping_source"],
        "table_materials": tables["materials"],
        "nuclear": None if nuc is None else _strip(nuc, drop),
        "lookups": [dict(lk) for lk in lookups],
    }


def models_block(effective_summary: Mapping[str, Any]) -> dict[str, Any]:
    """Model identities; ``evidence_at_run`` is copied verbatim from the run, never recomputed."""
    st = effective_summary.get("stopping") or {}
    ph = effective_summary["physics"]
    nuc = effective_summary.get("nuclear") or {}
    return {
        "stopping_model_id_effective": st.get("model_id_effective"),
        "stopping_version": st.get("version"),
        "stopping_spec_sha256": st.get("spec_sha256"),
        "evidence_at_run": st.get("evidence_at_run"),
        "straggling_model": ph["straggling_model"] if ph["straggling"] else None,
        "mcs_model": ph["mcs_model"] if ph["multiple_scattering"] else None,
        "end_of_range_model": ph["end_of_range_model_effective"],
        "delta_electrons": ph["delta_electrons"],
        "nuclear_table_id": nuc.get("table_id"),
    }


# ---------------------------------------------------------------------------------------------
# privacy: the frozen path set P

#: Regular expressions over JSON pointers: the only places that may hold an absolute path.
PATH_POINTERS: tuple[str, ...] = (
    r"^/run/output_path$",
    r"^/run/cache_dir$",
    r"^/input_manifest/base_dir$",
    r"^/input_manifest/config_file/path$",
    r"^/input_manifest/files/\d+/path$",
    r"^/effective/(stopping|nuclear)(/.*)?/cache_dir$",
)
_P_PATTERNS = tuple(re.compile(p) for p in PATH_POINTERS)
#: Blocks in which a bare string equal to the host or user name is a leak (an identity block, not an
#: enumerated vocabulary such as ``mcs_mapping.evaluation = "root"``).
_IDENTITY_POINTERS = re.compile(r"^/(provenance/(code|environment)|run/device)(/|$)")


def _walk_strings(obj: Any, ptr: str = "") -> Iterator[tuple[str, str]]:
    if isinstance(obj, str):
        yield ptr, obj
    elif isinstance(obj, Mapping):
        for k, v in obj.items():
            yield from _walk_strings(v, f"{ptr}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk_strings(v, f"{ptr}/{i}")


def _identity_names() -> set[str]:
    names: set[str] = set()
    for fn in (getpass.getuser, socket.gethostname):
        try:
            names.add(fn())
        except Exception:  # noqa: BLE001 - no user database entry in a container
            pass
    names.add(os.uname().nodename)
    names.discard("")
    return names


def privacy_violations(
    result_doc: Mapping[str, Any],
    *,
    names: Sequence[str] | None = None,
    home: str | None = None,
) -> list[str]:
    """Pointers (with a reason) that break the path-set rule of decision 0045 section 2.

    Outside P no string contains the home path (substring); a string that is a path (it has a
    ``/`` or ``\\``) has no component equal to the user or host name (whole components only, never
    substrings); in the identity blocks no bare string equals the user or host name."""
    ident = set(names) if names is not None else _identity_names()
    home_s = home if home is not None else str(Path.home())
    out: list[str] = []
    for ptr, s in _walk_strings(result_doc):
        if any(p.match(ptr) for p in _P_PATTERNS):
            continue
        if home_s not in ("", "/") and home_s in s:
            out.append(f"{ptr}: holds the home path")
            continue
        if "/" in s or "\\" in s:
            if ident & set(re.split(r"[/\\]", s)):
                out.append(f"{ptr}: a path component equals the user or host name")
        elif s in ident and _IDENTITY_POINTERS.match(ptr):
            out.append(f"{ptr}: equals the user or host name")
    return out


# ---------------------------------------------------------------------------------------------
# host-runner expectation (row PROV-7; the assertion is run by the orchestrator on the host)


def host_runner_identity_expectation(
    code_sha: str, package_tree_oid: str | None = None
) -> dict[str, Any]:
    """Expected ``provenance.code`` fields of a bundle written under the controlled host runner at
    the clean committed ``code_sha`` (``package_tree_oid`` = ``git rev-parse <sha>:src/ionmc``)."""
    exp: dict[str, Any] = {
        "source": "git",
        "status": "ok",
        "git_dirty": False,
        "git_commit": code_sha,
    }
    if package_tree_oid is not None:
        exp["package_tree_oid"] = package_tree_oid
    return exp


def host_runner_identity_mismatches(
    code: Mapping[str, Any], expected: Mapping[str, Any]
) -> list[str]:
    """Fields of ``code`` that differ from :func:`host_runner_identity_expectation`."""
    return [f"{k}: {code.get(k)!r} != {v!r}" for k, v in expected.items() if code.get(k) != v]
