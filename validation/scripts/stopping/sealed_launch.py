"""Sealed launcher of the amendment I3-FU-1 diagnostics (stdlib only; Codex REVIEW-cf6fa397).

Usage (the only accepted launch of ``i3_band_diagnostic.py`` and ``i3fu2_evaluate.py``)::

    .venv/bin/python -I -S -B -X pycache_prefix="$(mktemp -d)" \
        validation/scripts/stopping/sealed_launch.py <target script> <args...>

The target must be one of :data:`ALLOWED_TARGETS` (exact repository-relative path); any other
path is refused (Codex REVIEW-e0d37216).

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
(e) refuses if any ``numpy`` / ``ionmc`` module, or any module whose file lies under the
    repository's ``src/ionmc`` or ``validation/scripts`` (other than this launcher), is already
    in ``sys.modules`` (no wrapper can pre-import and monkeypatch them); records the module names
    present at launch (``sys_modules_at_launch``);
(f) reads the target's bytes once, records ``target`` (repository-relative), ``target_path``
    and ``target_sha256`` BEFORE executing, and executes exactly those bytes (compiled with the
    target's file name) as ``__main__`` in a fresh module (``__file__`` = the target path), so
    the hashed and the executed bytes are the same.

(0) before anything else (Codex REVIEW-7dd89694): the launcher must BE the interpreter's
    ``__main__`` (its ``__file__`` and ``sys.argv[0]`` are this file, ``__name__`` is
    ``__main__``), its module frame has no caller frame (it was not run in-process by a wrapper),
    and every loaded module with a ``__file__`` lies in the stdlib directory (built-in and frozen
    modules have none), so no wrapper can have run Python code in this process before it;
    ``builtins.compile`` / ``exec`` / ``__import__`` / ``open`` are the genuine built-in functions
    (captured as the very first statements and checked by type, owner module and name;
    ``io.open is builtins.open``; ``sys.modules['builtins']`` / ``['sys']`` are the genuine
    modules);
(g) authenticates the target BEFORE execution against the committed bytes: ``git rev-parse
    HEAD``, ``git status --porcelain --untracked-files=no -- <target>`` empty and ``git show
    HEAD:<target>`` equal to the bytes read from disk (list-form subprocess, no shell); records
    ``head_sha`` and ``target_committed``; compiles and executes exactly those bytes with the
    captured, verified ``compile`` / ``exec``.

The git calls run in a sanitized environment (PATH, ``LANG=C``, a non-existent ``HOME``, no
inherited ``GIT_*`` variable, system and global git configuration disabled) and are bound to the
launcher's own repository (``rev-parse --show-toplevel`` must equal it; the git common dir is
recorded); any git redirection variable (:data:`GIT_REDIRECT_VARIABLES`) in the launcher's own
environment is refused by name (Codex REVIEW-d8ecd044), except ``GIT_WORK_TREE`` whose real path
equals the launcher's repository root (set by the host-runner service; recorded as
``git_work_tree_env``; removed from the git subprocess environment like every ``GIT_*``
variable; Codex REVIEW-427b5254). The bytecode-cache prefix may be absent (a directory that does
not exist cannot hold stale bytecode; with -B nothing is written there) or an existing empty
directory; the state is recorded (``pycache_prefix_state``). The launcher refuses to be imported
as a module (ImportError before any function is defined), so a wrapper cannot reuse its
functions.

Any failed check: a message on stderr and exit status 2; the target is not run.

Trust model: the trusted base is the interpreter binary (hashed in the environment seal and
verified by the diagnostics after start), the git binary used for the bootstrap comparison,
the operating-system kernel and libc. Everything else is either checked before the target runs
(the flags, the empty bytecode prefix, the import machinery, the builtins, the absence of any
earlier Python code, the target's committed bytes) or recorded and verified afterwards against
the seal (the stdlib modules imported before the launcher's checks).

Binding, honestly stated: the interpreter start-up and the stdlib modules imported before this
launcher's own checks (including those it imports itself) are verified only after the fact:
they are recorded here (``stdlib_modules_at_launch``) and checked against the environment seal
by the diagnostics (record-and-detect). The target and everything it imports are bound before
execution: the target by its hash, its imports through the fresh empty bytecode prefix, the
explicit ``sys.path`` and the sealed sources checked by the diagnostics.
"""

from __future__ import annotations

# First statements (Codex REVIEW-7dd89694): capture the builtins used to read, compile and run the
# target and whether this module's frame has a caller, before anything else can run here.
import builtins as _builtins_at_start

_COMPILE = _builtins_at_start.compile
_EXEC = _builtins_at_start.exec
_OPEN = _builtins_at_start.open
_IMPORT = _builtins_at_start.__import__
_HAS_CALLER_FRAME = _IMPORT("sys")._getframe(0).f_back is not None
if __name__ != "__main__":  # Codex REVIEW-d8ecd044: never usable as a library from a wrapper
    raise ImportError("sealed_launch.py must be run as a script; it refuses to be imported")

import builtins  # noqa: E402
import hashlib  # noqa: E402
import importlib.machinery  # noqa: E402
import io  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import sysconfig  # noqa: E402
import types  # noqa: E402
import zipimport  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

RECORD_NAME = "__ionmc_sealed_launch__"
ORIGIN = "validation/scripts/stopping/sealed_launch.py"
FILEFINDER_HOOK = (
    "_frozen_importlib_external",
    "FileFinder.path_hook.<locals>.path_hook_for_FileFinder",
    "<frozen importlib._bootstrap_external>",
)
REPO = Path(__file__).resolve().parents[3]
ALLOWED_TARGETS = (
    "validation/scripts/stopping/i3_band_diagnostic.py",
    "validation/scripts/stopping/i3fu2_evaluate.py",
)


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
    state = prefix_state()
    if state == "unset":
        out.append("sys.pycache_prefix is not set (-X pycache_prefix=<fresh empty directory>)")
    elif state == "not-a-directory":
        out.append("the pycache_prefix path exists but is not a directory")
    elif state == "non-empty":
        out.append("the pycache_prefix directory is not empty")
    return out


def prefix_state() -> str:
    """``absent`` (a directory that does not exist cannot hold stale bytecode; with -B nothing is
    written there), ``empty``, ``non-empty``, ``not-a-directory`` or ``unset``."""
    prefix = getattr(sys, "pycache_prefix", None)
    if not prefix:
        return "unset"
    path = Path(prefix)
    if not path.exists():
        return "absent"
    if not path.is_dir():
        return "not-a-directory"
    return "non-empty" if any(p.is_file() for p in path.rglob("*")) else "empty"


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


def bootstrap_problems() -> list[str]:
    """(0) No wrapper ran Python code in this process before the launcher, and the builtins
    used to read, compile and run the target are genuine (Codex REVIEW-7dd89694)."""
    out = []
    own = Path(__file__).resolve()
    main = sys.modules.get("__main__")
    main_file = getattr(main, "__file__", None)
    if __name__ != "__main__" or main_file is None or Path(main_file).resolve() != own:
        out.append("the launcher is not the interpreter's __main__ module")
    if not sys.argv or Path(sys.argv[0]).resolve() != own:
        out.append("sys.argv[0] is not the launcher")
    if _HAS_CALLER_FRAME:
        out.append("the launcher's module frame has a caller (run in-process by a wrapper)")
    base = stdlib_dir()
    foreign = sorted(
        name
        for name, mod in list(sys.modules.items())
        if name != "__main__"
        and getattr(mod, "__file__", None)
        and not Path(str(mod.__file__)).resolve().is_relative_to(base)
    )
    if foreign:
        out.append(f"non-stdlib modules loaded before the launcher: {foreign}")
    want = {
        "compile": (_COMPILE, builtins),
        "exec": (_EXEC, builtins),
        "__import__": (_IMPORT, builtins),
        "open": (_OPEN, sys.modules.get("_io")),
    }
    for name, (func, owner) in want.items():
        genuine = (
            type(func) is types.BuiltinFunctionType
            and getattr(func, "__self__", None) is owner
            and getattr(func, "__name__", None) == name
            and getattr(builtins, name) is func
        )
        if not genuine:
            out.append(f"builtins.{name} is not the genuine built-in function")
    if sys.modules.get("builtins") is not builtins or sys.modules.get("sys") is not sys:
        out.append("sys.modules['builtins'] / ['sys'] are not the genuine modules")
    if io.open is not _OPEN:
        out.append("io.open is not builtins.open")
    return out


#: git variables that redirect the repository, index or object store: their presence in the
#: launcher's own environment is refused (fail closed, the variable named)
GIT_REDIRECT_VARIABLES = (
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_COMMON_DIR",
    "GIT_NAMESPACE",
)


def git_environment_problems(environ: Any = None) -> list[str]:
    """(g) no git redirection variable may be set in the launcher's own environment, except
    ``GIT_WORK_TREE`` whose real path equals the launcher's repository root (the host-runner
    service sets it to the run's checkout; Codex REVIEW-427b5254)."""
    env = os.environ if environ is None else environ
    out = []
    for v in GIT_REDIRECT_VARIABLES:
        if v not in env:
            continue
        if v == "GIT_WORK_TREE" and Path(env[v]).resolve() == REPO:
            continue
        if v == "GIT_WORK_TREE":
            out.append(
                f"GIT_WORK_TREE={env[v]!r} is not the launcher's repository {REPO} (refused)"
            )
        else:
            out.append(f"{v} is set in the launcher's environment (refused)")
    return out


def git_env() -> dict[str, str]:
    """The sanitized environment of the git calls: PATH, LANG=C, a non-existent HOME, no
    inherited ``GIT_*`` variable; system and global git configuration disabled."""
    return {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "LANG": "C",
        "HOME": "/nonexistent-ionmc-sealed-launch-home",
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_GLOBAL": os.devnull,
    }


def committed_target(rel: str) -> tuple[bytes, dict[str, str]]:
    """(g): the target's bytes committed at HEAD of the launcher's own repository, via the git
    binary (list-form subprocess, no shell, :func:`git_env`). The repository is bound: ``git -C
    <repo> rev-parse --show-toplevel`` must be the repository three levels above the launcher
    and ``--git-common-dir`` an existing directory (recorded); the target must have no
    uncommitted change. Returns (bytes, {head_sha, git, toplevel, git_common_dir});
    RuntimeError otherwise."""
    problems = git_environment_problems()
    if problems:
        raise RuntimeError("; ".join(problems))
    git = shutil.which("git", path=git_env()["PATH"])
    if git is None:
        raise RuntimeError("git not found (needed to authenticate the target)")

    def run(*args: str) -> bytes:
        proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
            [git, "-C", str(REPO), *args], capture_output=True, check=False, env=git_env()
        )
        if proc.returncode != 0:
            raise RuntimeError(f"git {' '.join(args)} failed: {proc.stderr.decode().strip()}")
        return proc.stdout

    toplevel = Path(run("rev-parse", "--show-toplevel").decode().strip()).resolve()
    if toplevel != REPO:
        raise RuntimeError(f"git toplevel {toplevel} is not the launcher's repository {REPO}")
    common = Path(
        run("rev-parse", "--path-format=absolute", "--git-common-dir").decode().strip()
    ).resolve()
    if not common.is_dir():
        raise RuntimeError(f"git common dir {common} is not a directory")
    head = run("rev-parse", "HEAD").decode().strip()
    if run("status", "--porcelain", "--untracked-files=no", "--", rel).strip():
        raise RuntimeError(f"the target {rel} has an uncommitted change")
    binding = {
        "head_sha": head,
        "git": git,
        "git_toplevel": str(toplevel),
        "git_common_dir": str(common),
    }
    return run("show", f"HEAD:{rel}"), binding


def preimport_problems(modules: dict[str, Any] | None = None) -> list[str]:
    """(e): no numpy / ionmc / repository-script module may be loaded before the target."""
    mods = sys.modules if modules is None else modules
    own = Path(__file__).resolve()
    roots = (REPO / "src" / "ionmc", REPO / "validation" / "scripts")
    bad = []
    for name, mod in list(mods.items()):
        if name in ("numpy", "ionmc") or name.startswith(("numpy.", "ionmc.")):
            bad.append(name)
            continue
        f = getattr(mod, "__file__", None)
        if f:
            path = Path(f).resolve()
            if path != own and any(path.is_relative_to(r) for r in roots):
                bad.append(name)
    return [f"modules already imported before the launch: {sorted(bad)}"] if bad else []


def target_rel(target: Path) -> str:
    """The repository-relative path of an allowed target; anything else raises RuntimeError."""
    try:
        rel = target.resolve().relative_to(REPO).as_posix()
    except ValueError:
        rel = ""
    if rel not in ALLOWED_TARGETS:
        raise RuntimeError(f"target {target} is not one of {ALLOWED_TARGETS}")
    return rel


def launch_record(target: Path, blob: bytes, binding: dict[str, str]) -> dict[str, Any]:
    """Run checks (a)-(e); raise RuntimeError on any problem; return the record."""
    problems = flag_problems() + machinery_problems() + preimport_problems()
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
        "target": target_rel(target),
        "target_path": str(target.resolve()),
        "target_sha256": hashlib.sha256(blob).hexdigest(),
        **binding,
        "target_committed": True,
        "sys_modules_at_launch": sorted(sys.modules),
        "python_executable": str(Path(sys.executable).resolve()),
        "flags": {
            "isolated": sys.flags.isolated,
            "no_site": sys.flags.no_site,
            "dont_write_bytecode": sys.flags.dont_write_bytecode,
        },
        "pycache_prefix": sys.pycache_prefix,
        "pycache_prefix_state": prefix_state(),
        "git_work_tree_env": os.environ.get("GIT_WORK_TREE"),
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
        problems = (
            bootstrap_problems()
            + flag_problems()
            + machinery_problems()
            + git_environment_problems()
        )
        if problems:
            raise RuntimeError("; ".join(problems))
        rel = target_rel(target)
        with _OPEN(target, "rb") as fh:  # the verified builtin open
            blob = fh.read()  # hashed, authenticated and executed bytes are the same
        committed, binding = committed_target(rel)
        if blob != committed:
            raise RuntimeError(f"the target {rel} differs from the committed bytes at HEAD")
        record = launch_record(target, blob, binding)
        code = _COMPILE(blob, str(target), "exec", dont_inherit=True)
    except (RuntimeError, OSError, SyntaxError) as exc:
        print(f"sealed_launch: refused: {exc}", file=sys.stderr)
        return 2
    setattr(builtins, RECORD_NAME, record)
    sys.argv = [str(target), *argv[1:]]
    module = types.ModuleType("__main__")
    module.__file__ = str(target)
    module.__builtins__ = builtins  # type: ignore[attr-defined]
    sys.modules["__main__"] = module
    _EXEC(code, module.__dict__)  # the verified builtin exec; the authenticated target
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
