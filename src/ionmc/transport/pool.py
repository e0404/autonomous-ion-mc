"""Process-pool CPU parallelism for the ``python`` and ``warp-cpu`` backends.

Workers are ``spawn``-context ``Process`` objects connected by one-way ``Pipe`` connections
(not ``ProcessPoolExecutor``: a failed or hung worker must be terminated, not waited for).
Worker ``w`` of ``W`` transports the histories ``[floor(w N / W), floor((w + 1) N / W))`` with
private accumulators; its payload is the pickled :class:`~ionmc.config.EffectiveConfig`, and
its answer a :class:`~ionmc.transport.tally.PartialTransport` (exact tally expansions,
counters, deposit grids, diagnostics) that the parent merges in worker order.

Fail closed: if any worker raises, dies without an answer, or the run exceeds
``RunOptions.worker_timeout_s``, every worker is terminated (killed if it ignores the
terminate signal) and :class:`~ionmc.errors.TransportWorkerError` is raised; no partial result
is ever returned. The parent compiles and caches the Warp kernel before spawning (workers
inherit ``WARP_CACHE_PATH``); each worker reports its own kernel load time (``compile_s``), so a
cache miss is visible in the result's ``transport_report``.

Spawned workers re-import the ``__main__`` module: scripts that start a pool need the usual
``if __name__ == "__main__":`` guard.
"""

from __future__ import annotations

import copyreg
import faulthandler
import multiprocessing as mp
import os
import pickle
import time
import traceback
from multiprocessing.connection import Connection, wait
from types import MappingProxyType
from typing import Any

from ionmc.config import EffectiveConfig
from ionmc.errors import TransportWorkerError
from ionmc.transport.tally import PartialTransport


def _read_only_mapping(content: dict[Any, Any]) -> MappingProxyType[Any, Any]:
    return MappingProxyType(content)


copyreg.pickle(MappingProxyType, lambda m: (_read_only_mapping, (dict(m),)))
"""Immutable mappings (materials, table identity) pickle as their dictionary content."""

_KILL_GRACE_S = 5.0


def history_ranges(n_histories: int, workers: int) -> list[tuple[int, int]]:
    """Contiguous ranges ``[floor(w N / W), floor((w + 1) N / W))`` of the ``W`` workers."""
    return [
        ((w * n_histories) // workers, ((w + 1) * n_histories) // workers) for w in range(workers)
    ]


def _worker_main(
    worker: int, h0: int, h1: int, payload: bytes, conn: Connection, fault: str | None
) -> None:
    """Worker entry point: unpickle the configuration, transport ``[h0, h1)``, send the result."""
    faulthandler.enable()  # a fatal signal (SIGSEGV, ...) leaves a Python stack in stderr
    try:
        t0 = time.perf_counter()
        eff: EffectiveConfig = pickle.loads(payload)
        unpickle_s = time.perf_counter() - t0
        if fault == "raise" and worker == 1:
            raise RuntimeError("injected worker failure")
        if fault == "exit" and worker == 1:
            os._exit(3)
        if fault == "hang" and worker == 1:
            time.sleep(3600)
        from ionmc.transport.run import run_range

        part = run_range(eff, h0, h1)
        part.meta.update({"worker": worker, "pid": os.getpid(), "unpickle_s": unpickle_s})
        conn.send(("ok", part))
    except BaseException:
        conn.send(("error", traceback.format_exc()))
    finally:
        conn.close()


def _terminate_all(procs: list[mp.process.BaseProcess]) -> None:
    for p in procs:
        if p.is_alive():
            p.terminate()
    deadline = time.monotonic() + _KILL_GRACE_S
    for p in procs:
        p.join(max(0.0, deadline - time.monotonic()))
    for p in procs:
        if p.is_alive():
            p.kill()
            p.join(_KILL_GRACE_S)


def run_pool(
    eff: EffectiveConfig,
    *,
    h_range: tuple[int, int] | None = None,
    _fault: str | None = None,
) -> list[PartialTransport]:
    """Transport all histories of ``eff`` in ``cpu_workers`` spawned processes.

    ``h_range`` restricts the run to the histories ``[h0, h1)`` (a part of a larger run; the
    batch of history ``h`` stays ``h mod B``); the workers tile that range.
    ``_fault`` ("raise", "exit" or "hang") injects a failure into worker 1; it exists only for
    the fail-closed tests.
    """
    run = eff.requested.run
    workers = run.cpu_workers
    if run.backend == "warp-cuda":
        raise TransportWorkerError("a process pool is not available for the warp-cuda backend")
    parent_compile_s = 0.0
    if run.backend == "warp-cpu":
        from ionmc.transport.warp_driver import load_kernel

        _kernel, parent_compile_s, _regs = load_kernel(eff, "cpu")  # compile once, before spawning
    try:
        payload = pickle.dumps(eff, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:
        raise TransportWorkerError(
            f"the effective configuration cannot be sent to worker processes: {exc}"
        ) from exc
    ctx = mp.get_context("spawn")
    lo, hi = h_range if h_range is not None else (0, run.n_histories)
    if not 0 <= lo < hi <= run.n_histories or workers > hi - lo:
        raise TransportWorkerError(f"invalid history range {lo}..{hi} for {workers} workers")
    ranges = [(lo + a, lo + b) for a, b in history_ranges(hi - lo, workers)]
    procs: list[mp.process.BaseProcess] = []
    conns: dict[int, Connection] = {}
    results: dict[int, PartialTransport] = {}
    t_start = time.monotonic()
    deadline = None if run.worker_timeout_s is None else t_start + run.worker_timeout_s
    try:
        for w, (h0, h1) in enumerate(ranges):
            parent_conn, child_conn = ctx.Pipe(duplex=False)
            proc = ctx.Process(
                target=_worker_main,
                args=(w, h0, h1, payload, child_conn, _fault),
                daemon=True,
                name=f"ionmc-worker-{w}",
            )
            proc.start()
            child_conn.close()
            procs.append(proc)
            conns[w] = parent_conn
        pending = dict(conns)
        while pending:
            timeout = None if deadline is None else max(0.0, deadline - time.monotonic())
            if timeout == 0.0:
                raise TransportWorkerError(
                    f"{len(pending)} of {workers} workers exceeded worker_timeout_s="
                    f"{run.worker_timeout_s}; all workers were terminated, no result"
                )
            by_obj: dict[Any, tuple[int, str]] = {}
            for w, c in pending.items():
                by_obj[c] = (w, "conn")
                by_obj[procs[w].sentinel] = (w, "proc")
            for obj in wait(list(by_obj), timeout):
                w, kind = by_obj[obj]
                if w not in pending:
                    continue
                conn = pending[w]
                if kind == "proc" and not conn.poll():
                    procs[w].join(1.0)
                    raise TransportWorkerError(
                        f"worker {w} (pid {procs[w].pid}) exited with code {procs[w].exitcode} "
                        "without a result; all workers were terminated, no result"
                    )
                try:
                    status, value = conn.recv()
                except EOFError as exc:
                    procs[w].join(1.0)
                    raise TransportWorkerError(
                        f"worker {w} (pid {procs[w].pid}) closed its connection without a "
                        f"result (exit code {procs[w].exitcode}); all workers were terminated, "
                        "no result"
                    ) from exc
                if status != "ok":
                    raise TransportWorkerError(
                        f"worker {w} failed; all workers were terminated, no result:\n{value}"
                    )
                results[w] = value
                del pending[w]
    finally:
        _terminate_all(procs)
        for c in conns.values():
            c.close()
    parts = [results[w] for w in range(workers)]
    elapsed = time.monotonic() - t_start
    for p in parts:
        p.meta["pool_wall_s"] = elapsed
        p.meta["parent_compile_s"] = parent_compile_s  # the compile before the workers started
    return parts
