import os, sys, time, resource, numpy as np
t_import0 = time.perf_counter()
import warp as wp
from toy_transport import build, tables
prec = sys.argv[1]; N = int(sys.argv[2]); nthreads = int(sys.argv[3]) if len(sys.argv) > 3 else 1
if N % nthreads != 0:
    raise SystemExit(f"N={N} must be divisible by threads={nthreads} so that exactly N histories are simulated")
real = wp.float32 if prec == "f32" else wp.float64; rnp = np.float32 if prec == "f32" else np.float64
wp.init(); t_init = time.perf_counter() - t_import0
kern = build(real)
nx, ny, nz = 120, 120, 300
mat = np.zeros((nx, ny, nz), np.int32); mat[:, :, 100:120] = 1; mat[:60, :, 150:160] = 2
rho = np.ones((nx, ny, nz), rnp); rho[:, :, 100:120] = 1.8; rho[:60, :, 150:160] = 0.3
(rt, it, st), (lE0, idE, lR0, idR, nb) = tables(rnp)
NB = 10
dev = "cpu"
A = dict(mat=wp.array(mat, dtype=wp.int32, device=dev), rho=wp.array(rho, dtype=real, device=dev),
         rt=wp.array(rt, dtype=real, device=dev), it=wp.array(it, dtype=real, device=dev), st=wp.array(st, dtype=real, device=dev))
nsv = (nx // 2) * (ny // 2) * (nz // 2)
def mk():
    return (wp.zeros((NB, nsv), dtype=wp.float32, device=dev), wp.zeros((NB, nsv), dtype=wp.float32, device=dev), wp.zeros(5, dtype=wp.int64, device=dev))
def launch(n, off, d, l, c):
    wp.launch(kern, dim=n, inputs=[wp.uint32(7), off, NB, real(150.0), A["mat"], A["rho"], A["rt"], A["it"], A["st"],
                                   real(lE0), real(idE), real(lR0), real(idR), nb, d, l, c, 100000], device=dev)
d, l, c = mk()
t0 = time.perf_counter(); wp.load_module(kern.module, device=dev) if hasattr(wp, "load_module") else None
t_comp = time.perf_counter() - t0
launch(1, 0, d, l, c); wp.synchronize()
d, l, c = mk()
if nthreads == 1:
    t0 = time.perf_counter(); launch(N, 0, d, l, c); wp.synchronize(); dt = time.perf_counter() - t0
    edep = float(d.numpy().astype(np.float64).sum()); cnt = c.numpy()
else:
    from concurrent.futures import ThreadPoolExecutor
    bufs = [mk() for _ in range(nthreads)]; chunk = N // nthreads
    t0 = time.perf_counter()
    with ThreadPoolExecutor(nthreads) as ex:
        list(ex.map(lambda k: launch(chunk, k * chunk, *bufs[k]), range(nthreads)))
    wp.synchronize(); dt = time.perf_counter() - t0
    edep = sum(float(b[0].numpy().astype(np.float64).sum()) for b in bufs); cnt = sum(b[2].numpy() for b in bufs)
print(f"prec={prec} N={N} simulated={N} threads={nthreads} wall={dt:.2f}s hist/s={N/dt:.0f} edep/primary={edep/N:.2f}MeV "
      f"particles={cnt[0]} truncated={cnt[1]} escaped={cnt[2]} stack_overflow={cnt[3]} steps/primary={cnt[4]/N:.0f} init={t_init:.2f}s load_module={t_comp:.2f}s "
      f"host_peak_MiB={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f}")
