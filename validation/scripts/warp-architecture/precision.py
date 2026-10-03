import numpy as np
rng = np.random.default_rng(1)
print("== float32 sequential accumulation into one voxel (worst case: hottest Bragg-peak voxel) ==")
for n in [10**5, 10**6, 3*10**6, 10**7]:
    x = rng.uniform(0.5, 8.0, n)                     # MeV per deposit
    exact = np.sum(x.astype(np.float32).astype(np.float64))   # exact sum of the f32-rounded inputs
    seq32 = np.add.accumulate(x.astype(np.float32), dtype=np.float32)[-1]   # sequential f32 (atomicAdd-like)
    print(f"n={n:>9d} rel.err f32 sequential={abs(seq32-exact)/exact:.2e}  (MC rel. stat. err of mean ~ {np.std(x)/np.mean(x)/np.sqrt(n):.1e})")
print("== same, split over 40 batch accumulators, then f64 sum of batches ==")
n = 10**7; x = rng.uniform(0.5, 8.0, n).astype(np.float32); exact = x.astype(np.float64).sum()
b = sum(np.float64(np.add.accumulate(x[i::40], dtype=np.float32)[-1]) for i in range(40))
print(f"n={n} 40 batches: rel.err={abs(b-exact)/exact:.2e}")

print("== position/direction accumulation over 300 steps: float32 vs float64 ==")
errs = []
for h in range(2000):
    r = np.random.default_rng(h)
    p64 = np.array([60.3, 60.7, 0.0]); u64 = np.array([0.0, 0.0, 1.0])
    p32 = p64.astype(np.float32); u32 = u64.astype(np.float32)
    for k in range(300):
        s = r.uniform(0.05, 1.0); th = r.normal(0, 0.003); ph = r.uniform(0, 2*np.pi)
        for (p, u, dt) in ((p64, u64, np.float64), (p32, u32, np.float32)):
            p += (dt(s) * u).astype(dt)
            ax = np.cross(u, [1.0, 0, 0]); ax = ax / np.linalg.norm(ax)
            v = (u * dt(np.cos(th)) + dt(np.sin(th)) * (np.cos(ph) * ax + np.sin(ph) * np.cross(u, ax))).astype(dt)
            u[:] = v / np.linalg.norm(v)
    errs.append(np.abs(p32.astype(np.float64) - p64).max())
errs = np.array(errs)
print(f"final |p32-p64| over 2000 paths: median={np.median(errs):.2e} mm max={errs.max():.2e} mm; float32 ulp at 300 mm={np.spacing(np.float32(300.0)):.1e} mm")

print("== energy near end of range: E from residual range R (Bragg-Kleeman R=a E^p, water) ==")
a, pw = 0.0225, 1.77
for E in [150.0, 10.0, 1.0, 0.1]:
    R = a * E**pw
    for ds in [1.0, 1e-2, 1e-4]:
        if ds >= R: continue
        e64 = ((R - ds) / a) ** (1/pw); dE64 = E - e64
        R32 = np.float32(a) * np.float32(E) ** np.float32(pw)
        e32 = ((R32 - np.float32(ds)) / np.float32(a)) ** np.float32(1/pw); dE32 = np.float32(E) - e32
        print(f"E={E:6.1f} MeV ds={ds:.0e} mm: dE64={dE64:.4e} dE32={float(dE32):.4e} rel.err={abs(float(dE32)-dE64)/dE64:.1e}")
