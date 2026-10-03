# Quantify state-space overlap of Warp's built-in RNG (32-bit state, state=pcg_hash(state)).
import warp as wp, numpy as np, time, sys
wp.config.quiet = True; wp.init()

@wp.kernel
def mark(seed: int, L: int, bitmap: wp.array(dtype=wp.uint32), hits: wp.array(dtype=wp.int64)):
    i = wp.tid()
    s = wp.rand_init(seed, i)
    nhit = wp.int64(0)
    for k in range(L):
        u = wp.randu(s)          # advances state; s is now the new state
        w = int(s >> wp.uint32(5)); b = s & wp.uint32(31)
        old = wp.atomic_or(bitmap, w, wp.uint32(1) << b)
        if (old >> b) & wp.uint32(1) == wp.uint32(1):
            nhit += wp.int64(1)
    wp.atomic_add(hits, 0, nhit)

N = int(sys.argv[1]); L = int(sys.argv[2])
bitmap = wp.zeros(1 << 27, dtype=wp.uint32, device="cpu")
hits = wp.zeros(1, dtype=wp.int64, device="cpu")
t = time.perf_counter()
wp.launch(mark, dim=N, inputs=[1234, L, bitmap, hits], device="cpu"); wp.synchronize()
h = int(hits.numpy()[0]); tot = N * L
print(f"N={N} L={L} draws={tot:.3e} fraction_of_2^32={tot/2**32:.3f} repeated_states={h:.3e} ({h/tot:.3%}) t={time.perf_counter()-t:.1f}s")
