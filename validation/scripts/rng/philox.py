import warp as wp, numpy as np, time
wp.config.quiet = True; wp.init()
M0 = wp.constant(wp.uint64(0xD2511F53)); M1 = wp.constant(wp.uint64(0xCD9E8D57))
W0 = wp.constant(wp.uint32(0x9E3779B9)); W1 = wp.constant(wp.uint32(0xBB67AE85))
MASK = wp.constant(wp.uint64(0xFFFFFFFF))

@wp.func
def philox4x32_10(c: wp.vec4ui, k: wp.vec2ui) -> wp.vec4ui:
    c0 = c[0]; c1 = c[1]; c2 = c[2]; c3 = c[3]; k0 = k[0]; k1 = k[1]
    for r in range(10):
        p0 = M0 * wp.uint64(c0); p1 = M1 * wp.uint64(c2)
        hi0 = wp.uint32(p0 >> wp.uint64(32)); lo0 = wp.uint32(p0 & MASK)
        hi1 = wp.uint32(p1 >> wp.uint64(32)); lo1 = wp.uint32(p1 & MASK)
        n0 = hi1 ^ c1 ^ k0; n2 = hi0 ^ c3 ^ k1
        c0 = n0; c1 = lo1; c2 = n2; c3 = lo0
        k0 = k0 + W0; k1 = k1 + W1
    return wp.vec4ui(c0, c1, c2, c3)

@wp.func
def u01_24(x: wp.uint32) -> float:   # (0,1): never 0 -> safe for log()
    return (float(x >> wp.uint32(8)) + 0.5) * (1.0 / 16777216.0)

@wp.kernel
def kat(cs: wp.array(dtype=wp.vec4ui), ks: wp.array(dtype=wp.vec2ui), out: wp.array(dtype=wp.vec4ui)):
    i = wp.tid(); out[i] = philox4x32_10(cs[i], ks[i])

cs = np.array([[0,0,0,0],[0xffffffff]*4,[0x243f6a88,0x85a308d3,0x13198a2e,0x03707344]], dtype=np.uint32)
ks = np.array([[0,0],[0xffffffff]*2,[0xa4093822,0x299f31d0]], dtype=np.uint32)
exp = [[0x6627e8d5,0xe169c58d,0xbc57ac4c,0x9b00dbd8],[0x408f276d,0x41c83b0e,0xa20bc7c6,0x6d5451fd],[0xd16cfe09,0x94fdcceb,0x5001e420,0x24126ea1]]
o = wp.zeros(3, dtype=wp.vec4ui)
wp.launch(kat, dim=3, inputs=[wp.array(cs, dtype=wp.vec4ui), wp.array(ks, dtype=wp.vec2ui), o], device="cpu")
print("Random123 KAT pass:", (o.numpy() == np.array(exp, dtype=np.uint32)).all())

def philox_py(c, k):   # pure-Python reference adapter (ints)
    c0,c1,c2,c3 = c; k0,k1 = k
    for _ in range(10):
        p0 = 0xD2511F53*c0; p1 = 0xCD9E8D57*c2
        c0, c1, c2, c3 = ((p1>>32) ^ c1 ^ k0) & 0xFFFFFFFF, p1 & 0xFFFFFFFF, ((p0>>32) ^ c3 ^ k1) & 0xFFFFFFFF, p0 & 0xFFFFFFFF
        k0 = (k0 + 0x9E3779B9) & 0xFFFFFFFF; k1 = (k1 + 0xBB67AE85) & 0xFFFFFFFF
    return (c0,c1,c2,c3)
print("python adapter KAT pass:", all(philox_py(tuple(map(int,c)), tuple(map(int,k))) == tuple(e) for c,k,e in zip(cs,ks,exp)))
# bit-equality on random inputs
rng = np.random.default_rng(0); n = 20000
C = rng.integers(0, 2**32, (n,4), dtype=np.uint64).astype(np.uint32); K = rng.integers(0, 2**32, (n,2), dtype=np.uint64).astype(np.uint32)
o = wp.zeros(n, dtype=wp.vec4ui); wp.launch(kat, dim=n, inputs=[wp.array(C, dtype=wp.vec4ui), wp.array(K, dtype=wp.vec2ui), o], device="cpu")
on = o.numpy(); t=time.perf_counter()
ok = all(philox_py(tuple(map(int,C[i])), tuple(map(int,K[i]))) == tuple(map(int,on[i])) for i in range(n))
tp = (time.perf_counter()-t)/n
print(f"kernel==python on {n} random blocks: {ok}; python adapter {tp*1e6:.1f} us/block")

# cost: N threads x D uniforms
@wp.kernel
def cost_philox(D: int, seed: wp.uint32, out: wp.array(dtype=float)):
    i = wp.tid(); acc = float(0.0)
    key = wp.vec2ui(seed, wp.uint32(0))
    for j in range(D // 4):
        r = philox4x32_10(wp.vec4ui(wp.uint32(i), wp.uint32(j), wp.uint32(0), wp.uint32(0)), key)
        acc += u01_24(r[0]) + u01_24(r[1]) + u01_24(r[2]) + u01_24(r[3])
    out[i] = acc
@wp.kernel
def cost_pcg(D: int, seed: int, out: wp.array(dtype=float)):
    i = wp.tid(); acc = float(0.0); s = wp.rand_init(seed, i)
    for j in range(D):
        acc += wp.randf(s)
    out[i] = acc
N, D = 20000, 2000
out = wp.zeros(N, dtype=float)
for name, kk, args in [("philox", cost_philox, [D, wp.uint32(7), out]), ("warp-pcg", cost_pcg, [D, 7, out])]:
    wp.launch(kk, dim=1, inputs=args, device="cpu")  # compile/warm
    t = time.perf_counter(); wp.launch(kk, dim=N, inputs=args, device="cpu"); wp.synchronize()
    dt = time.perf_counter()-t
    print(f"{name}: {dt/(N*D)*1e9:.2f} ns/uniform (Warp CPU, 1 thread), mean={out.numpy().mean()/D:.5f}")
