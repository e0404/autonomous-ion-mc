import warp as wp, numpy as np
wp.config.quiet = True
wp.init()

@wp.func
def mix32(x: wp.uint32) -> wp.uint32:
    b = x * wp.uint32(747796405) + wp.uint32(2891336453)
    c = ((b >> ((b >> wp.uint32(28)) + wp.uint32(4))) ^ b) * wp.uint32(277803737)
    return (c >> wp.uint32(22)) ^ c

@wp.func
def u01(x: wp.uint32) -> float:
    return float(x >> wp.uint32(8)) * (1.0 / 16777216.0)

@wp.func
def mulhilo(a: wp.uint32, b: wp.uint32) -> wp.vec2ui:
    p = wp.uint64(a) * wp.uint64(b)
    return wp.vec2ui(wp.uint32(p >> wp.uint64(32)), wp.uint32(p & wp.uint64(0xFFFFFFFF)))

@wp.kernel
def k(xs: wp.array(dtype=wp.uint32), out: wp.array(dtype=wp.uint32), outf: wp.array(dtype=float), hl: wp.array(dtype=wp.vec2ui)):
    i = wp.tid()
    out[i] = mix32(xs[i])
    outf[i] = u01(out[i])
    hl[i] = mulhilo(xs[i], wp.uint32(0xD2511F53))

xs = np.array([0, 1, 12345, 0xFFFFFFFF, 0x80000000, 987654321], dtype=np.uint32)
a = wp.array(xs, dtype=wp.uint32); o = wp.zeros(len(xs), dtype=wp.uint32); of = wp.zeros(len(xs), dtype=float); hl = wp.zeros(len(xs), dtype=wp.vec2ui)
wp.launch(k, dim=len(xs), inputs=[a, o, of, hl], device="cpu")
for i, x in enumerate(xs):
    try:
        py = mix32(wp.uint32(int(x)))
        pf = u01(py)
        ph = mulhilo(wp.uint32(int(x)), wp.uint32(0xD2511F53))
        print(int(x), "kernel", int(o.numpy()[i]), "py", int(py.value if hasattr(py,'value') else py), type(py).__name__,
              "f", of.numpy()[i], pf, "hl", hl.numpy()[i], ph)
    except Exception as e:
        print("pyscope error", type(e).__name__, e)

# rand_init from python scope?
for fn in ("rand_init",):
    try:
        print("rand_init py:", wp.rand_init(1, 2))
    except Exception as e:
        print("rand_init py error:", type(e).__name__, e)

# array indexing in python scope
@wp.func
def look(t: wp.array(dtype=float), i: int) -> float:
    return t[i]
try:
    print(look(wp.array([1.0,2.0],dtype=float), 1))
except Exception as e:
    print("array index pyscope error:", type(e).__name__, e)

# local fixed-size stack via vector with dynamic index (python scope + kernel)
vec16 = wp.types.vector(length=16, dtype=wp.float32)
@wp.func
def stack_demo(n: int) -> float:
    s = vec16()
    top = int(0)
    for j in range(n):
        s[top] = float(j) * 2.0
        top = top + 1
    acc = float(0.0)
    while top > 0:
        top = top - 1
        acc = acc + s[top]
    return acc
@wp.kernel
def ks(out: wp.array(dtype=float)):
    out[0] = stack_demo(10)
o2 = wp.zeros(1, dtype=float); wp.launch(ks, dim=1, inputs=[o2], device="cpu")
print("stack kernel", o2.numpy()[0])
try:
    print("stack pyscope", stack_demo(10))
except Exception as e:
    print("stack pyscope error", type(e).__name__, e)
