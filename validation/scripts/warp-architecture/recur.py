import warp as wp
wp.init()
@wp.func
def f(n: int) -> int:
    if n <= 0:
        return 0
    return n + f(n - 1)
@wp.kernel
def k(o: wp.array(dtype=int)):
    o[0] = f(5)
try:
    o = wp.zeros(1, dtype=int); wp.launch(k, dim=1, inputs=[o], device="cpu"); print("recursion result", o.numpy())
except Exception as e:
    print("recursion error:", type(e).__name__, str(e)[:200])
@wp.struct
class Particle:
    pos: wp.vec3
    dir: wp.vec3
    e: float
    species: wp.int32
@wp.kernel
def ks(ps: wp.array(dtype=Particle)):
    i = wp.tid(); p = ps[i]; p.pos = p.pos + p.dir * p.e; ps[i] = p
ps = wp.zeros(4, dtype=Particle); wp.launch(ks, dim=4, inputs=[ps], device="cpu"); print("array-of-struct with vec3 members OK; sizeof bytes:", wp.types.type_size_in_bytes(Particle))
