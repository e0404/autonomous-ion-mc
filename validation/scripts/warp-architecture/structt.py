import warp as wp
wp.init()
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
