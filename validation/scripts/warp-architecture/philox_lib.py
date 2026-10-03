import warp as wp
def make_philox(real):
    M0 = wp.constant(wp.uint64(0xD2511F53)); M1 = wp.constant(wp.uint64(0xCD9E8D57))
    W0 = wp.constant(wp.uint32(0x9E3779B9)); W1 = wp.constant(wp.uint32(0xBB67AE85))
    MASK = wp.constant(wp.uint64(0xFFFFFFFF))
    @wp.func
    def philox4x32_10(c: wp.vec4ui, k: wp.vec2ui) -> wp.vec4ui:
        c0 = c[0]; c1 = c[1]; c2 = c[2]; c3 = c[3]; k0 = k[0]; k1 = k[1]
        for r in range(10):
            p0 = M0 * wp.uint64(c0); p1 = M1 * wp.uint64(c2)
            n0 = wp.uint32(p1 >> wp.uint64(32)) ^ c1 ^ k0; n2 = wp.uint32(p0 >> wp.uint64(32)) ^ c3 ^ k1
            c0 = n0; c1 = wp.uint32(p1 & MASK); c2 = n2; c3 = wp.uint32(p0 & MASK)
            k0 = k0 + W0; k1 = k1 + W1
        return wp.vec4ui(c0, c1, c2, c3)
    @wp.func
    def u01(x: wp.uint32) -> real:
        return (real(x >> wp.uint32(8)) + real(0.5)) * real(1.0 / 16777216.0)
    return philox4x32_10, u01
