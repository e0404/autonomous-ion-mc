import warp as wp, numpy as np, time
wp.config.log_level = 40
wp.init()
F = wp.float64; V = wp.types.vector(3, wp.float64)
@wp.func
def h(x: F) -> F:
    return wp.log(x) + wp.sqrt(x) * wp.cos(x)
@wp.func
def v(a: V) -> V:
    return wp.normalize(a)
r = h(F(3.0)); print("float64 py-scope:", repr(r), "numpy f64:", np.log(3.0)+np.sqrt(3.0)*np.cos(3.0))
print("vec3d py-scope:", v(V(1.0, 2.0, 2.0)))
@wp.func
def h32(x: float) -> float:
    return wp.log(x) + wp.sqrt(x) * wp.cos(x)
t=time.perf_counter()
for i in range(20000): h(F(3.0))
t1=time.perf_counter()
for i in range(20000): h32(3.0)
t2=time.perf_counter()
print(f"py-scope call cost: float64-typed {1e6*(t1-t)/20000:.1f} us, float {1e6*(t2-t1)/20000:.1f} us")
