import warp as wp, numpy as np
wp.init()
@wp.func
def f(x: float) -> float:
    return x + 1.0e-10
@wp.func
def g(x: wp.float64) -> wp.float64:
    return x + wp.float64(1.0e-10)
@wp.func
def h(x: float) -> float:
    return wp.log(x)
r = f(1.0); print("f(float) py-scope:", repr(r), type(r), "==1.0?", r == 1.0)
r = g(wp.float64(1.0)); print("g(float64) py-scope:", repr(r), type(r))
try:
    r = g(1.0); print("g(python float) :", repr(r), type(r))
except Exception as e: print("g(python float) error:", e)
r = h(3.0); print("h log(3) py-scope:", repr(r), "f32 log:", np.log(np.float32(3)), "f64 log:", np.log(3.0))
@wp.func
def v(a: wp.vec3) -> wp.vec3:
    return wp.normalize(a)
print("vec3 py-scope:", v(wp.vec3(1.0, 2.0, 2.0)), type(v(wp.vec3(1.0,2.0,2.0))))
