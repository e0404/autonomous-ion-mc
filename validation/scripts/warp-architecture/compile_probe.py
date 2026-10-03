import sys, time, warp as wp
eb = sys.argv[1] == "1"
wp.config.enable_backward = eb
wp.init()
from toy_transport import build
k = build(wp.float32)
t = time.perf_counter(); wp.load_module(k.module, device="cpu"); print(f"enable_backward={eb}: CPU cold compile {time.perf_counter()-t:.2f}s")
