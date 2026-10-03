"""Cold CPU compile time of the toy transport module with a fresh, empty kernel cache."""
import os, sys, tempfile, time
cache = tempfile.mkdtemp(prefix="ionmc-warp-cold-")
os.environ["WARP_CACHE_PATH"] = cache  # must be set before importing warp
import warp as wp
eb = sys.argv[1] == "1"
wp.config.enable_backward = eb
wp.config.log_level = wp.LOG_WARNING
wp.init()
kdir = wp.config.kernel_cache_dir  # warp appends a version subdirectory below WARP_CACHE_PATH
assert kdir.startswith(cache), f"unexpected kernel cache {kdir}"
assert not any(os.scandir(kdir)) if os.path.isdir(kdir) else True, "cache not fresh"
from toy_transport import build
k = build(wp.float32)
t = time.perf_counter(); wp.load_module(k.module, device="cpu")
print(f"enable_backward={eb}: CPU cold compile {time.perf_counter() - t:.2f}s (fresh cache {cache}, warp {wp.__version__})")
