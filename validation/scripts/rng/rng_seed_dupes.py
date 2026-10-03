# Exact duplicate histories between two seeds: rand_init(s,i)=pcg(s+pcg(i)); duplicates iff s1+pcg(i)==s2+pcg(j).
import numpy as np
def pcg(x):
    x = x.astype(np.uint64)
    b = (x * 747796405 + 2891336453) & 0xFFFFFFFF
    c = (((b >> ((b >> 28) + 4)) ^ b) * 277803737) & 0xFFFFFFFF
    return ((c >> 22) ^ c) & 0xFFFFFFFF
N = 1_000_000
p = pcg(np.arange(N, dtype=np.uint64))
for s1, s2 in [(1, 2), (1000, 1001), (42, 43)]:
    a = (p + s1) & 0xFFFFFFFF; b = (p + s2) & 0xFFFFFFFF
    print(f"seeds {s1},{s2}: identical start states between batches (N={N}):", np.intersect1d(a, b).size)
# within one seed: are start states distinct?
print("distinct pcg(i) for i<1e6:", np.unique(p).size)
