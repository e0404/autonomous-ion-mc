"""Operating characteristics of candidate V7-R gate designs (binomial model; paired gate exact given the true
coverage, single gate approximated as CP on n_E i.i.d. hits without the reference box, which only makes the real
gate stricter). kappa error e = realised true coverage minus target at the calibration margins."""
import math
from binom_gate import accept_set, p_accept, p_le, LO, HI, NOM

def oc(n, a):
    hmin, hmax = accept_set(n, a)
    f = lambda c: p_accept(n, hmin, hmax, c)
    return f

def cal_pass_prob(rate, nmax, ncal=3200):  # P(observed <= nmax) at the true rate
    return p_le(ncal, nmax, rate)

designs = [  # name, pairs, single-eval replicates, alpha
    ("D0 as frozen", 3600, 5400, 0.04),
    ("D1 alpha 0.025", 3600, 5400, 0.025),
    ("D2 alpha 0.02", 3600, 5400, 0.02),
    ("D3 2x replicates, alpha 0.02", 7200, 10800, 0.02),
    ("D4 2x replicates, alpha 0.015", 7200, 10800, 0.015),
]
print("design | paired FA@margin, power | single FA@margin, power | joint power (indep.) | P(cal clauses pass) for kappa error e = 0 / 0.0007 / 0.0034")
for name, npair, nE, a in designs:
    fp, fs = oc(npair, a), oc(nE, a)
    pfa, sfa = max(fp(LO), fp(HI)), max(fs(LO), fs(HI))
    pp, sp = fp(NOM), fs(NOM)
    jp = pp * sp
    outs = []
    for e in (0.0, 0.0007, 0.0034):
        pfe, sfe = fp(LO + e), fs(HI - e)
        # gate-alone clause: <=109/3200; joint clause: <=50/3200 (joint <= min(gate) rate * conditional; bound: use 0.5*min as planning value)
        p_alone = cal_pass_prob(pfe, 109) * cal_pass_prob(sfe, 109)
        p_joint_ub = cal_pass_prob(min(pfe, sfe), 50)  # pessimistic: joint as large as the smaller gate rate
        p_power = 1 - p_le(3200, 2908, jp)  # >= 2909 accepted at the nominal
        outs.append(f"{p_alone * p_power:.3f}(joint>={p_joint_ub:.2f})")
    print(f"{name:32s} | {pfa:.4f}, {pp:.4f} | {sfa:.4f}, {sp:.4f} | {jp:.4f} | " + " / ".join(outs))
for nk in (20000, 200000, 1000000):
    print(f"kappa-solve n={nk}: coverage SE {math.sqrt(LO*(1-LO)/nk):.5f}")
