"""EXPLORATORY (planner, Revision 3): the section 3.8 selection rule applied to the measured kappa=1 coverages.
P(power clause) = P(Bin(3200, pw) >= 2909), pw = OC_paired(c_p) * OC_single(c_s, box shrink 0.0019*n_E)."""
import sys; sys.path.insert(0, sys.argv[1] + "/v7r-scripts")
from fastbin import oc, cdf, kmin_pow
kp = kmin_pow(3200)
for name, cp, cs in (("archive_fitted", 0.68396, 0.68664), ("rehearsal_block_sums", 0.68305, 0.68083)):
    for npair in (5400, 7200):
        nE = npair * 3 // 2
        for sh in (int(round(0.0019 * nE)), 0):
            pw = oc(npair, 0.005, cp) * oc(nE, 0.005, cs, sh)
            print(f"{name} {npair}/0.005 shrink {sh}: power {pw:.4f}  P(power clause) {1 - cdf(3200, kp - 1, pw):.4f}")
