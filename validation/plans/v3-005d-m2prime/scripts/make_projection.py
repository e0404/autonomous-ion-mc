# Builds ../M2PRIME-projected.md from m2prime.out / m2prime2.out / ../M2PRIME.md (scratch tool).
import re, sys
PIN={'C-12':{100:222.0,150:211.0,200:210.0},'O-16':{100:293.0,150:279.0,200:276.0}}  # Barashenkov p_in (ionmc.nuclear.bgg); not printed
rx=re.compile(r"^(C-12|O-16) (\d+) (.+?)\s+sTot\s+([\d.]+) \(([+-][\d.]+)%\) sR\s+([\d.]+) sEl\s+([\d.]+) el/Bar ([\d.]+) qmin (\S+) P\S*\s+([\d.]+)")
rows=[]
for fn,src in (('m2prime.out','m2prime.py'),('m2prime2.out','m2prime2.py')):
    for ln in open(fn):
        m=rx.match(ln.strip())
        if not m: continue
        t,E,lab,st,dev,sr,se,eb,qm,p=m.groups(); E=int(E); lab=lab.strip()
        if re.match(r"\((c|d|e)\)",lab):
            rows.append(f"| {t} | {E} | {lab} | fitted (n/a) | withheld | withheld | {qm} | {p} | {src} |")
        else:
            pf='pass' if abs(float(dev))<=10 else 'fail'
            is_m2=lab in ('(a) M2','M2')
            srr='withheld (M2 reference)' if is_m2 else f"{float(sr)/PIN[t][E]:.3f}"
            rows.append(f"| {t} | {E} | {lab} | {pf} | {srr} | {eb} | {qm} | {p} | {src} |")
fac=[ln.strip() for ln in open('m2prime2.out') if 'factor at peak' in ln]
txt=open('../M2PRIME.md').read()
quotes=txt.split("## Reading of rho1, rho2")[1].split("E_F(0.17)")[0].strip()
mc=txt.split("## R-C verification")[1].strip()
out=f"""# M2PRIME - Git-safe projection (research agent, 2026-10-18)

Projection of the scratch-only M2PRIME.md. No absolute cross section (sigma_tot, sigma_R, sigma_el) of any construction and no
measured value is stated. Measured n+A total cross sections are read at run time from the acquisition cache (EXFOR entries
13569 and 13753) by verify/exfor_tot.py; each model sigma_tot is reported only as pass/fail against +-10 % of the EXFOR mean
(mean over the entries of the per-entry mean of the n+natC / n+natO (N,TOT) points within +-3 % of E). For the constructions
fitted to the measured sigma_tot ((c) alpha scale, (d) chi scale, (e) Pauli + alpha scale) the ratios are withheld too, and
for the M2 reference (a) sigma_R/Barashenkov is withheld (its sigma_tot deviations are of record in Amendment 20; with the
committed Barashenkov arrays sigma_R + sigma_el would recover the measured sigma_tot).
Residual: the ratios times the public Barashenkov values (Geant4 G4BarashenkovData) reconstruct the model sigma_tot of the
non-fitted rows; with the pass/fail flags this does not give the measured value.

Scripts (scratch): verify/m2prime.py, verify/m2prime2.py (both import verify/exfor_tot.py), verify/cv_mc.py (no data input),
verify/make_projection.py (this file).
Inputs: GWU SAID SM16 sigma_NN, alpha_NN; beta = (1+alpha^2) sigma/(16 pi) per channel; de Vries HO charge densities C-12
(a 1.692 fm, alpha 1.082), O-16 (a 1.833 fm, alpha 1.544), proton size unfolded analytically (r_p 0.8414 fm); Barashenkov
sigma_el and p_in from ionmc.nuclear.bgg; theta_2mm lab 13.033 / 9.108 / 7.098 deg at 100 / 150 / 200 MeV (D-E0 kernels, D-E5).
Constructions: (a) M2 = free-NN optical limit; (b) = R-B; (b2) Bertulani-De Conti Eq.(14) with rho~ = 2 rho/rho0;
R-A Eq.(14) with rho1 = 0, rho2 = rho_tot; R-B Eq.(14) with rho1 = rho_p, rho2 = rho_n (report-only sensitivity);
R-C Clementel-Villi single Fermi sphere, local density, bare lab E (adopted M2' definition).

## Results

| target | E MeV | construction | sigma_tot within +-10 % of EXFOR mean | sigma_R / Barashenkov p_in | sigma_el / Barashenkov sigma_el | q_min fm^-1 | P(theta_lab > theta_2mm) | script |
|---|---|---|---|---|---|---|---|---|
""" + "\n".join(rows) + """

Pauli factor (sigma_eff/sigma_free) at peak point density (m2prime2.py):
""" + "\n".join("    "+f for f in fac) + """
Geant4 CHIPS P(theta_lab > theta_2mm), O-16, D-E4 cases 3/4 (information only, not like-for-like): 0.513 / 0.571 (150 / 200 MeV).
q_H (Helm, q^4-cumulant s): C-12 1.912, O-16 1.661 fm^-1.

## Reading of rho1, rho2 - quoted sentences (Bertulani & De Conti, PRC 81 064603, arXiv:1004.2096)
""" + quotes + """

## R-C verification
""" + mc + "\n"
leak=re.compile(r"(?<![\d.])(49[0-4]|34[0-5]|29[0-8]|64[0-9]|44[0-4]|37[4-6])(\.\d+)?(?![\d])")
hits=leak.findall(out)
if hits or re.search(r"[+-]\d+\.\d%",out) or re.search(r"\bsTot\b|\bsR\b|\bsEl\b",out):
    sys.exit(f"leak check failed: {hits}")
m2rows=[r for r in rows if r.split('|')[3].strip() in ('(a) M2','M2')]
if len(m2rows)!=12 or any('withheld' not in r.split('|')[5] for r in m2rows):
    sys.exit("leak check failed: M2 sigma_R ratio not withheld")
open('../M2PRIME-projected.md','w').write(out); print(len(rows),"rows; leak check passed")
