# M2PRIME - Git-safe projection (research agent, 2026-10-18)

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
| C-12 | 100 | (a) M2 | pass | withheld (M2 reference) | 0.905 | nan | 0.4853 | m2prime.py |
| C-12 | 100 | (b) BDC Pauli LDA rho~=rho/rho0 | fail | 1.160 | 0.717 | 1.507 | 0.4750 | m2prime.py |
| C-12 | 100 | (b2) BDC rho~=2rho/rho0 [alt reading] | fail | 1.064 | 0.555 | 1.426 | 0.4592 | m2prime.py |
| C-12 | 100 | (c) alpha x 0.893 | fitted (n/a) | withheld | withheld | 1.479 | 0.4832 | m2prime.py |
| C-12 | 100 | (d) chi x 0.935 | fitted (n/a) | withheld | withheld | nan | 0.4891 | m2prime.py |
| C-12 | 100 | (e) BDC + alpha x 1.284 | fitted (n/a) | withheld | withheld | nan | 0.4798 | m2prime.py |
| C-12 | 150 | (a) M2 | fail | withheld (M2 reference) | 1.336 | nan | 0.6080 | m2prime.py |
| C-12 | 150 | (b) BDC Pauli LDA rho~=rho/rho0 | pass | 1.051 | 1.099 | 1.565 | 0.6011 | m2prime.py |
| C-12 | 150 | (b2) BDC rho~=2rho/rho0 [alt reading] | pass | 0.996 | 0.937 | 1.529 | 0.5939 | m2prime.py |
| C-12 | 150 | (c) alpha x 0.498 | fitted (n/a) | withheld | withheld | 1.535 | 0.6042 | m2prime.py |
| C-12 | 150 | (d) chi x 0.795 | fitted (n/a) | withheld | withheld | 1.651 | 0.6172 | m2prime.py |
| C-12 | 150 | (e) BDC + alpha x 0.863 | fitted (n/a) | withheld | withheld | 1.545 | 0.6001 | m2prime.py |
| C-12 | 200 | (a) M2 | fail | withheld (M2 reference) | 1.454 | 1.596 | 0.6734 | m2prime.py |
| C-12 | 200 | (b) BDC Pauli LDA rho~=rho/rho0 | pass | 0.981 | 1.235 | 1.578 | 0.6686 | m2prime.py |
| C-12 | 200 | (b2) BDC rho~=2rho/rho0 [alt reading] | pass | 0.944 | 1.096 | 1.560 | 0.6645 | m2prime.py |
| C-12 | 200 | (c) alpha x 0.386 | fitted (n/a) | withheld | withheld | 1.566 | 0.6714 | m2prime.py |
| C-12 | 200 | (d) chi x 0.828 | fitted (n/a) | withheld | withheld | 1.632 | 0.6798 | m2prime.py |
| C-12 | 200 | (e) BDC + alpha x 0.804 | fitted (n/a) | withheld | withheld | 1.566 | 0.6680 | m2prime.py |
| O-16 | 100 | (a) M2 | pass | withheld (M2 reference) | 0.851 | 1.285 | 0.4117 | m2prime.py |
| O-16 | 100 | (b) BDC Pauli LDA rho~=rho/rho0 | fail | 1.121 | 0.699 | 1.306 | 0.4013 | m2prime.py |
| O-16 | 100 | (b2) BDC rho~=2rho/rho0 [alt reading] | fail | 1.040 | 0.560 | 1.269 | 0.3865 | m2prime.py |
| O-16 | 100 | (c) alpha x 0.953 | fitted (n/a) | withheld | withheld | 1.285 | 0.4104 | m2prime.py |
| O-16 | 100 | (d) chi x 0.969 | fitted (n/a) | withheld | withheld | 1.298 | 0.4135 | m2prime.py |
| O-16 | 100 | (e) BDC + alpha x 1.294 | fitted (n/a) | withheld | withheld | nan | 0.4082 | m2prime.py |
| O-16 | 150 | (a) M2 | fail | withheld (M2 reference) | 1.407 | 1.380 | 0.5401 | m2prime.py |
| O-16 | 150 | (b) BDC Pauli LDA rho~=rho/rho0 | pass | 1.020 | 1.182 | 1.370 | 0.5331 | m2prime.py |
| O-16 | 150 | (b2) BDC rho~=2rho/rho0 [alt reading] | pass | 0.972 | 1.024 | 1.349 | 0.5261 | m2prime.py |
| O-16 | 150 | (c) alpha x 0.494 | fitted (n/a) | withheld | withheld | 1.349 | 0.5349 | m2prime.py |
| O-16 | 150 | (d) chi x 0.784 | fitted (n/a) | withheld | withheld | 1.422 | 0.5505 | m2prime.py |
| O-16 | 150 | (e) BDC + alpha x 0.827 | fitted (n/a) | withheld | withheld | 1.359 | 0.5314 | m2prime.py |
| O-16 | 200 | (a) M2 | fail | withheld (M2 reference) | 1.622 | 1.399 | 0.6115 | m2prime.py |
| O-16 | 200 | (b) BDC Pauli LDA rho~=rho/rho0 | pass | 0.961 | 1.396 | 1.387 | 0.6066 | m2prime.py |
| O-16 | 200 | (b2) BDC rho~=2rho/rho0 [alt reading] | pass | 0.928 | 1.252 | 1.381 | 0.6024 | m2prime.py |
| O-16 | 200 | (c) alpha x 0.234 | fitted (n/a) | withheld | withheld | 1.381 | 0.6084 | m2prime.py |
| O-16 | 200 | (d) chi x 0.799 | fitted (n/a) | withheld | withheld | 1.430 | 0.6199 | m2prime.py |
| O-16 | 200 | (e) BDC + alpha x 0.693 | fitted (n/a) | withheld | withheld | 1.381 | 0.6052 | m2prime.py |
| C-12 | 100 | M2 | pass | withheld (M2 reference) | 0.905 | nan | 0.4853 | m2prime2.py |
| C-12 | 100 | R-A rho1=0,rho2=rho_tot | fail | 0.582 | 0.163 | 1.607 | 0.5082 | m2prime2.py |
| C-12 | 100 | R-B rho1=rho_p,rho2=rho_n | fail | 1.160 | 0.717 | 1.507 | 0.4750 | m2prime2.py |
| C-12 | 100 | R-C Clementel-Villi LDA | fail | 1.090 | 0.595 | 1.442 | 0.4638 | m2prime2.py |
| C-12 | 150 | M2 | fail | withheld (M2 reference) | 1.336 | nan | 0.6080 | m2prime2.py |
| C-12 | 150 | R-A rho1=0,rho2=rho_tot | fail | 0.510 | 0.229 | 1.671 | 0.6274 | m2prime2.py |
| C-12 | 150 | R-B rho1=rho_p,rho2=rho_n | pass | 1.051 | 1.099 | 1.565 | 0.6011 | m2prime2.py |
| C-12 | 150 | R-C Clementel-Villi LDA | pass | 1.010 | 0.976 | 1.540 | 0.5959 | m2prime2.py |
| C-12 | 200 | M2 | fail | withheld (M2 reference) | 1.454 | 1.596 | 0.6734 | m2prime2.py |
| C-12 | 200 | R-A rho1=0,rho2=rho_tot | fail | 0.470 | 0.245 | 1.698 | 0.6912 | m2prime2.py |
| C-12 | 200 | R-B rho1=rho_p,rho2=rho_n | pass | 0.981 | 1.235 | 1.578 | 0.6686 | m2prime2.py |
| C-12 | 200 | R-C Clementel-Villi LDA | pass | 0.953 | 1.129 | 1.566 | 0.6655 | m2prime2.py |
| O-16 | 100 | M2 | pass | withheld (M2 reference) | 0.851 | 1.285 | 0.4117 | m2prime2.py |
| O-16 | 100 | R-A rho1=0,rho2=rho_tot | fail | 0.579 | 0.169 | 1.423 | 0.4359 | m2prime2.py |
| O-16 | 100 | R-B rho1=rho_p,rho2=rho_n | fail | 1.121 | 0.699 | 1.306 | 0.4013 | m2prime2.py |
| O-16 | 100 | R-C Clementel-Villi LDA | fail | 1.061 | 0.595 | 1.281 | 0.3908 | m2prime2.py |
| O-16 | 150 | M2 | fail | withheld (M2 reference) | 1.407 | 1.380 | 0.5401 | m2prime2.py |
| O-16 | 150 | R-A rho1=0,rho2=rho_tot | fail | 0.507 | 0.258 | 1.469 | 0.5623 | m2prime2.py |
| O-16 | 150 | R-B rho1=rho_p,rho2=rho_n | pass | 1.020 | 1.182 | 1.370 | 0.5331 | m2prime2.py |
| O-16 | 150 | R-C Clementel-Villi LDA | pass | 0.984 | 1.062 | 1.354 | 0.5280 | m2prime2.py |
| O-16 | 200 | M2 | fail | withheld (M2 reference) | 1.622 | 1.399 | 0.6115 | m2prime2.py |
| O-16 | 200 | R-A rho1=0,rho2=rho_tot | fail | 0.471 | 0.287 | 1.492 | 0.6324 | m2prime2.py |
| O-16 | 200 | R-B rho1=rho_p,rho2=rho_n | pass | 0.961 | 1.396 | 1.387 | 0.6066 | m2prime2.py |
| O-16 | 200 | R-C Clementel-Villi LDA | pass | 0.936 | 1.286 | 1.381 | 0.6035 | m2prime2.py |

Pauli factor (sigma_eff/sigma_free) at peak point density (m2prime2.py):
    factor at peak density 0.179: R-B 0.617  R-C 0.444  R-A 0.213
    factor at peak density 0.179: R-B 0.745  R-C 0.629  R-A 0.257
    factor at peak density 0.179: R-B 0.808  R-C 0.722  R-A 0.280
    factor at peak density 0.165: R-B 0.637  R-C 0.473  R-A 0.220
    factor at peak density 0.165: R-B 0.758  R-C 0.649  R-A 0.262
    factor at peak density 0.165: R-B 0.819  R-C 0.737  R-A 0.283
Geant4 CHIPS P(theta_lab > theta_2mm), O-16, D-E4 cases 3/4 (information only, not like-for-like): 0.513 / 0.571 (150 / 200 MeV).
q_H (Helm, q^4-cumulant s): C-12 1.912, O-16 1.661 fm^-1.

## Reading of rho1, rho2 - quoted sentences (Bertulani & De Conti, PRC 81 064603, arXiv:1004.2096)
(appended 2026-10-18; script verify/m2prime2.py)
Paper (arXiv:1004.2096, text dl_bdc/p.txt) defining sentences:
 - "The Fermi momenta kF1,F2 are related to the proton and neutron densities by means of the zero temperature density
   approximation, kFi = (3 pi^2 rho_i)^(1/3). For finite nuclei, one usually replaces rho_i by the local densities"
 - Eq.(7) context: "an expression was obtained for the geometrical Q operator for nucleon scattering in asymmetric nuclear
   matter, involving two Fermi momentum spheres, one for the proton and another for the neutron" (Ref.[21]); Eq.(7) averages
   BOTH colliding nucleons over Fermi spheres displaced by k0 (bound-bound geometry, Hussein et al. Fig.30).
 - Eq.(14): "where E is the laboratory energy in MeV, rho~ = (rho1 + rho2)/rho0, with rho0 = 0.17 fm^-3"; fit to isotropic sigma.
 - Application: "we have obtained an effective local density for protons and neutrons by taking the point along the impact
   parameter direction where the two densities (one from the projectile and the other from the target) cross each other.
   This effective density was then used in Eqs. (14), (15), and (16)"
 - Nucleon-nucleus: "Clementel and Villi [20] ... obtained an analytical expression for the scattering of a nucleon on a nucleon
   Fermi gas. By using the local density approximation, their work have been widely used to describe Pauli-blocking in
   nucleon-nucleus scattering." The paper does not treat the free-nucleon-on-matter limit with Eq.(14).

## R-C verification
(appended 2026-10-18; script verify/cv_mc.py)
Formula: x = E_F/E; f = 1 - 7x/5 (E >= 2E_F); f = 1 - 7x/5 + (2x/5)(2 - E/E_F)^(5/2) (E_F < E < 2E_F); f = 0 (E <= E_F).
No open text source quoting the equation was obtained (HRB 1991 / Kikuchi-Kawai not acquired; Kelly nucl-th/9607001 and
Chen-Sammarruca-Bertulani arXiv:1304.6096 checked, formula not present). Verified instead by independent Monte Carlo derivation
(4e6 samples/point): nonrelativistic, constant isotropic sigma, target momenta uniform in the Fermi sphere, both final momenta
> k_F, flux-weighted f = <(|k0-k2|/k0) Q>:
 E/E_F 1.00 1.05 1.20 1.50 1.80 2.00 2.50 3.00 4.00 6.00
 form  0.0000 0.0018 0.0241 0.1138 0.2262 0.3000 0.4400 0.5333 0.6500 0.7667
 MC    0.0000 0.0018 0.0241 0.1137 0.2260 0.3004 0.4399 0.5336 0.6498 0.7668
 (unweighted <Q> differs: 0.2711 at E=2E_F -> the formula is the flux-weighted, per-path-length definition, which is the
 one appropriate for the Glauber attenuation rho*sigma per unit length.)
Continuity at E = 2E_F: both branches 0.3000. Vanishes at E = E_F (f ~ (E-E_F)^n, 0.0018 at 1.05 E_F).
E_F convention: E_F = (hbar c)^2 k_F^2/(2 m_p c^2), k_F = (3 pi^2 rho_tot/2)^(1/3), rho_tot = rho_p + rho_n (point, fm^-3),
 hbar c = 197.3269804 MeV fm, m_p = 938.27208816 MeV: E_F(0.17) = 38.392 MeV, k_F = 1.3602 fm^-1 (m_N average 938.9187 -> 38.366).
Low branch on the pinned table grid (1fd24cff, 400 nodes 1-250 MeV): interior peak point density C 0.179 / O 0.165 fm^-3 ->
 2E_F 79.5 / 75.2 MeV, E_F 39.7 / 37.6 MeV. Low branch exercised somewhere in the nucleus for E < ~75-80 MeV (C 122, O 120 nodes
 in [5 MeV, 2E_F)); interior fully blocked (f=0) for E <= ~38-40 MeV (C 95, O 91 nodes in [5, E_F]). High branch only for
 E >= 2E_F(peak): C 43, O 45 nodes. All validation energies 100/150/200 MeV use the high branch only.
