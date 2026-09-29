# Ion-specific stopping-power and range data sources (research note, 2026-09-30)

Source: physics-researcher subagent. NIST endpoint and G4EMLOW contents were
checked first-hand; items marked **unverified** were not.

## NIST PSTAR/ASTAR (checked by live query)
- `POST https://physics.nist.gov/cgi-bin/Star/ap_table.pl` (multipart form):
  `prog=PSTAR|ASTAR`, `matno=276` (liquid water; 277 water vapour; elements
  001 H, 006 C, 007 N, 008 O; tissues 103 adipose, 119 compact bone, 120
  cortical bone, ...), `ShowDefault=on`, `GraphType=None`, optional `Energies`.
- Output: HTML table — T (MeV); electronic/nuclear/total stopping power
  (MeV cm²/g); CSDA and projected range (g/cm²); detour factor. The HTML is
  not byte-stable (scripts); hash the parsed table.
- PSTAR: 1 keV–10 GeV. ASTAR: T is the **total** alpha energy, 1 keV–1000 MeV
  (≤ 250 MeV/u); cannot cover He above 250 MeV/u.
- Spot checks: PSTAR water 100 MeV electronic 7.286 MeV cm²/g, CSDA 7.718 g/cm²;
  200 MeV CSDA 25.96 g/cm²; ASTAR water 100 MeV electronic 86.45, CSDA 0.6409.
- Evaluation: ICRU 49 (1993), water I = 75 eV. Access: NIST SRD 124, freely
  accessible; record the access basis, do not redistribute (unverified terms).
- ICRU 90 addendum "Update ESTAR, PSTAR and ASTAR databases":
  https://www.nist.gov/document/update-estar-pstar-and-astar-databases
  (17-page PDF, 1,589,504 bytes, sha256 e8340e3d…b474); location of the
  machine-readable ICRU 90 tables **unverified**.

## Geant4 G4EMLOW 8.8 (downloaded and inspected, then discarded)
- https://cern.ch/geant4-data/datasets/G4EMLOW.8.8.tar.gz (350,393,595 bytes,
  MD5 328330009df633f7e9b3a9f445745298 matching Geant4 v11.4.2 CMake, sha256
  b60cfd63176f5d16107e2a25b35b235155032d1735d749670ca50fede12624cf).
- `ion_stopping_data/icru73` (Z = 3–80; Z = 3–18 in 70 targets incl.
  G4_WATER and tissues; 0.025–1000 MeV/u) and `icru90` (Z = 3–18 in water,
  air, graphite, H/C/N/O; 0.00025–1000 MeV/u); units MeV cm²/mg; both are PASS
  (Schinner/Sigmund binary theory) output, not verbatim ICRU tables; agree
  within ≈0.05% above 10 MeV/u for C/O, 3–7% apart below 1 MeV/u.
- **Licence restriction:** README states the data are not for commercial use
  and must be used within Geant4. Decision: not used for construction or
  evaluation in this project; the local scratch copy was deleted.
- Geant4 uses these tables only below 2.5 MeV/u (`G4IonICRU73Data.cc`).
- ICRU 90 proton/alpha arrays for air/water/graphite are hard-coded in
  `source/materials/src/G4ICRU90StoppingData.cc` (Geant4 Software License):
  water protons 100 MeV 7.250 MeV cm²/g (0.5% below PSTAR); alphas 100 MeV 85.93.

## Other sources (unverified)
- IAEA-NDS "Electronic Stopping Power of Matter for Ions" (Paul/Schinner)
  measured compilation, https://www-nds.iaea.org/stopping/ — the main
  **measured** He/C/O data, mostly below ≈10 MeV/u, few liquid-water cases.
- MSTAR (Paul & Schinner) scales ASTAR to Z = 3–18; libdedx bundles tables
  (licence unchecked); SRIM freeware, semi-empirical fit to the Paul data.
- arXiv:1812.07877 on ICRU 90 impact for carbon-ion stopping-power ratios.

## Mean excitation energy of water
ICRU 49/PSTAR: 75 eV. ICRU 90 (2016): 78 eV (graphite 81, air 85.7).
Effect ≈ −0.4…−0.5% stopping, +0.4…+0.5% range (≈ +1–1.5 mm at 25–30 cm),
similar relative shift for carbon at equal velocity.

## Recommendation summary
Construction: PSTAR/ASTAR (ICRU 49) for p and He ≤ 250 MeV/u, or ICRU 90 p/α
arrays; for C, O and He > 250 MeV/u, no licence-clean open table exists —
build Bethe-based tables (shell, Barkas, Bloch, Lindhard–Sørensen, effective
charge) with a single consistent I; these are a *related model* and need
ion-specific held-out evidence (measured ranges/Bragg peaks, Paul compilation).
Shared lineage: PSTAR/ASTAR ≡ Geant4 PSTAR/ASTAR arrays; ICRU 73/90 heavy-ion
tables ≡ PASS; ICRU 90 p/α share BEST/Bethe formalism with PSTAR.
