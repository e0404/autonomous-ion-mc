# Operator directive (2026-11-07): validation target for V3-008 changed to agreement with TOPAS

Received 2026-11-07 (UTC) as the operator's response to intervention IR-20261102-010435-B59B2D (task V3-008).
Recorded verbatim by the orchestrator; processed as a version-controlled amendment of the V3-008 acceptance plan
(revision 2) and decision 0043. The text between the rules is the operator's message, unchanged.

---

For IR-20261102-010435-B59B2D, I am changing the validation target for this stage. Do not pursue experimental measurement acquisition as a prerequisite. Use agreement with TOPAS as the acceptance reference for depth-dose distributions and lateral profiles.
The final TOPAS reference configuration should use:
sv:Ph/Default/Modules = 7
    "g4em-standard_opt4"
    "g4h-phy_QGSP_BIC_HP"
    "g4decay"
    "g4h-elastic_HP"
    "g4stopping"
    "g4ion-QMD"
    "g4radioactivedecay"
Autonomously choose and justify representative validation cases, materials, energies, observables, statistical precision and acceptance criteria. You may configure simpler TOPAS physics for intermediate validation, including electromagnetic-only comparisons. Retain useful complementary checks such as libdedx comparisons.
The proton-equivalent extrapolation is your responsibility as a model-design decision. Assess its adequacy through the resulting ion predictions against the chosen reference; do not make obtaining high-energy proton measurements an operator dependency.
Record this instruction as an explicit, version-controlled amendment and process the intervention response accordingly. Experimental validation is deferred; describe achieved validation as agreement with the specified TOPAS configuration. Preserve existing results and their historical verdicts. Continue autonomous work under the amended target.

---

Orchestrator reading (2026-11-07): (1) the intervention's unavailable inputs (measured lateral data, independent-physics
single-scattering reference, measured C/O ranges) are no longer prerequisites: experimental validation is deferred;
(2) the acceptance reference for depth-dose distributions and lateral profiles is agreement with OpenTOPAS 4.3 /
Geant4 11.4.2 in the seven-module configuration above (the "specified TOPAS configuration"); EM-only or otherwise
simplified TOPAS configurations may be used for intermediate validation; (3) the case set, materials, energies,
observables, statistical precision and acceptance criteria are chosen and justified autonomously in a reviewed plan
revision before any result; (4) complementary checks (libdEdx I3, Molière–Hanson I11, V5-SS) are retained where
useful; (5) the proton-equivalent (f_dM) extrapolation is a model-design decision assessed through the ion
predictions against the TOPAS reference, not an operator dependency; (6) every existing result and historical verdict
(P1 Pb FAILED then cleared by attribution, P6 EVIDENCE GAPs, V3 EVIDENCE GAP under revision 1, the SEEN rows) is
preserved as history; achieved validation is described as agreement with the specified TOPAS configuration.
