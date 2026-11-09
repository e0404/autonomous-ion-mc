# P9: freeze record of the V3-008 acceptance plan revision 1 and decision 0043 (D-E8 pattern)

Prepared 2026-11-01 by the planner at task head `f694b747be015abad0602f1d886f8a91cc8a604e` (clean, pushed,
Codex-passed). This record binds the state against which revision 1 of
`validation/plans/v3-008-acceptance.md` and decision 0043 are frozen. It contains hashes, identifiers, row
lists and rules only: no stopping, range or scattering reference value, and no transport or TOPAS output.

## 0. Freeze conditions (read first)

The plan freezes only "after every pre-freeze check of section P is closed" (plan header, §P order, §F gate 1).
The audit in plan §P "Freeze audit (revision 1)" finds two items that the plan lists as required before freeze and
that have no closing record at `f694b747`:

| Id | Item | Plan text that requires it | State at `f694b747` | What disposes it |
|---|---|---|---|---|
| FC-1 | **P1 completion (order step 5):** Pb attribution on ICRU 49 ranges; a water/compound θ_Hanson check if a source exists | P1 row ("P1 completion (order step 5)"); §P order step 5; §F gate 1 | not executed; no record (`p1-moliere-hanson-reference.md` "OPEN (after P1b)" items 2 and 7 unchanged) | (a) execute it and add a dated P1 addendum, or (b) a dated orchestrator amendment making it optional, recorded "not done"; either way P1 stays FAILED for Pb unless cleared |
| FC-2 | **P2 follow-up residuals and the P4 material printout:** pin `G4BetheBlochModel.cc/.hh` and settle its ICRU 90 branch; X0 (RadL) of G4_WATER; λ charge scaling and lateral displacement for ions; the `G4VMscModel.hh` initialisers; the icru90 z6/z8 silent-fallback check; the G4_WATER I-value and RadL from the archived smoke printout | P2 row ("Follow-ups, required before freeze"); P4 row (r5 sentence: the smokes archive "the G4_WATER material table (RadL)"); P5 row ("Oxygen/carbon lineage checks (P2 follow-up)") | partly done (F0–F5 of `p2-geant4-ion-models.md`; `G4BetheBlochModel.cc` hashed but not acquired and its ICRU 90 branch "not settled", `p3b-ls-default.md` §2 and §5); the rest has no record | (a) bounded source reads / a read of the archived smoke logs (no new run), recorded in a dated P2/P4 addendum, or (b) a dated orchestrator amendment declaring the residuals as lineage items (no row depends on them; see plan §P audit) |

**Dated note (2026-11-01, orchestrator) — FC-1 DISPOSED.** Record `validation/plans/v3-008/p1c-completion.md` (sha256 a1a24cd34c664a5dbb74b732bebff9244a55b4c12bfa30f0bff9c0cb8ebed0a5; computed in the scratchpad from the committed reference `moliere_hanson.py` 5a99b3b8… unchanged, variants as runtime monkeypatches only). Outcome: **Pb CLEARED by attribution** — the P1b Pb failure (+1.10/+1.32 % at x/R1 ≥ 0.9) is reproduced by the range-table input (the ionmc Bethe Pb table is +0.85 % longer than Gottschalk's ICRU 49 R1 and +2.6–3.1 % in residual range at 20–41 MeV, amplified by d ln θ/d ln x = 2.3–5.7 near the end of range); on an ICRU 49 node path the committed reference reproduces Gottschalk's Pb θ_Hanson to ≤ 0.044 % (≤ 0.024 % on a PSTAR path), controls Be/Al/Cu ≤ 0.105 %; the Fano form (Δb 0.05–0.06 for Pb), the compound rule (n/a) and the quadrature (≤ 1.2e-6) do not explain the failure. The committed test (ionmc Bethe ranges) still fails at x/R1 ≥ 0.9, so its xfail(strict=True) stays as committed (no bound file changes before the freeze); P1 is CLOSED with the attribution and the proton validation statement reads 'passes Be/Al/Cu; Pb cleared by attribution to the range-table input'. Water/compound check: Gottschalk 1993 (NIM B 74, 467) is paywalled and lists no water target; arXiv:1610.01279 tabulates measured θ0 and θ_Hanson for polystyrene, Lexan, Lucite, Nylon, Teflon and C (no uncertainties) — the committed reference agrees with Gottschalk's θ_Hanson (LOOKUP) within −0.57 to −0.04 % (polystyrene −0.16/−0.20 %, Lexan C16H14O3 −0.15/−0.10 %), which closes OPEN item 7 for H/C/N/O/F compounds and discriminates against the electron- and mass-weighted Fano averages (+4.1 % / +1.8 %); **water remains not validated directly** (the only water θ_Hanson value available has two digits; used as a check on b: ours low by 0.040 ± 0.002, which moves θ by ≤ 0.16 % at I11-like B). Side finding (outside the V3-008 domain, recorded for high-Z use): ionmc Bethe Pb ranges are 0.8–3 % longer than ICRU 49 at 20–160 MeV. TO-VERIFY: Gottschalk 1993's measurement uncertainties; the exact LOOKUP compositions. The P1 proton tests skip unless IONMC_CACHE_DIR names the objects store (test default cache path; recorded, not changed).

**Dated note (2026-11-02, orchestrator) — intervention filed.** IR-20261102-010435-B59B2D (inaccessible-information) filed at 655d16a0 and committed verbatim at cc2e4e26 (`reference_requests/interventions/IR-20261102-010435-B59B2D/`, README.md 86eba08d…, request.json 068cafd7…), as the plan section "Intervention at the freeze" requires (Codex REVIEW-100e7fb6 finding 1). FC-1's calculation is committed with provenance in the same freeze commit (finding 2; `validation/scripts/ion_em/p1c_attribution.py`, `tests/ionmc/test_ion_em_p1c.py`, inputs pinned: arXiv:1610.01279 PDF 934fd316…, NIST PSTAR Be 0fcbb24d…, Al 3fe6459e…, Cu 986e8577…, Pb 11dbcd1b…, matno 219 8974efe7…, 223 1e13828f…, 226 dc65ded9…; evaluation role). The §1 diff rule is amended: the freeze commit may add the FC-1 calculation module and test and the intervention record, and may change `src/ionmc/physics/stopping.py` docstrings only (Barkas-form attribution, Codex finding 3); the §2 hashes are regenerated at the freeze commit's parent by the orchestrator's script and listed in a dated §2 note.

**Rule.** The freeze commit is made only after FC-1 and FC-2 are disposed. Any record or amendment that disposes them
is committed before or in the freeze commit; its sha256 is added to §2 below by a dated note, and §5 is re-checked.
Until then revision 1 is a freeze *candidate*, and nothing below establishes a frozen state.

**Dated note 2026-11-01 (planner, at 6bdf3897): FC-2 DISPOSED** by the executed record
`validation/plans/v3-008/p2b-follow-ups.md` (sha256 `9c70eb1b9d49f0375af7976232817dbda925576eeb4b64e31aee61b33dfb5457`, computed after pre-commit), option (a)
of the FC-2 row: bounded source reads (`G4BetheBlochModel.cc/.hh`, `G4NistMaterialBuilder.cc`, the pinned P2/P10
sources) and a read of the archived smoke logs, no new run. Its sha256 is added to §2.4 and the new sources to §3
by dated notes; the plan carries "Amendment (2026-11-01) — FC-2 disposition" (Amendments item 6) with pointer notes
in rows P2, P4, P5 and the freeze audit, and decision 0043 §6 the matching lineage amendment. Residuals declared as
lineage items, none a freeze condition: the `G4ICRU90StoppingData` index/span (INFERRED), the possible z = 1 double
count of ΔL_LS, the Urban λ charge scaling and `dtrl` for ions (inline headers not pinned), and the absent read
evidence for `icru90/z{6,8}_G4_WATER.dat`. This note concerns FC-2 only.

## 1. Commit lineage

**Preparation head:** `f694b747be015abad0602f1d886f8a91cc8a604e`.

**Required ancestors of the freeze commit** (`git merge-base --is-ancestor <c> <freeze>` must succeed for each):

| Commit | Full SHA | Content |
|---|---|---|
| 2632a46 | `2632a46560af7f607e590896c7727406c9aeb61c` | DRAFT r3 with decision 0043 (first scientific commit, a draft) |
| 6e155b5 | `6e155b56468c6c4a146dcf9daeebc59076a4da77` | P8, P7, P2, P1/P1b records (added by the planner to the orchestrator's list: these records are bound below) |
| 4927958 | `4927958d01be0e95b8291f023010f6443225aefd` | DRAFT r5; REVIEW-98315d4e fixes |
| 3d2547a | `3d2547adcafcdd268fdef33395e4099359b2e163` | freeze-record hashes regenerated |
| 070b07f | `070b07f07ce2280e1d24f8c25326743a8438b3ea` | P2 follow-up, P10, P11 (mapping, tests and first I11 evaluation together) |
| 00af9a6 | `00af9a6c181a5c3bd32a303c4c93cd4467452542` | I11-rel at equal velocity; P5, P6 |
| b3eef90 | `b3eef90b494d94267d155d873bb07d956a1156cd` | P5b, P6-lateral; P11b |
| 41a1f8a | `41a1f8ac53bc23d600bb0874cfee67f83cc84cbc` | DRAFT r6; P4/P10 smoke cases |
| a66867d | `a66867dedb66daa142ce25a277e268f61532bc0a` | P3 closed (first time); V5-SS branch B; orchestrator decisions |
| ce4d3d6 | `ce4d3d6a2d1b95323d30f1a9c79eb5d5c970aa48` | P4 CLOSED; 1 mm SS smoke case |
| ea220a8 | `ea220a883c74ecf3883f1e31a7f35899b827a440` | round V3-008-4 dispositions; `ss_budget.py` |
| daa1520 | `daa1520a52b745c772a3acfc53d6808252f7bbb6` | LS option (default off by amendment); V5-SS TOST, 6e4 design |
| 683480b | `683480b19f6afb18aaa49a0ac40f8246cd5916ab` | exact-Ω mapping in src; I0 fields |
| 7fac30a | `7fac30add8804c7c9b2cdcc7512be17429e3a3b9` | recorded hashes regenerated |
| f694b74 | `f694b747be015abad0602f1d886f8a91cc8a604e` | I0 implementation hashes rebound |

**Clean-tree and diff rule.**
- The freeze commit's tree is clean (`git status --porcelain` empty at commit and during review).
- `git diff --name-only f694b747 <freeze>` lists only: `validation/plans/v3-008-acceptance.md`,
  `decisions/0043-ion-em-transport.md`, `validation/plans/v3-008/review.md`, this file, and the records or amendments
  that dispose FC-1/FC-2 (each added to §2 by a dated note; FC-1: `validation/plans/v3-008/p1c-completion.md`). Any other path is a STOP for the freeze.
- Every sha256 in §2 recomputed at the freeze commit (`sha256sum <path>`) equals the value below.

## 2. Bound repository files (sha256 at `f694b747`; `sha256sum`, git-independent)

### 2.1 Molière–Hanson reference, P11 mapping script, P3 and V5-SS calculations

| File | sha256 | Bound by |
|---|---|---|
| `validation/scripts/ion_em/__init__.py` | `7962c4af19580a10f51a4bfc703ce45b67a0f96f010aaaf878aef6919b2ed4ea` | P1 |
| `validation/scripts/ion_em/moliere_hanson.py` | `5a99b3b8748a6865584fcbe1337206222ddd7fe197e95f25656eac2324acf164` | P1 (REVIEW-98315d4e table) |
| `validation/scripts/ion_em/exact_omega.py` | `f5bc5f2c684e060f550f26b1d5e621e2ef5ecf01d4e9c35ed51d326d928c41b3` | P11 |
| `validation/scripts/ion_em/lindhard_sorensen.py` | `3bfeb288df6b586756e9420b32af809e744cc783eae498547cace47625114684` | P3 |
| `validation/scripts/ion_em/p3_driver.py` | `50db856a5341826e8c915eb11a6df82db38108e2dfbb4b6ab738764c20b17129` | P3 |
| `validation/scripts/ion_em/ss_budget.py` | `13158795c32094609188b22ec64b2d612e3037a1d782b392801f5aa53af4e009` | P9 (iv), V5-SS |

### 2.2 Tests

| File | sha256 | Bound by |
|---|---|---|
| `tests/ionmc/test_ion_em_moliere_hanson.py` | `64a408298bff0a66f6d43dbebb63dadf971c32378f36c16f43aab9548146f748` | P1 |
| `tests/ionmc/test_ion_em_exact_omega.py` | `12cbdba28530db9561ac469acba5d45f7533690f5587dd74d44b2143b542baee` | P11 (I11-map/I11-rel, reference script) |
| `tests/ionmc/test_ion_em_lindhard_sorensen.py` | `d35ca4e34610e153a6aa31f672a180cab4092572c8c839dc4087e606024dc6ef` | P3 |
| `tests/ionmc/test_ion_em_ss_budget.py` | `ca49735b5b00831aeaabd59d63faecfc23b98fe0d4247024ca963d7e9a77f097` | P9 (iv) |
| `tests/ionmc/test_ion_em_ls_option.py` | `c67966c8d2c3d9564032e33e666321dd11bd8d2386780afa0cb74850455ca84d` | P9 (iii) |
| `tests/ionmc/test_ion_em_mapping_src.py` | `865d3143892f8c6ebc548a86400bae198d455c3ad8c3bec2df074ebeb0210175` | I11-map (src), I0, I1 |
| `tests/ionmc/test_v3_008_seed_families.py` | `5a11ca5128709fe23618a8179dab0288c2d01fbe0431f5b328ca45d5107e9114` | P7, §C |

### 2.3 Implementation already on the branch (src)

| File | sha256 | Bound by |
|---|---|---|
| `src/ionmc/physics/lindhard_sorensen.py` | `800d1f0a13524af5bad3946032f894479f606ba9b1482bf3ea190eeb9b4d3a6a` | P9 (iii) |
| `src/ionmc/physics/exact_omega.py` | `6471c388fb75c94e9f2f4f1d19ac94be044f73444cf1e7b4bfb1a62e53ecf57b` | I0/I11-map notes |
| `src/ionmc/physics/stopping.py` | `89db1de28e5cd1f03efc1622669b58dd766fc88d5d0a91405af190e96094753e` | P9 (iii) |
| `src/ionmc/config.py` | `a24aa6b44b72a26d2771d50ebfafa447fbc7aa141f70b5be324d67fec5dfc411` | I0 note |
| `pyproject.toml` | `0e73bd64aea5ef153e198c00a67c20efbc8ef1a5061f1f77d2bdab555dda9352` | lint exceptions for the two mapping files (added by the planner; the only other changed file outside the plan set) |

### 2.4 Pre-freeze records

| File | sha256 |
|---|---|
| `validation/plans/v3-008/p1-moliere-hanson-reference.md` | `a6c901975cbaea7252a531600cd98c5caf4987781f7325863785ffaf1a02f751` |
| `validation/plans/v3-008/p2-geant4-ion-models.md` | `4ef361abfcb51928516e5a2dded48fdd6b6af6eed6c9e35aeb8b5af951ed13b1` |
| `validation/plans/v3-008/p3-lindhard-sorensen.md` | `a83bd2046443b6cbfac1592b38fedd53256732e4123434935c8e3a50b2d0cecd` |
| `validation/plans/v3-008/p3-output.json` | `491355f267b55fb742b5963f31ebe1ad6e43e0f65256c212136e0947caa8609c` |
| `validation/plans/v3-008/p3b-ls-default.md` | `01f75993b259a3b1c2d86a3ee88b345c048c1750f3963d005ed7f2dac7614a0b` |
| `validation/plans/v3-008/p5-g4emlow-ion-stopping.md` | `1a44dd50e4e59c8b74418a6fce54b296f289673ec2631247672d2a35496ea9b9` |
| `validation/plans/v3-008/p5b-libdedx-icru73.md` | `1aea3bd74fe39b8dbc8f808eecf39a6455e09fb7d6b6da5011869c4166a4917f` |
| `validation/plans/v3-008/p6-measured-data-access.md` | `60ad27f264f50dd9bd4ff75e59a4fd479268b48f64ad09367828ef9922e6f83a` |
| `validation/plans/v3-008/p6-lateral-data-search.md` | `d149550b8dc5684ec4cd7530a88e50b75bcf083dbc385934aa73e33a708212ef` |
| `validation/plans/v3-008/p8-repository-facts.md` | `18639c6bc5ff8152c39adc971652e14fb35f318499b9821d45de03bf9b374478` |
| `validation/plans/v3-008/p10-single-scattering-availability.md` | `a7f539a0b71a39939aed2644a9d7c2969133251ceb1676a3e48572a9292f7fda` |
| `validation/plans/v3-008/p11-exact-omega-mapping.md` | `b40530f48391d6eca98656d5bdb06a5e4272bf22ee0a469e08357ebcd457e330` |
| `validation/plans/v3-008/research/research-models.md` | `d194a516967f290e9aa65f693c3c414282c7d45cde5c1226ab08464a46b396cd` |
| `validation/plans/v3-008/research/research-data-sources.md` | `2bf40c44d73dd1d5ce6afe6d05d7891128cb3b4468af2afbf196f4b99afa0732` |

**Dated note 2026-11-01 (FC-2 record, bound by the §0 rule):**

| File | sha256 |
|---|---|
| `validation/plans/v3-008/p2b-follow-ups.md` | `9c70eb1b9d49f0375af7976232817dbda925576eeb4b64e31aee61b33dfb5457` |

Not bound by hash, because the freeze commit itself changes them: the plan, decision 0043, `review.md` and this
file. They are bound by the freeze commit SHA.

**Dated note (2026-11-01, FC-1):** `validation/plans/v3-008/p1c-completion.md` `a1a24cd34c664a5dbb74b732bebff9244a55b4c12bfa30f0bff9c0cb8ebed0a5` (added by the freeze commit; §7 counts: +1 repository file).

### 2.5 P4/P10 scorer-free smoke cases (execution only; seeds from the rehearsal family 2060xxxx)

| Case directory (`validation/reference_cases/topas/`) | Seed | `case.json` sha256 | `input.txt` sha256 |
|---|---|---|---|
| `he-water-150mevu-em-smoke-alpha` | 20600401 | `a41b419e28687415d1c01558b82b6e110c85a4c0a611d0e6ec49162ead8232f8` | `89a71039a43a905a3228d46566b22e3b4d3590eef294f403161f4cf75b8cef74` |
| `he-water-150mevu-em-smoke-genericion` | 20600402 | `ce0306a6e0d9c3fcbb2762de3cbad207a68887121b8149c45e85981421ca5c3d` | `b04edbddbab28f97f2f6fd94572c331e4bfcbee8c82749a00f81669bb0073c18` |
| `carbon-water-290mevu-em-smoke` | 20600403 | `d5f217bcea57d1385b0528dfdd350c467fe72cc9abc40cf151c3f36265baebee` | `9a8048cd34dc56796fc628e6a33c0b7e4a80a059abba807f0782d87e23aca9e1` |
| `oxygen-water-350mevu-em-smoke` | 20600404 | `95afc1c9f7f320810b965d67ceed8410737aa2a9c811df87fc817a5cc5282c08` | `7d91ee72412f7d55389e2901d3e441cb8ad2649f7dc592e9540dbf7d900da28a` |
| `proton-water-150mev-em-smoke-p4` | 20600405 | `0f4bc5fcdf694d2c1c9cdcaf6f8d708a725a4886d89b09f7c9f6f1dea34a3a23` | `1e8de716eaa27b7f20303a21c92bbbf6445861d0d7f16efd10b9f56989ab6b71` |
| `carbon-water-290mevu-ss-smoke` | 20600406 | `3b028844a4122c1d6baf43964644cff3dc8230ce5c3851585550ee4728a07466` | `8f3ae1470696d78f883842712fb11a091136e3df8b49327bf9cf0e7f0834ad21` |
| `proton-water-150mev-ss-smoke` | 20600407 | `3d53f66509a03bcbb69e3b2de0f5dee5f28f0da43162966fcd4bbab73782091f` | `79b0c77bb1730afbd09f696cb5a1d862e7476af6c2824950d662d3eeb04bcd19` |
| `he-water-150mevu-em-smoke-alpha-verbose` | 20600408 | `830c4c3477baf3a6a9b9fbd6da3cb78f6b3a6d1a69fe88fe4d6f6bada6177128` | `c0099ce0b96dfa3c91271043ea54871630470aefdd4f43b083b449f133f07c39` |
| `he-water-150mevu-em-smoke-genericion-verbose` | 20600409 | `6330239e13a564a4215ac13dcc777a36eea008975faa818c2cdf1f41e8f426a6` | `77f243cbc5c45495361e1cbbb018a6cd46641deaba3a8bdb0a9856daddec8bdf` |
| `carbon-water-290mevu-ss-slab1mm-smoke` | 20600410 | `40ce2cf3bceaf60b5633f9d2e4f593418275c5bbfa6e2c84aaacf72062303ba4` | `5e5487c1dba22aaa69334bcab4e0ef362e595916fa665b6412830e0247aa486b` |

`validation/reference_cases/README.md`: `2c160fb585dc571d26753aa4a1a43e44279ae4f41dfbeba480f96452a07b491e` (revision-1 value; **superseded in revision 2**, current entry in R2.2).

**Smoke run ids** (as recorded in the plan P4/P10 rows; only the last is recorded in full in the V3-008 records):
proton opt4 `REF-188fa1d8…`; alpha `REF-ed13a238…`; GenericIon(2,4) `REF-c764466c…`; carbon `REF-5910656d…`;
oxygen `REF-8031592a…`; alpha verbose `REF-c98e198f…`; GenericIon(2,4) verbose `REF-7537d7e4…`; carbon SS 300 mm
`REF-61e64b1b…`; proton SS `REF-285075c5…`; carbon SS 1 mm `REF-35df82fcd7723f6a9a0c-db361e98`. The full ids of the
first nine are to be copied from `list_reference_artifacts` into this table by a dated note before the freeze
commit (not established by the planner). All are execution-only: no dose, IDD or phase-space output was read; the
placeholder scorer output was deleted unread.


**Dated note (2026-11-01, orchestrator) — full run ids of the ten smokes** (controlled engine tool, OpenTOPAS 4.3 / Geant4 11.4.2; nine at 41a1f8ac/a66867de, the slab smoke at ce4d3d6a):

| case | run id |
|---|---|
| `proton-water-150mev-em-smoke-p4` | `REF-188fa1d8600f598345fd-bb3772ee` |
| `he-water-150mevu-em-smoke-alpha` | `REF-ed13a2383682b780a086-7c1daae9` |
| `he-water-150mevu-em-smoke-genericion` | `REF-c764466cfe377d8eaf69-f82cd01a` |
| `carbon-water-290mevu-em-smoke` | `REF-5910656de47adfc817ef-e043426a` |
| `oxygen-water-350mevu-em-smoke` | `REF-8031592a201b7ddf5a48-7a2c2ae5` |
| `carbon-water-290mevu-ss-smoke` | `REF-61e64b1b80f37f99f9e5-02b74b6c` |
| `proton-water-150mev-ss-smoke` | `REF-285075c563439dac3279-432c3126` |
| `he-water-150mevu-em-smoke-alpha-verbose` | `REF-c98e198face1e783fde0-4f13dcfb` |
| `he-water-150mevu-em-smoke-genericion-verbose` | `REF-7537d7e4c5d8e4065d3e-cac2340a` |
| `carbon-water-290mevu-ss-slab1mm-smoke` | `REF-35df82fcd7723f6a9a0c-db361e98` |


**Dated note (2026-11-02, orchestrator) — files added by the freeze commit (sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/p1c_attribution.py` | `5ec7482a0604b78129e38ca59bdd0b29acccdbd1c423aa398ba33093661b500a` |
| `tests/ionmc/test_ion_em_p1c.py` | `d95597c79d74e17e1cf5ef20e2436f719716a359527239435afc9f8f6b07f13e` |
| `validation/plans/v3-008/p1c-aggregates.json` | `b603c3714c5d701c5560d6e348bab43526aaecb962a4e6dc016229fdd46306d4` |
| `validation/plans/v3-008/p2b-follow-ups.md` | `9c70eb1b9d49f0375af7976232817dbda925576eeb4b64e31aee61b33dfb5457` |
| `reference_requests/interventions/IR-20261102-010435-B59B2D/request.json` | `e81a7efa104614e87188c3ef083ef3fb061eb018d2731c46be2cb4a88e1ea7f3` |
| `reference_requests/interventions/IR-20261102-010435-B59B2D/README.md` | `22712340fcd38179571ac4f7f056b8364b99136e9586040b234cce0040ca451b` |


**Dated note (2026-11-02, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `tests/ionmc/test_ion_em_moliere_hanson.py` 789fd4b1… → 64a40829…; `src/ionmc/physics/stopping.py` 951f29ac… → 89db1de2…. All other entries unchanged (recomputed).

**Dated note (2026-11-02, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `validation/scripts/ion_em/p1c_attribution.py` 5b61e333… → 5ec7482a…; `validation/plans/v3-008/p1c-aggregates.json` dd2a2d0f… → b603c371…. All other entries unchanged (recomputed).

**Dated note (2026-11-02, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `tests/ionmc/test_ion_em_p1c.py` 031d45f9… → d95597c7…. All other entries unchanged (recomputed).

**Dated note (2026-11-07, orchestrator) — files added by the freeze commit (sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/gen_r2_bundles.py` | `32d7c55ea9d44faea595789de6dfac275516848182209d338b3733589a8c80a0` |
| `tests/ionmc/test_r2_bundles.py` | `aaa08b42abd3ed660530ba8897cb980fec3e9b3977dfddf7c7a699bc6bcb51ac` |
| `reference_requests/operator-directives/2026-11-07-topas-validation-target.md` | `ca079e9f0af990cdff9d8efde60d38ac55e2d0fc0fb9ef22355a81f4a2ab0b26` |
| `reference_requests/interventions/IR-20261102-010435-B59B2D/response.json` | `dc6a610cb67dee50fa0bb96327026eea0b81bb9506e7b898c2e9e44322f69bd0` |
| `validation/plans/v3-008/r2-topas-target.md` | `dc6b80d2a00de765f4458e6c7416ab31437a42ec47e2953111ab00a205316efe` |


**Dated note (2026-11-07, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `validation/plans/v3-008/p10-single-scattering-availability.md` bea86277… → a7f539a0…; `validation/plans/v3-008/p11-exact-omega-mapping.md` 6116648a… → f78e95f5…; `reference_requests/interventions/IR-20261102-010435-B59B2D/request.json` 068cafd7… → e81a7efa…. All other entries unchanged (recomputed).

**Dated note (2026-11-07, Codex REVIEW-7a7186a2 round R2-TOPAS-2) — revision-2 hashes regenerated; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` 90732599… → bfd1e333… (non-disclosing smokes); `tests/ionmc/test_r2_bundles.py` dd1a7af0… → d95f23f7…; `validation/plans/v3-008/r2-topas-target.md` f3995f57… → 9b253154… (appendix A, smoke paragraphs); `validation/reference_cases/README.md` rebound to e16cc9d3… (R2.2; revision-1 value 2c160fb5… superseded). The rows above carry the current values (one current entry per path); the bundle-set digest is `4a091ff6…`.

**Dated note (2026-11-07, Codex REVIEW-4c5b6803, R2-TOPAS-3: smoke scorers restored; blinding by ordering) — supersedes the R2-TOPAS-2 note above for the smoke bundles.** The 13 smoke bundles (26 files) are byte-identical again to their b8c9836f versions (representative, configuration-specific scorers; the R2-TOPAS-2 placeholder variant is withdrawn); the bundle-set digest returns to `9bd6f8023b4aa4b0bebbb21c5e237ea5a71c136dbc5243dbef5ef772c54a212f` (the value `4a091ff6…` of the R2-TOPAS-2 note is superseded, kept as history). The rows above carry the current values (one current entry per path): `gen_r2_bundles.py`, `test_r2_bundles.py` and `r2-topas-target.md` rebound, appendix-A hashes identical to b8c9836f.


**Dated note (2026-11-07, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `validation/scripts/ion_em/gen_r2_bundles.py` bfd1e333… → 0ebdf1b3…; `tests/ionmc/test_r2_bundles.py` d95f23f7… → 76d445bf…; `validation/plans/v3-008/r2-topas-target.md` 9b253154… → edac3266…. All other entries unchanged (recomputed).

**Dated note (2026-11-07, Codex REVIEW-e637178c, R2-TOPAS-4: ancestry-based run validity; accurate smoke wording; intervention status) — supersedes the bundle digest of the R2-TOPAS-3 note.** (1) The run-validity rule against the seal is restated in `r2-topas-target.md` §6 and the plan by git ancestry and controlled records (`code_sha` has the seal commit as ancestor; the seal commit's review finished before `reference_started`; evaluator `validation/scripts/ion_em/r2_evaluate.py` fails closed). (2) The 13 smoke bundles (26 files; the 192 evidence bundles are byte-identical to 274de574) carry accurate post-seal archival wording instead of "deleted unread"; the bundle-set digest is now `a8a1c27d0177d38e763344248753f0cb9792bbf2277158622b41bcbf9d8f423a` (`9bd6f802…` of the R2-TOPAS-3 note is superseded, kept as history); `test_r2_bundles.py` additionally rejects "deleted unread", "unread", "non-disclosing" and "SmokeEdep" in every bundle. (3) The intervention `README.md` has a status-update section appended (original lines verbatim): resolved 2026-11-07 20:37 UTC per `request.json`/`response.json`; the revision-1 entry (`86eba08d…`) is superseded history and the table row carries the current value. Rebound in the rows above: `gen_r2_bundles.py`, `test_r2_bundles.py`, `r2-topas-target.md`, the intervention `README.md`.


**Dated note (2026-11-07, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `reference_requests/interventions/IR-20261102-010435-B59B2D/README.md` 86eba08d… → 22712340…; `validation/scripts/ion_em/gen_r2_bundles.py` 0ebdf1b3… → fe122218…; `tests/ionmc/test_r2_bundles.py` 76d445bf… → 1cc41128…; `validation/plans/v3-008/r2-topas-target.md` edac3266… → 5336266b…. All other entries unchanged (recomputed).

**Dated note (2026-11-07, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `tests/ionmc/test_r2_bundles.py` 1cc41128… → 54af24fa…. All other entries unchanged (recomputed).

**Dated note (2026-11-08, implementer for the orchestrator) — post-freeze implementation commits.** (a) The implementation commits after the revision-2 freeze commit are `53424cfa` (light-ion species transport, CI rows, r2_evaluate chronology checks, revision-2 step definitions), `d235e2d3` (Codex REVIEW-5f5ac50d fixes) and the forthcoming commit of the 2026-11-08 working tree (Codex REVIEW-d0f41eff fixes and the mapping tabulation; dated note "2026-11-08, Amendment" at the end of this record). The freeze commit `0245158b` and the R2.1 rule are unchanged: the R2.1 whitelist governs the freeze commit's diff against `d0829601`; these commits are the post-freeze phase foreseen by the plan (A-rows → seal). (b) Rebound in the rows of this record (current values in place; old → new, short hashes; computed from the working tree by `sha256sum`): `pyproject.toml` 99ff83af… → 0e73bd64…; `src/ionmc/config.py` 2ead8d09… → f8b1c447…; `src/ionmc/physics/exact_omega.py` e8427d23… → 1b67ae5f…; `tests/ionmc/test_ion_em_mapping_src.py` e11b7b65… → d81e0915…; `validation/plans/v3-008/p11-exact-omega-mapping.md` f78e95f5… → b40530f4… (dated note 2026-11-08 appended; both its §2 row and its R2.2 row carry the new value); `validation/plans/v3-008/r2-topas-target.md` 5336266b… → 474480e7… (dated §8 note 2026-11-08). Newly bound implementation paths (15):

| File | sha256 |
|---|---|
| `src/ionmc/physics/ion_mcs.py` | `42194ccf6ede781b6254a138e8c64aacd8a72cf5385232c9e2758ea82edb9c81` |
| `src/ionmc/simulation.py` | `4d490ba490f521169d52040476662994b5e06eb24d3fb66f4367f290be7c6c28` |
| `src/ionmc/species.py` | `19df88c722ad93ca65e636b870a4b56b1a80df2e215bee39a7623f441daa2322` |
| `src/ionmc/transport/kernels.py` | `db4ea1f096e4919303b5e1d1413f25d78dd7ced149c331983016172b0b099266` |
| `src/ionmc/transport/reference.py` | `f3e3dec9ed5c7bffed827213783af19eee2bbe833f47120550c2c88a0a6b9bb3` |
| `src/ionmc/transport/run.py` | `bbfae2a4c1c2168397ab898204ede3bc2d7c92e1968a08b8156a03326217436b` |
| `src/ionmc/transport/tables.py` | `45542d46fb208ba1ca9ccdcb9e8baea41e47f86a7429b46f8566ad12b1b2619d` |
| `src/ionmc/transport/tally.py` | `be202f7310e3991bead49fd40e24cef1f5493992044c72aff6c5f8e8826b547c` |
| `src/ionmc/transport/warp_driver.py` | `cb32245e48d6bffd03092d09fb7421d8b11a8191158a311fd0b9a92e568aeb86` |
| `validation/scripts/ion_em/r2_evaluate.py` | `d1da57b57d8d55f9a545544eb773ff3c0a019e07c6663ae092779754248f5a05` |
| `validation/scripts/transport/steps_v8.py` | `e8b73006921b479797e762c1319adf7e1e4edb95594ffac4109827b4ee185a7d` |
| `tests/ionmc/test_ion_mcs_table.py` | `b3fb9992f866238b11243e032f394583775e292cec49601798ef44b54efcb749` |
| `tests/ionmc/test_v3_008_rows_ci.py` | `ffea9f20976ec2d721067ccf52d5cb9797fad51a5bda7dd224948829c1cade0e` |
| `tests/ionmc/test_r2_evaluate.py` | `aa3e28327f0ff5502ecdffb281cf356484890104333418d7ed7bc68633d13e53` |
| `tests/ionmc/test_v3_008_steps_v8.py` | `ab4922ea8842751e22c9d1882e9e1dec0ae88a25b55174c99586c4802868b731` |

(c) R2.5 counts, superseded by this note (the R2.5 text is kept as written): repository paths 475 + 15 = **490**; total full sha256 values 519 + 15 = **534**, plus the bundle-set digest (unchanged, `a8a1c27d…`). The count in R2.6 item (3) (475) is superseded accordingly by 490 for the acceptance tests of the seal.

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-7, Codex REVIEW-d3a17ebf: per-table certification of the tabulation bound; attest_verify bootstrap migration-only; manifest sentinel) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` 1b67ae5f… → fbeaaf18…; `src/ionmc/config.py` f8b1c447… → 0b9b34a3…; `validation/plans/v3-008/r2-topas-target.md` 474480e7… → 7aa4d9ee…; `src/ionmc/simulation.py` 783fa561… → cfb8f52d…; `tests/ionmc/test_ion_mcs_table.py` ce87db6f… → 95b02229…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts of the 2026-11-08 post-freeze note (490 paths / 534 values) are unchanged: no path added or removed.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-8, Codex REVIEW-9529fd35: analytic domain-wide certification of the tabulation bound with the relative-contract scope TABLE_MIN_LOSS = 1e-4) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` fbeaaf18… → 53bbc7f9…; `src/ionmc/config.py` 0b9b34a3… → 39bddffb…; `validation/plans/v3-008/r2-topas-target.md` 7aa4d9ee… → a90dc6df…; `src/ionmc/simulation.py` cfb8f52d… → 1d4d701b…; `tests/ionmc/test_ion_mcs_table.py` 95b02229… → 03e1925f…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts (490 paths / 534 values) are unchanged.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-9, Codex REVIEW-ee53867c: the table demoted to an accelerator of the exact root solve on every step; no interpolation-bound claim) — rebound entries (old → new): `tests/ionmc/test_ion_em_mapping_src.py` d81e0915… → 865d3143…; `src/ionmc/physics/exact_omega.py` 53bbc7f9… → 5f4aa165…; `src/ionmc/config.py` 39bddffb… → 34be1169…; `validation/plans/v3-008/r2-topas-target.md` a90dc6df… → 4722d837…; `src/ionmc/physics/ion_mcs.py` 709c61b1… → 36f4f3a7…; `src/ionmc/simulation.py` 1d4d701b… → b855273a…; `src/ionmc/transport/kernels.py` 73bcaf87… → db4ea1f0…; `src/ionmc/transport/reference.py` 103672f5… → f3e3dec9…; `src/ionmc/transport/tally.py` debbe3ed… → 7a746c99…; `tests/ionmc/test_ion_mcs_table.py` 03e1925f… → 62d495bc…; `tests/ionmc/test_v3_008_rows_ci.py` f40b7006… → efce4972…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts (490 paths / 534 values) are unchanged.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-9b, Codex REVIEW-4bf4f5bc minor: the `mcs_table_bracket_miss` descriptions now state that the counter counts every miss of the initial tight bracket, widened or full; documentation only) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` 5f4aa165… → e43eecc3…; `src/ionmc/simulation.py` b855273a… → 2958dbd3…; `src/ionmc/transport/tally.py` 7a746c99… → b941a46b…. All other entries unchanged (recomputed against the tree). Counts unchanged (490 paths / 534 values).**

**Dated note (2026-11-08, ionmc sample chunk plan appended to r2-topas-target.md section 8; records only) — rebound entries (old → new): `validation/plans/v3-008/r2-topas-target.md` 4722d837… → 346ad951…. Counts unchanged (490 paths / 534 values).**

**Dated note (2026-11-08, steps_v8 writes the sample npz with np.savez_compressed (lossless; the lateral grids are mostly zero, ~385 MB → tens of MB per chunk; chosen before any sealed sample, the single chunk run so far is repeated) — rebound entries (old → new): `validation/scripts/transport/steps_v8.py` 8d2617cd… → e8b73006…. Counts unchanged (490 paths / 534 values).**

**Dated note (2026-11-09, Codex REVIEW-16c63514 of the sample seal commit 84b8b2cd, decision D1 mapping-domain rule; orchestrator) — hashes regenerated at the working tree; entries that changed since the 2026-11-08 post-freeze notes:** `src/ionmc/physics/exact_omega.py` e43eecc3… → 1c8f3945…; `src/ionmc/config.py` 34be1169… → 513c32c5…; `validation/plans/v3-008/r2-topas-target.md` 346ad951… → 2e41ca0e…; `src/ionmc/physics/ion_mcs.py` 36f4f3a7… → 42194ccf…; `src/ionmc/simulation.py` 2958dbd3… → 6efadd7f…; `src/ionmc/transport/tally.py` b941a46b… → be202f73…; `tests/ionmc/test_ion_mcs_table.py` 62d495bc… → 7a86e44e…; `tests/ionmc/test_v3_008_rows_ci.py` efce4972… → ffea9f20…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, Codex REVIEW-f3361441 round, MCS-off gating; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `src/ionmc/physics/exact_omega.py` 1c8f3945… → 6471c388…; `src/ionmc/config.py` 513c32c5… → a24aa6b4…; `src/ionmc/simulation.py` 6efadd7f… → 4d490ba4…; `tests/ionmc/test_ion_mcs_table.py` 7a86e44e… → b3fb9992…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, seal-note commit; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/plans/v3-008/r2-topas-target.md` 2e41ca0e… → d4d177e3… (seal cross-reference note; an intermediate value 17673937… was superseded within the same working tree by the repair of a shell-mangled sentence). All other entries unchanged (recomputed).

**Dated note (2026-11-09, sizing-note and evaluator commit; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/plans/v3-008/r2-topas-target.md` d4d177e3… → eb7c2622…; `validation/scripts/ion_em/r2_evaluate.py` be1eadcc… → 6f92b41a…; `tests/ionmc/test_r2_evaluate.py` 5b996ebb… → dd9ab813…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, sizing-note and evaluator commit; orchestrator) — newly bound implementation and evidence paths (sha256 as committed; the evidence records index will be rebound as evidence runs are added):**

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/r2_io.py` | `827170cfa994ffedd709fd998cdca32f76004851bcd9e5752c3e504595e81e48` |
| `validation/scripts/ion_em/r2_statistics.py` | `fbb19370ae0bb20ed1cc5289864c3331ba63bdcb88315b33f06b6d710c79ed6b` |
| `validation/results/transport/v8r2-topas/smokes.json` | `9f5c1b24bfe278869b30cc89c03137786fc683ed85d6b6509c48cdffa964f621` |
| `validation/results/transport/v8r2-topas/records/review.json` | `87988397d4c14304d75f987117e4c22199d4452a0b48030e131a2bd52edb1199` |
| `validation/results/transport/v8r2-topas/records/index.json` | `264fba78aa314c04ac72574eb02e79a9e58798238de77ca5a85adaa9da209b98` |


**Dated note (2026-11-09, R2-SIZE-1 bundle regeneration; orchestrator) — supersedes the bundle-set digest `a8a1c27d…` of the rows above.** The bundle set is now 209 bundles / 418 files with bundle-set digest `e51d3b3d1bf00e51daadf889f6541bf8bd6ba8aea6453eee132110e07ff0c4d6` (`gen_r2_bundles.py --check`; per-file sha256 of the 168 changed or new files in `r2-topas-target.md` appendix A, 2026-11-09 listing; the 250 unchanged files keep their appendix-A hashes). Changed: 48 LAT-FULL evidence bundles (`CutForElectron = 1 km`, texts), 32 carbon/oxygen IDD-FULL bundles (`NumberOfHistoriesInRun` 3333 / 2529). **(2026-11-09, later; Codex REVIEW-6b4731a7)** every evidence bundle and the four new smokes additionally request `b:Ph/ListProcesses = "True"` (runtime evidence of the executed module set, g4decay included); bundle-set digest now `6d1ae3a06211831cfb7e8924c65d6d8fe31fae7781443ff84c077204045ea790` (appendix A 2026-11-09 "ListProcesses" listing; the thirteen smokes already run are byte-identical). New execution-only smokes (sha256 of `case.json` | `input.txt`):

| Bundle | Seed | case.json | input.txt |
|---|---|---|---|
| `he-water-250mevu-r2-lat-full-ecut1km-timing-smoke` | 20600424 | `6567ac23198f5b65a0405a6209200409daec0504c4907d40c203b0408369ddc6` | `88ebffbe8c7bbef064fc324adc4a71347f411582ce6111e55c836da578370bb5` |
| `carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke` | 20600425 | `94f255cd117727a4370bec056c4cef18944d670abcdca31ba8adc6794e0cdc16` | `f66cd544dcca03d43a90bbb1984a8e3e7267d611e613888e853caff50b2c2ddd` |
| `oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke` | 20600426 | `f81f382d3ee31d505d08ea30b63506a3055df72e231a5be88c7d2bcf7928ddac` | `21bbac3a4881ad9516eafd8cff5748dce5dbfb7f1ff93104e42d06e8fd8e6175` |
| `carbon-water-290mevu-r2-lat-em-steps-verbose-smoke` | 20600427 | `812c1335cd298ae102af825b80d781409b83d0dc612fd48c74b83637ec74aead` | `37f21fee0a6e64b3a90e9c67e0865a72bab4f6764234cce41db210872fa76ddb` |

**Dated note (2026-11-09, R2-SIZE-1 bundle-regeneration commit; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` fe122218… → f970715a…; `tests/ionmc/test_r2_bundles.py` 54af24fa… → f767104f…; `validation/plans/v3-008/r2-topas-target.md` eb7c2622… → 1999e5cf…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, evaluator review-fix commit (Codex REVIEW-01b6acd6: signed start record, executed module-list check); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/r2_evaluate.py` 6f92b41a… → c74faab2…; `tests/ionmc/test_r2_evaluate.py` dd9ab813… → 478dd540…; `validation/scripts/ion_em/r2_io.py` fa129516… → a4aae6c9…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, evaluator round-3 fixes (Codex REVIEW-6b4731a7) and ListProcesses bundles; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` f970715a… → 8a92d40b…; `tests/ionmc/test_r2_bundles.py` f767104f… → fd87a7dc…; `validation/plans/v3-008/r2-topas-target.md` 1999e5cf… → deade6b6…; `validation/scripts/ion_em/r2_evaluate.py` c74faab2… → 441a6a82…; `tests/ionmc/test_r2_evaluate.py` 478dd540… → 2c9d52b7…; `validation/scripts/ion_em/r2_io.py` a4aae6c9… → 4c2cabeb…; `he-water-250mevu-r2-lat-full-ecut1km-timing-smoke` 6567ac23/00428331… → 6567ac23/88ebffbe…; `carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke` 94f255cd/92c02974… → 94f255cd/f66cd544…; `oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke` f81f382d/95ca1a53… → f81f382d/21bbac3a…; `carbon-water-290mevu-r2-lat-em-steps-verbose-smoke` 812c1335/4e9e5179… → 812c1335/37f21fee…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, evaluator round-4 fixes (Codex REVIEW-a50374b1: evaluator pinned to its own checkout, committed-HEAD inputs, frozen module sets); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/r2_evaluate.py` 441a6a82… → 63ca633c…; `tests/ionmc/test_r2_evaluate.py` 2c9d52b7… → 773f2b4f…; `validation/scripts/ion_em/r2_io.py` 4c2cabeb… → 417d0de7…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, frozen bundle manifest; orchestrator) — newly bound path (machine-readable manifest of the 209 bundles, bundle-set digest 6d1ae3a0…, verified by `gen_r2_bundles.py --check`; every TOPAS evidence run is bound to it by the evaluator):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/frozen-bundles.json` | `f5e39fe6c7c1c846070e34970024717d86c36dd04b8b287b993a07947856400d` |


**Dated note (2026-11-09, evaluator round-5 fixes (Codex REVIEW-096af6a3: chronology inputs pinned to the run's code_sha, frozen bundle manifest); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` 8a92d40b… → 9924429d…; `tests/ionmc/test_r2_bundles.py` fd87a7dc… → b0e91e2b…; `validation/scripts/ion_em/r2_evaluate.py` 63ca633c… → ad7c45f1…; `tests/ionmc/test_r2_evaluate.py` 773f2b4f… → f733732e…; `validation/scripts/ion_em/r2_io.py` 417d0de7… → ec08eea5…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, re-sizing note and smoke record update; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/plans/v3-008/r2-topas-target.md` deade6b6… → 22fd4067…; `validation/results/transport/v8r2-topas/smokes.json` bae61e85… → 9f5c1b24…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, smoke run-validity amendment (Codex REVIEW-8429c50e); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/plans/v3-008/r2-topas-target.md` 22fd4067… → beb35820…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, R2-SMOKE-LP; orchestrator) — supersedes the bundle-set digest `6d1ae3a0…` above.** Ten ListProcesses smoke variants added (execution only): bundle set 219 bundles / 438 files, bundle-set digest `76d785667211bb39e96d517f46d636ff4b60be212b0d7fda9d968993ef47ed66` (`gen_r2_bundles.py --check`; appendix A 2026-11-09 "R2-SMOKE-LP" listing); all 209 earlier files byte-identical. New smoke bundles (sha256 of `case.json` | `input.txt`):

| Bundle | Seed | case.json | input.txt |
|---|---|---|---|
| `he-water-250mevu-r2-idd-em-lp-timing-smoke` | 20600428 | `33942503c5eb0370e08333c0ef9c68925b886abbf387b2899f087c00af807266` | `09acb762d27c0de2b1f3264d3359c89fe6176b52bf20baf8a7bd1dee9d7142fd` |
| `he-water-250mevu-r2-lat-em-lp-timing-smoke` | 20600429 | `f9641a407b6ccd44e5327a10f4374b57f479f78cb8c6a5c7e98ad5a0e405dbf1` | `4115bab225ef5d42ab1b35d5f4a324762de11e36c15ca0f36fa2c12be8d38902` |
| `he-water-250mevu-r2-idd-full-lp-timing-smoke` | 20600430 | `34a80eb639b62719911e368d892dbbaeebb496b07dda97c77d86fcf6369823ff` | `0deebb528d989e186723872e31e2ed07879e213347969a2fa6349c0cc9b19947` |
| `carbon-water-430mevu-r2-idd-em-lp-timing-smoke` | 20600431 | `660fa09e1e2cf1e656de05be99562c2ba541660a15a6547b6d0cf1f07dcc1bc7` | `21c5ab7d7df1a1582776e7bb3c1ba409ec744093e79ec4fbdc12a70c57fc8888` |
| `carbon-water-430mevu-r2-lat-em-lp-timing-smoke` | 20600432 | `14b1441542169823520c6b56e88743d35c01c26db065aa8bdb95da4f767fee54` | `f7e5a6ac015d368029a39ce14f33d9b4b6adfa05b54841522aa2a1ba3a981579` |
| `carbon-water-430mevu-r2-idd-full-lp-timing-smoke` | 20600433 | `c1ca492f2f929193d968952de0673a5a594281da7b0e171833ccf88a5324d311` | `dc358b90d13f34b368238f1f47596b8078efb37d296d3141ec79e04dd8b47cc7` |
| `oxygen-water-430mevu-r2-idd-em-lp-timing-smoke` | 20600434 | `074de84a1b73b82c9a4169b1af2892db79f0af6d6089eed24b3901b9214d72b5` | `8bd26ee79201156ce003efae5bb6a0a151f6a125906e107d3b55e719be435df4` |
| `oxygen-water-430mevu-r2-lat-em-lp-timing-smoke` | 20600435 | `b8183c78e3429b947ed4e1040ef4e5a4ff61fb2a1c4def3820ad350d330786fe` | `7d09eded2d6761b8820451e4ef9d43ec2448b0420e1127d8346e973d3c99e271` |
| `oxygen-water-430mevu-r2-idd-full-lp-timing-smoke` | 20600436 | `367db7e36c6304880ffd058e1f998710863db0f2ef09ef7ec301f5b46aedbc33` | `7814b35838e2635a7fc4de267793389ffdc21ebb73ef75b9349a7039d2f8fc79` |
| `carbon-water-290mevu-r2-idd-em-lp-steps-verbose-smoke` | 20600437 | `05b0c010132dd5e50934f46eee32cb5952345f4929394a4282dd34affc503b49` | `75be3a5823c98902d04fb47b65180fc88a827a94b667de703953ce83c62c10ec` |

**Dated note (2026-11-09, signed smoke re-runs; orchestrator) — newly bound paths: the check-run summary, the smoke-rerun summariser and its test, and the 14 signed smoke-run records (record.json + artifacts.json; sha256 as committed; the records also carry the reference-runner signature):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/smoke-reruns.json` | `4a3fb8af1163dea56e0c8d6868b127f4277ccb02d7a80b133bddf0eb3511ff58` |
| `validation/scripts/ion_em/r2_smoke_reruns.py` | `2411d900b7b2aac36629edd7829f3cfc44e53a3807b0677a77d8893e2c8d2cb6` |
| `tests/ionmc/test_r2_smoke_reruns.py` | `85541edab8e4a5a54a578178b911884c9e394d046411deb84b41925dc7f68144` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-127b2dd92a8dd609c560-ff8db64b/artifacts.json` | `b2fb14264923e3d9d4b49424b8df72fe6d1d060ac73ebd57b423a3f44fa18dd4` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-127b2dd92a8dd609c560-ff8db64b/record.json` | `a79ab4244232b39964db8fa4beb43f78cd7e0d5102cd1465f6979a3bf9f6e8f6` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-4e5b0681e5576d7ab451-7e640368/artifacts.json` | `c4bcf072e554a1e5882bc4ac06430b3fdf35a543ce15325a196cc21383d7cf8f` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-4e5b0681e5576d7ab451-7e640368/record.json` | `e47eb958c4457dfef7cddc84af00435fa71d6d15a69f3b7b8b4e8058cae77f98` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-7020741cc2707c255204-a567db01/artifacts.json` | `671a4c9ef68c77e24be5b977615ccc617f17e6a41a7534144ce05644ee21abba` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-7020741cc2707c255204-a567db01/record.json` | `c39c996349872d9d945199cbfd81673d7b1734def8b55c2c518ae56a5acf6035` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-84df1a6d3f6427a74474-1306ceb2/artifacts.json` | `5c6fbf932c0ba166304082a12aca4ec54881a36695559ba7b5059a9efd139000` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-84df1a6d3f6427a74474-1306ceb2/record.json` | `6ca770fd0a54141f666e06b755dac95353e30d26a0c6bc1459f00cfd32ad5bf4` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-878a53365185a3360106-6ceac04f/artifacts.json` | `99df27a033fa8b713872ff08045b93874a2b6beeceff9f5a86fc42601467009d` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-878a53365185a3360106-6ceac04f/record.json` | `c3d5e77532cce6ad1e8cbbe04d4a1e0d793edd20898ba80ab273d254702f4269` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-87d804e7bd55f38370b0-4e7af629/artifacts.json` | `5c5d7d59734a3a79c87451ecae0f18ef60899020def5c6f3c9eb91729d0de118` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-87d804e7bd55f38370b0-4e7af629/record.json` | `d5a683a2c48070fc34f974d4afdbc371c28c7de1030657c7d5f616693c9e1ce1` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-8bf4ca70fe3b29de6260-11223111/artifacts.json` | `2473b2f20c185d60ffadec2c8354967c65f4694d60c4046d73e2dddd7bfdb6a5` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-8bf4ca70fe3b29de6260-11223111/record.json` | `251466978adc72af93182cacd70d7a890ce8f7374f5c41615fc5154df2a774e4` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-95f5caff95494e521eba-c5edd170/artifacts.json` | `afe678f13ee66e61823251c82caf78c55d7f927ab60627baab154506905cd10b` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-95f5caff95494e521eba-c5edd170/record.json` | `a83bfcc6f758cb8340dd1b654261e6e293cd07ea7f2417cf2cabe93736bc405f` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-97622ace316bf969f405-947c0413/artifacts.json` | `f20307a85fbed43cacd48ba3a6d1f8c3fdd709f9af9a973b588665cf3c1e9920` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-97622ace316bf969f405-947c0413/record.json` | `9a5212e874cd3a2050e57790ae6068d8e24777be95ac9d98c4add0d967921279` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ac5f93bc5cc1507f3e27-8af01fff/artifacts.json` | `3ecff8321c6b14b4ced540fe72c253ff85ad09408691fd2cbddbb64dd67a7167` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ac5f93bc5cc1507f3e27-8af01fff/record.json` | `7988ed8de6f9411f73b350ee52f0aa233efdfa435448bc120c3a479f055957d1` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ba0afbd01523d5d68398-153f9b26/artifacts.json` | `71abf85d7e813a2b592599e3df664a87f64020997876ec53e871be5697acbf35` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ba0afbd01523d5d68398-153f9b26/record.json` | `61dac024afdcdd0cd44206eaafe13309b8a708508eb7d36143fc18a58362398e` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c32fd0fa07582b3fdb5b-717b2259/artifacts.json` | `78537b0cf0ab6b8644d2ce56ae3cd6faa6bc7f25411d0037e714666a7897b8aa` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c32fd0fa07582b3fdb5b-717b2259/record.json` | `df815a116dd789209e21ae6ebf05d699163d9400713dd76fdbcd0f18e54b55d7` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-d784c6d81d2a2ca6df30-313c4c31/artifacts.json` | `e29aa535344a4ddf3498ed33c25939d8490c4c20dd562c9af099cc4d365d42ce` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-d784c6d81d2a2ca6df30-313c4c31/record.json` | `ef6b12a766e4b002ee789f7afb37576a0e212377861ba447abe91afdf316921e` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-f1cb893d737b990068a3-bdb32d20/artifacts.json` | `efafa54d89c3f8df7564407e4d0ae730ac04de2a293e43a6dc572f9df7ea0113` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-f1cb893d737b990068a3-bdb32d20/record.json` | `641d76972b7417458e9bd45b4a2843e7ab6e4b31601406cdfb4566849270701c` |


**Dated note (2026-11-09, signed smoke re-runs / STOP / R2-SMOKE-LP commit; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` 9924429d… → 74ab35e3…; `tests/ionmc/test_r2_bundles.py` b0e91e2b… → e6431320…; `validation/plans/v3-008/r2-topas-target.md` beb35820… → 3b2a5a40…; `validation/scripts/ion_em/r2_evaluate.py` ad7c45f1… → 1aa72236…; `tests/ionmc/test_r2_evaluate.py` f733732e… → 34bdee01…; `validation/scripts/ion_em/r2_io.py` ec08eea5… → 621a4979…; `validation/results/transport/v8r2-topas/frozen-bundles.json` b3e27964… → fa490129…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, R2-SMOKE-LP correction commit (Codex REVIEW-8286be8c: lp name grammar, manifest-current confirmation set); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `tests/ionmc/test_r2_evaluate.py` 34bdee01… → 0de50eb4…; `validation/scripts/ion_em/r2_io.py` 621a4979… → e4ceb4b2…; `validation/scripts/ion_em/r2_smoke_reruns.py` ffee2f4f… → ffcd9d13…; `tests/ionmc/test_r2_smoke_reruns.py` 0b40e2ed… → d0eca2fa…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, confirmation-set records; orchestrator) — newly bound paths: the 14 signed smoke-run records of the R2-SMOKE-LP confirmation set (10 lp smokes, 4 manifest-current re-runs; record.json + artifacts.json, sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/smoke-runs/REF-022775f8c6a8530c48aa-a50056d1/record.json` | `c03a50cc1bfccee55ec3c725c6f64e379f13e2d68c52d859a811f05de3b96a68` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-022775f8c6a8530c48aa-a50056d1/artifacts.json` | `b52ec38d9a1092cb070ca6b4087d8e6971536e3cc3a2ca3e9b480da296e2aada` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-100a453268dacee7ea3f-80664e3e/record.json` | `1b4fa29d220e526a83e3716a8184601745276f68e0010d32130e4b13f800439d` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-100a453268dacee7ea3f-80664e3e/artifacts.json` | `8912090ab0146cd0e53616e7669c45b6940b536978def3075f9c4d952e3a58d1` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-1bee6c0ab49920cf6ff3-685816fc/record.json` | `0917666c242b4e412dd30cbd983eab415a7f7247a9e6c12cae4cadffafa8665a` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-1bee6c0ab49920cf6ff3-685816fc/artifacts.json` | `c585bf6e4f7aa102709895e9478efc7b01a6db290f36fd10ac72c44465707edb` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-2c6d31fa3f3ae83bce78-ea24d511/record.json` | `4a2a4d074ed4ee59e03db1f0bff6791562c3589c84be7493ddb4591f33769e72` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-2c6d31fa3f3ae83bce78-ea24d511/artifacts.json` | `cbd4ed3f28d2f12f7274129942f9b5cba4527a8808590cc5433e07df683b0c3e` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-494fafda858a0a1f2c2d-5e74e4ef/record.json` | `be6c8567bdd463a009ce77839643c81e7fc00696b46f33e0d1d12d002f28e670` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-494fafda858a0a1f2c2d-5e74e4ef/artifacts.json` | `f01b80d1daee1c51335bc817667ad0d349d8f948d64468e7d4b9c654306a7c21` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-706eae23ccef7407d213-b7cff919/record.json` | `2762ef3a6a39c296a5d809613e15af60a2251ee9718942b2625b4e96999226de` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-706eae23ccef7407d213-b7cff919/artifacts.json` | `6d5586cdffe6464d53880e70a14813002502038be72b810fdb928d1c6b069810` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-77d4d10efe5c035b544e-caf6cb7f/record.json` | `8802b8444dd184cd55452c013e891f8d15a7759e6cf15c03aa91be842bd29571` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-77d4d10efe5c035b544e-caf6cb7f/artifacts.json` | `fb399aaaf0c8d82b22865afcc26e236a89c75cd5fc12162a2300cead0a10b2b5` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-88218abcafa16a04e081-8930543c/record.json` | `d79f7cd545e93d847a0042a64d2d38fe43b106303bf073b1052b8a5fc167536c` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-88218abcafa16a04e081-8930543c/artifacts.json` | `45bc49c87c950cec1b437e553456908e72376753d04f35232629e527238dec7d` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-8823040474464055579c-29f6732b/record.json` | `867caf33b263e031f51314b35c2192a71a68ce9e987d61e6e3cdc8ec9ee99123` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-8823040474464055579c-29f6732b/artifacts.json` | `8652e67d320dcdfbf2fef8f75c100602d4697a82c2624cf5a017934d362839a1` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-964a865042b11a34b6bc-688112f3/record.json` | `76a85623fc3d32dee58bdb472f3bfa73fe2ae13d536f49673bc82f2d647b6eea` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-964a865042b11a34b6bc-688112f3/artifacts.json` | `e5477a82ef1b703e7fb5c4bf3e2c8d2ab5b9eb99c3497bb8f84c90d339aacd32` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ae7dbb3bd0e51a23a709-e348d838/record.json` | `abf7d9acac52e0f031c252a1932606a239510e17fe62851a995d945e5daa538b` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-ae7dbb3bd0e51a23a709-e348d838/artifacts.json` | `47bdb98f431c53fde6eb009d78e67767ce1ab65b7bf004f20a48d2f9879348b9` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c26d98d1bb1665012906-ef220531/record.json` | `06a64f19c3ae336e226f4e3d24e7fafd0aa139b7e3ab92842de04333ff152490` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c26d98d1bb1665012906-ef220531/artifacts.json` | `2874573b20ade948a293a90dc804c80813140c643825a1779f59c78a3685f5f2` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c7ab0adc969cd4034512-f5df42e7/record.json` | `c1bd2a4888cfb0c78632477474a9cf20100ccebc6a056e39862961da1e723c65` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-c7ab0adc969cd4034512-f5df42e7/artifacts.json` | `200e2de570613f0183ae992d57df782e911c892a8c4d2ae767f050d81fcb14fc` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-fda97285d2e8ab1adcb5-9699e0a1/record.json` | `0958e176ece17564893bbf6cce7e4498bdd0f20355fa74a0d9c3fac7012d461e` |
| `validation/results/transport/v8r2-topas/smoke-runs/REF-fda97285d2e8ab1adcb5-9699e0a1/artifacts.json` | `05c8cb218303c909fa069d38908b35c78275eafd7219174abef19c6c0300f91b` |


**Dated note (2026-11-09, confirmation + IDD-FULL re-sizing commit (C 3200 / O 2500 per part; bundle-set digest b5fa180fa5e8a7fee59ef22716048b1ad3e06e683897bd7962ea0c9f4e3afb75 supersedes 76d78566…; smoke runs bound per bundle entry, evidence runs to the complete manifest); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/gen_r2_bundles.py` 74ab35e3… → 32d7c55e…; `tests/ionmc/test_r2_bundles.py` e6431320… → aaa08b42…; `validation/plans/v3-008/r2-topas-target.md` 3b2a5a40… → dc6b80d2…; `validation/scripts/ion_em/r2_evaluate.py` 1aa72236… → d1da57b5…; `tests/ionmc/test_r2_evaluate.py` 0de50eb4… → aa3e2832…; `validation/scripts/ion_em/r2_io.py` e4ceb4b2… → 827170cf…; `validation/results/transport/v8r2-topas/frozen-bundles.json` fa490129… → f5e39fe6…; `validation/scripts/ion_em/r2_smoke_reruns.py` ffcd9d13… → 3fa0cb48…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, confirmation summary; orchestrator) — newly bound path: the check-run confirmation summary of the 14-run R2-SMOKE-LP set (all VALID at evaluator head 1d7b2394; sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/smoke-confirmation.json` | `9eb05abb2b850d10c473d4aaea6aabeb8a0d5710e335c05d6c64bebeece4ed1d` |


**Dated note (2026-11-09, exact confirmation-set gate (Codex REVIEW-27887e4b); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/r2_smoke_reruns.py` 3fa0cb48… → 33964577…; `tests/ionmc/test_r2_smoke_reruns.py` d0eca2fa… → 85541eda…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, confirmation summary regenerated at f1ef5a51 with the exact-set gate; orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/results/transport/v8r2-topas/smoke-confirmation.json` 15df91b7… → 9eb05abb…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, summariser docstring (Codex REVIEW-027cb008 minor); orchestrator) — hashes regenerated at the working tree; entries that changed since the note above:** `validation/scripts/ion_em/r2_smoke_reruns.py` 33964577… → b0a8e70e…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, orchestrator) — hashes regenerated at the freeze commit's parent; entries that changed since f694b747:** `validation/scripts/ion_em/r2_smoke_reruns.py` b0a8e70e… → 2411d900…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, orchestrator) — files added by the EM evidence-records commit of 2026-11-09 (the 72 EM evidence parts: IDD-EM parts 0–3 and LAT-EM parts 4–7 of cases 1–9, run at code_sha ead583f1 after REVIEW-5baeb516; signed run records and verbatim artifact listings; sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/records/runs/REF-04dae32f7e6e4b7c53d9-cd22c30a/record.json` | `a8aac6e6d8b04dc963c53a2aea05243da73e6c19f0e7cc00cf8bc266e0faaa05` |
| `validation/results/transport/v8r2-topas/records/runs/REF-04dae32f7e6e4b7c53d9-cd22c30a/artifacts.json` | `1af7490da5bad94536fb8ff5feaf295a6a8a7492c9ca8393af985793668a722a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-07e2fc253fd227f8865f-ef8acb6a/record.json` | `f6f7b13ee1355812703c08c1aac5c1148edd22bb242de8ce22762567d9884b55` |
| `validation/results/transport/v8r2-topas/records/runs/REF-07e2fc253fd227f8865f-ef8acb6a/artifacts.json` | `701294ec694456ac618731c634f2ef52edd6754c9cc6a7bd6afa858329c01b94` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0881c06b2c79cc46abe9-d9b3783d/record.json` | `990474c3fdca7ca0832190407af716262aae8d3f0335fbe313f8c979c05a0658` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0881c06b2c79cc46abe9-d9b3783d/artifacts.json` | `b90a3e4401c923bd1d4016d49d23378dfb3e8633ed8676650b05e0d7b1c451c8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0a1ffa6a5e64335d797a-6351fe0e/record.json` | `b9fa0d9deb187b3dbddd8d628f71e40020a92d4edf3467bb1ecaf74d18b9044f` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0a1ffa6a5e64335d797a-6351fe0e/artifacts.json` | `76d060bd03785c825c03b62581ff3e481185299f8e5ef12032885c339cf580a6` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0d6cbfab1d9b4bef9647-df1501a1/record.json` | `f98e843088cfab9e9c63a6b3e23993fb4a207642849194391e882b25a064d3eb` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0d6cbfab1d9b4bef9647-df1501a1/artifacts.json` | `a2a08e91fbb00a9a940bba271c4cec71c25f8b992d8e1fd982905721cf4d900e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0f6fa87af643aa024933-64402f0a/record.json` | `3b627b6897c27291d5e296d60aec8b29db5c1bae8abb5cf25b83f498417efcc8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-0f6fa87af643aa024933-64402f0a/artifacts.json` | `cde645808cd3047d14ea62aed8b6b52602e92deb83fa993a341e90cc949c3227` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1128ee8ff3ebb9a57dae-5eb6035a/record.json` | `0fe0afdcd8300b16c627450ca4a32f85c6ad56070f0a41369f806d29139a1641` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1128ee8ff3ebb9a57dae-5eb6035a/artifacts.json` | `abb598f1511d20a636d7c6bd20d710e7c2f3df4b6fa0af05cd86a277a327c368` |
| `validation/results/transport/v8r2-topas/records/runs/REF-11da44e21bce0ef5627f-59121724/record.json` | `9ff034624d92f27ba8bfa2aeedd1902e717bafa5938f56c3445f24dce72164b8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-11da44e21bce0ef5627f-59121724/artifacts.json` | `fc792977617f28741a098c476f734a071df049ee2da4e63e0e288e7d71faec37` |
| `validation/results/transport/v8r2-topas/records/runs/REF-12b6bbf50b273f244387-ad7f9f2b/record.json` | `a0de80dd244284b3bcc65397271d44d8beea547289b514ba5249cba488ca19c2` |
| `validation/results/transport/v8r2-topas/records/runs/REF-12b6bbf50b273f244387-ad7f9f2b/artifacts.json` | `d7ab27e142ab27bcafcf6238b7bb22aeec67a2a8c45d5300d9b60014f3a19faa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-13caa50a25cdc3a9e6fc-7f169e4e/record.json` | `cbdcc8005e641a65ab3924c6976af29e5a9a8dd0a0e1e44fa1cae85d830536e4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-13caa50a25cdc3a9e6fc-7f169e4e/artifacts.json` | `9789f669c65f8a8330c540f1b1ae7526532cba426eafaa15ebb822cad8ce07f8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-13e2a8a44d1f00e3d7df-88894365/record.json` | `a8e2c3f780d179bbb0589b93c0ff8797d5534b47482a368a8428dc74713a01ec` |
| `validation/results/transport/v8r2-topas/records/runs/REF-13e2a8a44d1f00e3d7df-88894365/artifacts.json` | `4571463b9055518fb57cc66ecc0fb988a83b14365aa2f6dab8b1364c7fb453d5` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1c6adae4c80435460bbd-1a2f16fe/record.json` | `e85e415c5ac64ffcfbba5d19dc67309d9a7b732c95256afad599b1cd9d010d1b` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1c6adae4c80435460bbd-1a2f16fe/artifacts.json` | `f2c88c73b3a2b4d12a34782b9bb87529c901a889718d5f5265d703b6bcd28148` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1dbb2ed332a1135f22b9-6694ebcc/record.json` | `fb272a6ded4a16377b3b3497c845d60adeaba882b220cc9c85a63b19ee665cec` |
| `validation/results/transport/v8r2-topas/records/runs/REF-1dbb2ed332a1135f22b9-6694ebcc/artifacts.json` | `428d93e4aeae4a7f8c167f1df9a6c73f61848ed7caf305aa1a606d02adcee058` |
| `validation/results/transport/v8r2-topas/records/runs/REF-22b2769f0653860495b9-01882251/record.json` | `77d748b5aa54e67489c7ba6493c29b816567eb27aecbb84ee70c20c7c4ab6620` |
| `validation/results/transport/v8r2-topas/records/runs/REF-22b2769f0653860495b9-01882251/artifacts.json` | `7bcac72f2753340793819920cb4666a730d2fea3ae0a7eba54da77450d9bb06a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-260acb9d47e48a902a61-0c7cfb06/record.json` | `ad847cb39accf87c55c6c9c1870a618da1bf2e333b428a03e59815d216df429a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-260acb9d47e48a902a61-0c7cfb06/artifacts.json` | `5a4498e65fcd36c28ec41e2c0f628f869ab168269b7e3c3d95f61a37e5647988` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2a3f5a44cdfac62f921c-8182cea9/record.json` | `3550b119f8be2096d3abe054743db4d7446f378208a720bcb79610068798198e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2a3f5a44cdfac62f921c-8182cea9/artifacts.json` | `0f62f58d2c8872e2d1dccc7dd49eab0c98f6503d95f25396e4ed6903e82998ee` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2ae51e027f7910f26a08-b8f8c673/record.json` | `ac81b01aeda788a475d79999efb6dcbd163c21e8cc81c15589dad859d13302fb` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2ae51e027f7910f26a08-b8f8c673/artifacts.json` | `415b666209760e90a5693c8012401a94055ac0db7dd971f8e459fcbc29ed0c5d` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2e66a62b0de5aec1aae2-43c82a33/record.json` | `2fac380519602643093cdc6af6eca5730567c098db5f84ac25eed712553f3dba` |
| `validation/results/transport/v8r2-topas/records/runs/REF-2e66a62b0de5aec1aae2-43c82a33/artifacts.json` | `17b5ad97677eb1a2cdb46a9935a033e607b2d33efca863487c149adba494dbb8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-3263419233fca6fb8e42-120e4a53/record.json` | `6451cb69039ce03c39c5c39f4ce535db626bccc4e457706d2290079307dfd89f` |
| `validation/results/transport/v8r2-topas/records/runs/REF-3263419233fca6fb8e42-120e4a53/artifacts.json` | `7be67a6196e409c29205d802bf9d9f29a7858331c0ac1d4f6b195acfc6ee349d` |
| `validation/results/transport/v8r2-topas/records/runs/REF-33697a5201767cc8e211-956a6d1a/record.json` | `02cb14768b5052cc9c9a0a75f8fa258c1bfcef3d3e2a74237a194da73ae7c9f2` |
| `validation/results/transport/v8r2-topas/records/runs/REF-33697a5201767cc8e211-956a6d1a/artifacts.json` | `e75c328820f810de49ac8327e766c8823c15f39e59b32c3e6ee6f652adef786b` |
| `validation/results/transport/v8r2-topas/records/runs/REF-33c53d05d0754f3f7d78-a60f9ea8/record.json` | `cf2cce902eb866bedb130dd3d6853868d1024ba47b9936817cd9c910e9665ad0` |
| `validation/results/transport/v8r2-topas/records/runs/REF-33c53d05d0754f3f7d78-a60f9ea8/artifacts.json` | `a0a5633c68181da94adc01145e9137d8b507dc00044f54bbe8e58eee6d94a5ac` |
| `validation/results/transport/v8r2-topas/records/runs/REF-377d191363dfd225d31e-6c60fce6/record.json` | `233d41cee7abcf50549b620059bb8578b3926842bb08434cdb13ffa484bb0f8d` |
| `validation/results/transport/v8r2-topas/records/runs/REF-377d191363dfd225d31e-6c60fce6/artifacts.json` | `39d6b68e5ee7c8a80bdf252625da1e6c2503b46da07da898e77ca28fcd0b3161` |
| `validation/results/transport/v8r2-topas/records/runs/REF-37b0c6d0fb895d4f64c0-04266851/record.json` | `fcccf9bd65c0d3445723ee405f8ee0363324d84790bcdae1b328197d54875287` |
| `validation/results/transport/v8r2-topas/records/runs/REF-37b0c6d0fb895d4f64c0-04266851/artifacts.json` | `9557d6feead43e88044a5ad74f6b311a92679db850187226fd9c4f3f42cb8620` |
| `validation/results/transport/v8r2-topas/records/runs/REF-39731c5d7e8122f2c417-5c2db0ad/record.json` | `991f9b301515fb1e92ea80f6a4b9a84cfdd39053544d1729d5dc38a87ab7affa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-39731c5d7e8122f2c417-5c2db0ad/artifacts.json` | `3188f626464715fd1c3c0d6e231ee3dac956bcadf0a738b564164c3732da1346` |
| `validation/results/transport/v8r2-topas/records/runs/REF-468d75c0234179822ef6-03056a77/record.json` | `afb20fa6853b6dac8e69e2ed685ff570c4281a477eb8904bec7a9bbf5f7fee79` |
| `validation/results/transport/v8r2-topas/records/runs/REF-468d75c0234179822ef6-03056a77/artifacts.json` | `f40967889543ce4b3f9c51083c531f52df6f135839712f8b80c588efae88544d` |
| `validation/results/transport/v8r2-topas/records/runs/REF-498b9fa7f3d232ad0360-68b55a0d/record.json` | `fcd2bc8af064849de97ad686bbb16acace8415496d16634a73700b265bebf366` |
| `validation/results/transport/v8r2-topas/records/runs/REF-498b9fa7f3d232ad0360-68b55a0d/artifacts.json` | `bb7b036084adfc7a74265f9ad3232a64a1a8bfca540c30240e6ec9f3bfb891ae` |
| `validation/results/transport/v8r2-topas/records/runs/REF-56b37ef1fd8dbf0789c3-ec02c119/record.json` | `84b40d698d7368b9ebfb59dcb42885c275afc9efc85f38c852194c58bc104161` |
| `validation/results/transport/v8r2-topas/records/runs/REF-56b37ef1fd8dbf0789c3-ec02c119/artifacts.json` | `8d1ee9002fda6bb8800ad70195639434b5d9aa03e857728f1e51d2a6f80d1a06` |
| `validation/results/transport/v8r2-topas/records/runs/REF-5c309b0ee9905e6be306-efccf482/record.json` | `ac707adc9269594926432aa1495ec5b1c898d8f1daecdbda2cadda3839ec47c1` |
| `validation/results/transport/v8r2-topas/records/runs/REF-5c309b0ee9905e6be306-efccf482/artifacts.json` | `26aa20f5a5a36a1cded5b5f2d1080a84f05613e071915500d5696f30c1b86271` |
| `validation/results/transport/v8r2-topas/records/runs/REF-62dbefc64ddefddef3b4-5724a8ed/record.json` | `80c9e698862321ee060f830125a1073c17be8785fa4c4d8148480c18453df216` |
| `validation/results/transport/v8r2-topas/records/runs/REF-62dbefc64ddefddef3b4-5724a8ed/artifacts.json` | `5ade4e90567da05aeb144ad6b29918620f15fb473127de5ab9d7ad822929f07f` |
| `validation/results/transport/v8r2-topas/records/runs/REF-659638d1daa0c25ab0c4-1f757981/record.json` | `94f1ef124b0f854974a009c2de04e7ea59d9168ab2b2a5f45a459cf4ca1546da` |
| `validation/results/transport/v8r2-topas/records/runs/REF-659638d1daa0c25ab0c4-1f757981/artifacts.json` | `c5a2999b854e7095fb40ff6136d27f5162059d021f80bc9a30b46c2ea735c9bc` |
| `validation/results/transport/v8r2-topas/records/runs/REF-685038154a1f58026833-65882747/record.json` | `dd112f433848f47f13a610e1f940609c63996db0315a5dd4f85ca0633356022e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-685038154a1f58026833-65882747/artifacts.json` | `46d294cf2f84ed744fcc182a7348a599d043078886b2af90a0d2ea13a5244ebd` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6ad649bce2f2c9e9bc80-7b753c4a/record.json` | `39f6ca7441f8c7d971f1ab46bf0fde96453ba540f266dd1fea319c3907902202` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6ad649bce2f2c9e9bc80-7b753c4a/artifacts.json` | `5e235a11aaf427e871c9e259d02ab2d142136825e8a0bb6e943c04b06700ee14` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6b30b8a86a1edc12d348-8bd43e2a/record.json` | `8b73a94eb47bf74ab8c08c9323c85101d9846820389a59e0d426851d8b6ca2ce` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6b30b8a86a1edc12d348-8bd43e2a/artifacts.json` | `5848e8e3c1dcfe9b681a55d9de68823e4c231e6a0ad2d6475e25e2f814a58a04` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6f1fc1d8d8c16b921a7d-267b5a49/record.json` | `2bc93197830b33467e9ffe5438e702cc6a47c920fdfd05bfee5abe0ade0f93d7` |
| `validation/results/transport/v8r2-topas/records/runs/REF-6f1fc1d8d8c16b921a7d-267b5a49/artifacts.json` | `9adb2a1a4f75219b4d83324d5cff8096d94a3f2cccaf345b6d02b75d30fae32e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-72ef20c9ce60ae1dd228-6783b617/record.json` | `b4a46d56b8c6cfbd31be14734e6e3c0bbe8cab304d2fb36df1c7635382c032d2` |
| `validation/results/transport/v8r2-topas/records/runs/REF-72ef20c9ce60ae1dd228-6783b617/artifacts.json` | `a38d6b2be0005114cb2c2ecd1519b8e4809c736f213e5428f446c959bb49c0ad` |
| `validation/results/transport/v8r2-topas/records/runs/REF-750d9a34da3f616e101c-d65c932b/record.json` | `172e67ae7a8be791751829347a2ec455ddcbbcb78107a557c074c452a9062b20` |
| `validation/results/transport/v8r2-topas/records/runs/REF-750d9a34da3f616e101c-d65c932b/artifacts.json` | `86c66c3e1192c9336d2d4048233a3a4bdc9dd1c6a0e4e1901f0b85fc9bef8cc4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-75440efa37bc041bc727-0733d482/record.json` | `39d9943331d4baf04aadce54fb78041e90342f2da9f58edc79609a94d472e021` |
| `validation/results/transport/v8r2-topas/records/runs/REF-75440efa37bc041bc727-0733d482/artifacts.json` | `adbe8d9ab3141a54ff98f052e80f44e262e5179b6d354cd2e024c8e4fabc79b8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-79f238180437968711f6-f2cb2cb0/record.json` | `bb9eafc7472cf1576ffd4573c323e79968065edd3379d60dc2c86eb1736c3111` |
| `validation/results/transport/v8r2-topas/records/runs/REF-79f238180437968711f6-f2cb2cb0/artifacts.json` | `bba447fcc72a8d1a328c7d420ad85689d165b5b08a95e80abb3bef274ebf57f0` |
| `validation/results/transport/v8r2-topas/records/runs/REF-81587d7789d43f89f46f-5bec307a/record.json` | `67e4f596e5b9a20ca941138edfe85e8c49af9e05ed920d5854f2376eb6081f0f` |
| `validation/results/transport/v8r2-topas/records/runs/REF-81587d7789d43f89f46f-5bec307a/artifacts.json` | `68282a92fdf9424091a9a001bedd27a4f387eb443f6f9ffaf76c6ccb87764859` |
| `validation/results/transport/v8r2-topas/records/runs/REF-844b1cf2bfc4b3287667-fc8a609d/record.json` | `c9d836f39bcbc28e83fe623a62b4b1ef43300e099792995cac2296d963de0925` |
| `validation/results/transport/v8r2-topas/records/runs/REF-844b1cf2bfc4b3287667-fc8a609d/artifacts.json` | `c9a7c24f74bb7310b05abec669aa910972a032eba4cdbec5ca0284ee95d59b79` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8516878f8ed4839c0e40-a65debca/record.json` | `de028f4a6f047ec4462dbc0f57a46faf5abcfc5ead616aefd1e313d9b9614a9e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8516878f8ed4839c0e40-a65debca/artifacts.json` | `72a1e977513bcaa790cc6b19f2b3b2836f1a71f7822337d3aff5d685a5e0b283` |
| `validation/results/transport/v8r2-topas/records/runs/REF-85bf230318944e37467a-35b2cabb/record.json` | `83d8b8309f700bbec1d94bb52b6d1cf3dafca7a36edc5f2f9452c4aad4dab363` |
| `validation/results/transport/v8r2-topas/records/runs/REF-85bf230318944e37467a-35b2cabb/artifacts.json` | `cb4b4d6e1e2bacd85353dd9c106f97a0959f45f5ff6cf1c2a66939693e1d977b` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8ab4771dd1c3f998104e-85ff657d/record.json` | `fade3b6c39cf0302c9eb7265a59dc87cdcbe399f9c495a937e26ed9918af1ce1` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8ab4771dd1c3f998104e-85ff657d/artifacts.json` | `f9b554a01a39b10ad21fac2c48a203bd8e9dc961e643b5feb25521e7141cf686` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8b7009c123e89db33318-78bc11f0/record.json` | `558aca769025c6dad6fc920f8ab367a418bb0bfcc30b4f3574dfb4a7feb47120` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8b7009c123e89db33318-78bc11f0/artifacts.json` | `ea06b6779b5e845c8b718ede4fb23c41da152556fa9353605bd81475567bd128` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8e0b95d380d64b45cae0-70144471/record.json` | `f32699eb5770bb19911e9c7579b2670c833c733758e8f8e3a24181e3bd75eca0` |
| `validation/results/transport/v8r2-topas/records/runs/REF-8e0b95d380d64b45cae0-70144471/artifacts.json` | `4200dcacad65f40a3312b8101d14b42e667f4a64a27aee3855202cd331aa88bb` |
| `validation/results/transport/v8r2-topas/records/runs/REF-91a6979b5edd1995642f-93b9881a/record.json` | `fa7f5f9917027b8de963f8929fc770a1b5d36544ee591ed964b5e0e9ac183463` |
| `validation/results/transport/v8r2-topas/records/runs/REF-91a6979b5edd1995642f-93b9881a/artifacts.json` | `8f05dd909569547ce12d8ce5a6f352de421beb75e39323c6d2e8f810199a4681` |
| `validation/results/transport/v8r2-topas/records/runs/REF-933b2e2647ee26272deb-9b8e26de/record.json` | `ea49ad515a92669aee388112b75b9f01d636f71aac614a10128e451949b0a819` |
| `validation/results/transport/v8r2-topas/records/runs/REF-933b2e2647ee26272deb-9b8e26de/artifacts.json` | `65c628a2db6aa2a6be778245f481394d98b2916c3b91571a527edf54a3dda80b` |
| `validation/results/transport/v8r2-topas/records/runs/REF-97f42afba70d58702506-6d3f2983/record.json` | `ed7647e5beb80bd4df35bbfbdd222a7165586f68e67f38d54d65559db7ad4647` |
| `validation/results/transport/v8r2-topas/records/runs/REF-97f42afba70d58702506-6d3f2983/artifacts.json` | `a44010450a58ee8fead836f465250b971d4d7efeb9873e83f46bab7634e1977b` |
| `validation/results/transport/v8r2-topas/records/runs/REF-986a7d3ea8d26f84cd71-e55464e8/record.json` | `d1988264eee3aa3367ed287b6f824b09f60943e54e13ccfc48f65b88f083cd63` |
| `validation/results/transport/v8r2-topas/records/runs/REF-986a7d3ea8d26f84cd71-e55464e8/artifacts.json` | `e6a16c3953284e0cc0cf5098a218033039d122e7cb89dead1f83431ff0be8438` |
| `validation/results/transport/v8r2-topas/records/runs/REF-98e7a0959e5f90db9e8f-4d557788/record.json` | `a5ad5067a270e4baeac140936f246a19f15ff06e7678991fca8f4e7bc205ae9a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-98e7a0959e5f90db9e8f-4d557788/artifacts.json` | `67359267b69fb7dc5c35244efd5941096d71b61e54a4b024d9149c61ca06a7ad` |
| `validation/results/transport/v8r2-topas/records/runs/REF-999bc535b2538ab9591f-0c130fb3/record.json` | `924e9fd802769cd3b3ee2b7be1739c29af1fce03e494780a2fadc12e2710a45e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-999bc535b2538ab9591f-0c130fb3/artifacts.json` | `7ec483ec61aa7d62503276563eb450c21565c2aa71a965fe588abf359b6578c3` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a3a9cbc0bf63b1aa4f16-351e0c17/record.json` | `a9839e4c71a4b353114cef04a740a0073343a042e45e02593662ba5eac01334e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a3a9cbc0bf63b1aa4f16-351e0c17/artifacts.json` | `8baf5df22e751dd6aa233c196c1259374327ca3f1dd22be157b6739b05c1dbac` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a76389b59e272bee7b3d-c398fe21/record.json` | `1646f13e255826e6ded6432ae7c7b47704db4030aade7fafbdb531b5296a8062` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a76389b59e272bee7b3d-c398fe21/artifacts.json` | `45997127d89517e70cd7bc31710f37dc9c10203ec0a91a52dcb6828cf70437aa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a90e8e8eaccf867b5862-6b59925b/record.json` | `64a370d253dbed641b000bf9f2ec46b6a7d387d03c3c474b6ffd1f6cb4eac4b9` |
| `validation/results/transport/v8r2-topas/records/runs/REF-a90e8e8eaccf867b5862-6b59925b/artifacts.json` | `67dddf497312fc248ba6623d84208f2652640f84d8b27742cb03c5f1e645b7e7` |
| `validation/results/transport/v8r2-topas/records/runs/REF-af61a9d385173066f95b-3cf34e11/record.json` | `e978820cfd5c0b2033b0bd2ec4bccc7b08c482eabf7c9312a4d152921c61296a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-af61a9d385173066f95b-3cf34e11/artifacts.json` | `fe6bb9f7429cbf35e20197d93dfc8caeb71d8c4aedcf81ff7c812df335d5c8f1` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b235c3220e66b464555a-85711fcc/record.json` | `4fae9f3222d62941cb684ded92e77a602de9f01276c4b4c2fe54e16fecc9e1ed` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b235c3220e66b464555a-85711fcc/artifacts.json` | `ca357367ab02990b42008e79301b2822c44a22edb2cb01c60024bb9769dda452` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b3289bd56dfb73d3d2eb-87adddc9/record.json` | `add1b441e44e021c669ab4f160f78efbb5705e0748622fc66862bffb735f2124` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b3289bd56dfb73d3d2eb-87adddc9/artifacts.json` | `cef458ba4c62085de2974dea2c31356f3d1ebefd4ea26226e7242958ed15a9e4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b3507ec2481cbf21ebce-65cedd9b/record.json` | `2bb5944a0a4a62e3e2d85f83981b3bf26bde5bbb690254ad5de32c7bab814aab` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b3507ec2481cbf21ebce-65cedd9b/artifacts.json` | `5dfc3aec39b6da135dfb0aa2067e282599f45e62a3fd0da04769b91ced59c328` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b5b9ff08f408b4ca5317-7d7f168d/record.json` | `80a92e0172bfe43c2c5afd4deffcb0db5a652061194f4cb26ed737cbd68d07a4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b5b9ff08f408b4ca5317-7d7f168d/artifacts.json` | `17dc0b66a5e88877e2baed065893499b0367c87cccc5c3ea9ccc2cae7d1eaba0` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b73b308392747b407aba-884b2a64/record.json` | `9699f569714661aa1784eb3482ef172a2d227b4004c1972a573192f2260917d5` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b73b308392747b407aba-884b2a64/artifacts.json` | `9d460217e7922b1fb2f0e19a76284652d483db9fad1e8fb1e7ef279e45cbbaa8` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b8295b6f8d4592fa533e-b95784a2/record.json` | `408d7aea22369d9f6b5bc92210fb9afdd0afb2efa31245d5bda805e6e3e93441` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b8295b6f8d4592fa533e-b95784a2/artifacts.json` | `c245679e08e610913ba6cd6afd7ea6024d89add5f027e8a4fa4840cb96ae23df` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b8ecd7fdcf2b60c4f22a-b5a017b7/record.json` | `a460f3971648d13ba491b2626382f0bf3ab108d713468722d06ffaeb9ed05635` |
| `validation/results/transport/v8r2-topas/records/runs/REF-b8ecd7fdcf2b60c4f22a-b5a017b7/artifacts.json` | `702dc821e027faf047c007becb7e9a2c3af66ec18695a9e46e19a38be46c82be` |
| `validation/results/transport/v8r2-topas/records/runs/REF-badd5762e53b6bd276e3-21fd2fff/record.json` | `f1eead5c815b7f270e7712909cddcf6bc28926d537632af6bcbea5686193626c` |
| `validation/results/transport/v8r2-topas/records/runs/REF-badd5762e53b6bd276e3-21fd2fff/artifacts.json` | `ad1715897b815bb224169d911bc34a25fa7f2f30784517550d0081aab559c1da` |
| `validation/results/transport/v8r2-topas/records/runs/REF-c040317b99df027519b9-83b73a0f/record.json` | `cb701d04b9f928ee2924d2fd1cd857a7823d1da8306bb480f1a32085c5c40f43` |
| `validation/results/transport/v8r2-topas/records/runs/REF-c040317b99df027519b9-83b73a0f/artifacts.json` | `fad51c899e348d1177ec8f1f8d5f67a5995cad8f8306c650c30657a18160e32d` |
| `validation/results/transport/v8r2-topas/records/runs/REF-c431be25d1b37248cd7e-4bdb3417/record.json` | `bd3a5a3d0cde8c4f70d852d478721c24858c0bb7ad3bc6d701be414137ada1b5` |
| `validation/results/transport/v8r2-topas/records/runs/REF-c431be25d1b37248cd7e-4bdb3417/artifacts.json` | `75cc2d0411bc7c77ac34cb79d61aea8f6a7fac9c5007918e72288bab5872e226` |
| `validation/results/transport/v8r2-topas/records/runs/REF-d16fbbeb28af18b513b1-eedfbfd8/record.json` | `2ed2803a938ef4e51b5ab3cf605e01c32e92656137fcac80cf38bc786a0faaaa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-d16fbbeb28af18b513b1-eedfbfd8/artifacts.json` | `ce3f086e2ecccfceaa9e071a225dcdd009e7c726ff9c83a65aa9c69e438cf6ba` |
| `validation/results/transport/v8r2-topas/records/runs/REF-ddb4748e2128aaa99e06-48894b44/record.json` | `04154c327af4e3034ae329f7b121ed4fb267875656d959c112b3fd81e7d77fde` |
| `validation/results/transport/v8r2-topas/records/runs/REF-ddb4748e2128aaa99e06-48894b44/artifacts.json` | `e2c1e1f6bb0cf5ce8ae382335ad2cbc01a6a6e451f234a0d5d955da1a13eb1f2` |
| `validation/results/transport/v8r2-topas/records/runs/REF-e244c0801114dabea9df-af8f93df/record.json` | `30fc6e48d3014e2fcbc0e5fed8bf1777094c4a3c3fa092313bda1476f1e36f9e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-e244c0801114dabea9df-af8f93df/artifacts.json` | `0647735df3137dd748b693f6ae2f28e8a2618bfb6169f622e3b360c45118d229` |
| `validation/results/transport/v8r2-topas/records/runs/REF-e6de6cb81a11a5af603b-264be1c5/record.json` | `6749fa90e457d08e0f794859e59f15c6bbc4be935809da7be66afa630a2ce59e` |
| `validation/results/transport/v8r2-topas/records/runs/REF-e6de6cb81a11a5af603b-264be1c5/artifacts.json` | `77375b8144d3a7d2bf1772249773e403d5452484407a8cd64ea9de00e5bffcaa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-f085e4479d78b277a4c6-707fcf44/record.json` | `a6231cc46bcc6725f825a75966c7e0ff16912cc313f6e0f056b184c549476e79` |
| `validation/results/transport/v8r2-topas/records/runs/REF-f085e4479d78b277a4c6-707fcf44/artifacts.json` | `5ebedee6992723a1c2b66fdaa14b7e5055913a8e5461d1958a699b8ebd2b720a` |
| `validation/results/transport/v8r2-topas/records/runs/REF-f596373c28b3d26ab5d5-a88b001e/record.json` | `86d846dee1afebd7d72aac1a47e5a526a092e4a834fc043dd715fad211c1c2f4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-f596373c28b3d26ab5d5-a88b001e/artifacts.json` | `7fa5441703d6ef0a3d33597ebf2a94dd6fb68ce3877fbfae5cc0897c8c92c3fa` |
| `validation/results/transport/v8r2-topas/records/runs/REF-fc2f8c0c67ddd847ed75-f40e519a/record.json` | `2aecbc913e9c0ce4cb92d38acec6f5c0184ef0fcf50257d91ac6c7ff998900a4` |
| `validation/results/transport/v8r2-topas/records/runs/REF-fc2f8c0c67ddd847ed75-f40e519a/artifacts.json` | `a1c4fbd7c0ed52251e9b6ac05c9ac7a75df7f6b6c7aa4ed300fe57c5836e16d0` |


**Dated note (2026-11-09, orchestrator) — records/index.json filled with the 72 EM evidence parts (IDD-EM 0–3, LAT-EM 4–7, cases 1–9); entries that changed:** `validation/results/transport/v8r2-topas/records/index.json` 2e092ee3… → 264fba78…. All other entries unchanged (recomputed).

**Dated note (2026-11-09, orchestrator) — files added by the EM results commit (the 18 TD-EM/TL-EM grade reports of cases 1–9, graded once at evaluator HEAD 94bc948e; final, never re-graded; sha256 as committed):**

| File | sha256 |
|---|---|
| `validation/results/transport/v8r2-topas/reports/TD-EM-case1.json` | `d7e2f1c6f2bcf988c8ebf0b87f981b3ebed7e28f8b7f98f098157a620495841f` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case2.json` | `61c5dbb7d4dd8e8e5c275dc055bdef526e047e9e183d68c77c5dbefdfab3dd2d` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case3.json` | `d7a160df3b7c5a1a0bd7abe074890dec5e7590170c0b020615f1c94f52d914e4` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case4.json` | `9ea779c5377cc3432e669c5e6332b6ab178b66254c37e661812526365e094ac2` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case5.json` | `072054ddfa9eaabf8766a2a365789baade1c5f42f544098a1f425a01ad6b4949` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case6.json` | `f8cfc3cb0e962d2cf2dc60178efce0d5666cb00bc779ff914299d7d53d7fc14b` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case7.json` | `42eefce83f1345085ddbd7543709039cd8739bfb7d2e198d9269c608897b5038` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case8.json` | `70475da41a29b56969cad08b31203bee779d58e270a9375e213871c1e1b25667` |
| `validation/results/transport/v8r2-topas/reports/TD-EM-case9.json` | `4f4f670ee978456ba883b946b9a34bff2ba6ba648a9e205215c0a119391ab0bd` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case1.json` | `7c6726fc1b7df0dd744492fdb31e618e447cbee9a2192c2cafc9909ab7baed3e` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case2.json` | `159f6442979caa8286963ad02b08fc47f74eee7755a9f8ae6d30c4a3df58abb6` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case3.json` | `aa135a48e52d0e56e221f6849a369da319f77325d1282cec06bb21d0a9987c35` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case4.json` | `a1eb0d27f62c6f97445de1856657e21619054843fe1dbbcf3eff328d5e13f9a3` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case5.json` | `1a4bc36ff47461cee37e8499c901f1c2c73caee269a4f3544690e0288dd081b9` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case6.json` | `98e210e19333b8e4d8aaa2e5a13c94cca290cea15b1f2c5f2c294f4dbd09e553` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case7.json` | `032f1307c69c537349bb0f39ca297be0381814f6caa6ccc2284c488e743fedfa` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case8.json` | `66b4aea4ef24247c34cb44a761888c97254809cbda233ed9a4c1d7b48d9216a1` |
| `validation/results/transport/v8r2-topas/reports/TL-EM-case9.json` | `b8d07ea7b383b69e4da2f5e3c18c15006d8c02892f7605d05e1b068db16f18b5` |

## 3. Acquired and pinned external sources (hashes as recorded in the P-records; plan §G)

Roles: exploratory / lineage unless stated. Nothing below is redistributed in Git.

| §G | Item | sha256 | Record |
|---|---|---|---|
| 1 | [URB] `G4UrbanMscModel.cc` v11.4.2 | `22c5f52b3a1cfc8dd99dfd07929bed702cb2a35a0d3a7c727450bfbfc9b69558` | P2 §0 |
| 1 | [LS] `G4LindhardSorensenIonModel.cc` | `e16ed6cd1d132183684ef47948b6f8fe36526a48663a6dc1878d5bb8330eb12c` | P2 §0 |
| 1 | [OPT4] `G4EmStandardPhysics_option4.cc` | `04337300cc119a4fd893d515059cc75cef9f042849b183848643a4c4bf104c89` | P2 §0 |
| 1 | [ICRU73] `G4IonICRU73Data.cc` | `cf771a5a44ea936c56593b7dd2d67bb4157345eb316815f60880986404fe08d2` | P2 §0, F0 |
| 1b | [EMB] `G4EmBuilder.cc` | `40782f5c5bc0e66f9b4b635a08d8b0c236b8576ad5b885afca04b679bb18d3d1` | P2 F0 |
| 1b | [HMSC] `G4hMultipleScattering.cc` | `a664386984f22898ca8bcaf0f38326479e1e463de88446d41fd465dac83f50ea` | P2 F0 |
| 1b | [VMSC] `G4VMscModel.cc` | `87b9f203542f2412b3bf327b2cc41eb0c5295a92fabf703009b66f7d37b97a69` | P2 F0 |
| 1b | [IPM] `G4IonisParamMat.cc` | `990c6f315a3b6951253eb251024ba439f8edf4e151994c33991510645807849d` | P2 F0 |
| 1c | [SS] `G4EmStandardPhysicsSS.cc` | `8aeaf91ce73c520e09e482d64a2e6a3b6c090124a9c14451d07b0b6eed9d9767` | P10 §7.1 |
| 1c | [CS] `G4CoulombScattering.cc` | `db354d80c0a6efd278abf7793144bac1e5f3939a8d842ccd123014ea9ca0b754` | P10 §7.1 |
| 1c | [eCS] `G4eCoulombScatteringModel.cc` | `fa5ddeb923b221c6c086593d0a006fd1f21502b5d460692117ba7e52e2f8c9af` | P10 §7.1 |
| 1c | [WOKVI] `G4WentzelOKandVIxSection.cc` | `a64c94b1adf953beba1289575134514ca1261587e06b15b3286461418339a24c` | P10 §7.1 |
| 1c | [WOKVI-h] `G4WentzelOKandVIxSection.hh` | `3a3f1db069d55ae06bd1c2e84c84c15a64dac2362371ed4a7dc783362184090d` | P10 §7.1 |
| 1c | [ICS] `G4IonCoulombScatteringModel.cc` | `ef9b31d40e99d72545dd06d329e13d7453bc20dea0a7681ce6ecc929e28d9296` | P10 §7.4.1 |
| 1c | [ICS-h] `G4IonCoulombScatteringModel.hh` | `a1a9a4f5eb1dd0f323e5eac08873f435e87b739ff705c5b441caf4bb2262e4c0` | P10 §7.4.1 |
| 1c | [IXS] `G4IonCoulombCrossSection.cc` | `87989b93748a8026c7b7031098166fa24c2f6d9c00979671bd0d9649e321e2f4` | P10 §7.4.1 |
| 1c | [IXS-h] `G4IonCoulombCrossSection.hh` | `cd5d0b5cbb4dbc2963ededc5519d110ac1dd9042d38a9a8a0742051c46928b84` | P10 §7.4.1 |
| 1c | [EMP] `G4EmParameters.cc` | `aaa143c1cb34b544f0d67fa6a7242d2b3fbd86ebfeb1bf4bc5c9fcca2973909e` | P10 §7.4.1 |
| 1c | [eCS-h] `G4eCoulombScatteringModel.hh` | `68ce763e099bc31a8dfc67609ad39dcb7c070719a3a446f3971ed044255f109e` | P10 §7.4.1 |
| 2 | G4EMLOW 8.8 archive (lineage only; Geant4-restricted use) | `b60cfd63176f5d16107e2a25b35b235155032d1735d749670ca50fede12624cf` | P5 |
| 2 | member `icru73/z6_G4_WATER.dat` | `444be013f716f5701051a5d4acf1a5ee7ebbdff7b6e2ebc78ca07563f9820159` | P5 |
| 2 | member `icru73/z8_G4_WATER.dat` | `6d8255d47aca417bd10d19684350203c0276d3df828886fbdccaea8cdb19bd67` | P5 |
| 2 | member `icru90/z6_G4_WATER.dat` | `d9d4809dcbacdb74e5e6c069709ef848c9d715b097f4ff4db24ef9969e727d75` | P5 |
| 2 | member `icru90/z8_G4_WATER.dat` | `863898af3649947df532c6eeb244b2658503c0c7fa3882f81a504278793717f6` | P5 |
| 3 | libdEdx archive (GPL-3; I3 source, exploratory until frozen, evaluation after) | `b941a0c30793daddc4d119a2cc664bfbdd430554c3a4788be04e97ff6fcc0863` | P5b |
| 3 | member `data/raw/ICRU73_NEW.dat` | `749768c1c2de475c928cc2c50a354903577139c11debd4b5948b035ca3a3d87e` | P5b |
| 3 | member `data/raw/icru73_newEng.dat` | `8b5b809f5c0cc4f6018f9faa89e9f73c7b9b81a3c588f98bef556712d4c77dd1` | P5b |
| 3 | member `src/data/embedded/dedx_icru73new.h` | `f6d088b233e6584363ee09db6b6d828dde0a3c2fe62850e2a35285e83ff78459` | P5b |
| 4 | Sihver, Schardt, Kanai 1998 PDF (FAILED as a V3 reference) | `68567feb58bf3c2a669193e87d1f0ce41f86c21dcd36959f5ffd7c390bfcda55` | P6 |

**Dated note 2026-11-01 (FC-2; `p2b-follow-ups.md` §0):** acquired sources added.

| §G | Item | sha256 | Record |
|---|---|---|---|
| 1d | [BB] `G4BetheBlochModel.cc` v11.4.2 (acquired, exploratory, 2026-11-01; same hash as the earlier scratch download listed below) | `5dc2a9002bff25f239a6402d3976a1f285f071abf4b41585b76f275e4539e01a` | `p2b-follow-ups.md` §1 |
| 1d | [BB-h] `G4BetheBlochModel.hh` v11.4.2 (acquired, exploratory, 2026-11-01) | `cf05e98c6b56aa07e227800ccff169b74e2c724d56ce5753fcab886be8fac58e` | `p2b-follow-ups.md` §1 |
| 1d | [NMB] `G4NistMaterialBuilder.cc` v11.4.2 (acquired by the orchestrator, 2026-11-01; lineage) | `a1f8944dd63f4c9915c2205ee2ca5ec0ffade15ac0ca8afe38d53597c2671abc` | `p2b-follow-ups.md` §2 |

Also referenced, not acquired through `acquire_reference_data` (scratch downloads, lineage only; recorded so the
freeze names every source a record relies on): `G4BetheBlochModel.cc`
`5dc2a9002bff25f239a6402d3976a1f285f071abf4b41585b76f275e4539e01a`, `G4EmCorrections.cc`
`a357d07617816f4e3ef1c7558890cf4966ac55bf85ba4ad2392713b10b52d3c7`, `G4ionEffectiveCharge.cc`
`891854a97b97f038f9a0fa88b9f2f4ed139ea080e8edf9bd349cfd93031ae4c5`, `G4LindhardSorensenData.cc`
`70ae7c6ba3ef990b7005240bd2b60e4040a3fb60bb1c9e4dbf829c109180eda9` (all `p3b-ls-default.md` §2); the Tessonnier 2017
thesis `70df35938f57acbf6d88922593e8cb9612a0e0f83701499f181244fc9850428a` (P6-lateral). The Gottschalk arXiv:0908.1413
object used by the P1 proton validation is recorded only by its prefix `67fb1478…` in the V3-008 records (full hash
not established by the planner; the P1 tests resolve it from the cache). **Dated note 2026-11-01:** `G4BetheBlochModel.cc`
is now acquired (row [BB] above) and is counted there, no longer among the scratch-referenced files.


**Dated note (2026-11-02, FC-1 provenance) — inputs pinned for the committed P1 completion calculation (evaluation role):** arXiv:1610.01279 PDF `934fd31644eb97734a0f4c9a756a40d49012d3c126ef2b4391e088be7a092aa4`; NIST PSTAR POST responses Be `0fcbb24d9b7f07af7898ac3954802d08010d22339e195ee7cd18fe2638380757`, Al `3fe6459ec4abde6052bcb63f427ea856bbc35dd8506014191145b4eb9253db5e`, Cu `986e8577db37120c405f6686394c1250c3e7f1d6af7e88547b964de2e99e0063`, Pb `11dbcd1b014d926ef3e5b8bdd2bf9e539c153b5a4c8d51a62d81c34f5b0a05a7`, matno 219 `8974efe7242af76f1df4fb696a7b957c696f6ca7de711e01c3c68a2f40c04040`, 223 `1e13828ffc3dee1e5748d2d922e5ecf2c583380903bb787f84fbe2886dc549e9`, 226 `dc65ded99e0cb9b455482404642054c628f0639fd75ed8ec5ebc0e5e68f981c1` (POST body `prog=PSTAR&matno=<n>&ShowDefault=on&NumofEnergies=0&character=space&electronic=on&nuclear=on&total=on&csda=on&project=on&detour=on`).
## 4. Seed reservation (plan §C; P7)

- V3-008 block 20572004–20632003 (`V3_008_RESERVED = range(20572004, 20632004)`; ledger commit 06974d07 on the V3-005
  branch). Guard test `tests/ionmc/test_v3_008_seed_families.py` (sha256 in §2.2).
- lv8 base 20581004, hr8 base 20591004; row seed `base + 100·r_index + k`.
- Rehearsal family 2060xxxx; consumed so far by the ten smokes 20600401–20600410 (execution only; never evidence).
- TOPAS seeds `20619000 + 100·part + case` (case 1–40, part 0–19). V5-SS cases 22, 25–31 (branch B), parts 0–1 under
  the 6e4 design; cases 23, 24, 32, 33 retired unused.
- Re-run bases 20621004, then 20631004, each declared in an amendment before use.
- **No V3-008 qualification seed has been consumed at `f694b747`.**

## 5. Results status at `f694b747`

- **No transport or qualification row result exists**: no ionmc transport run, no LV/HR qualification step, no
  scored reference-engine (TOPAS) run, no V-row run. The only reference-engine runs are the ten scorer-free smokes
  of §2.5.
- Pre-freeze SEEN computations exist and are disclosed (plan header; decision 0043 §3 "Freezing"): P1/P1b tables;
  P3 and P3-mixed residuals; P11 I11-rel/I11-abs (070b07f0, regenerated 00af9a6c), P11b/I11-ext (b3eef90b); the
  V5-SS planning budget, bias, discrimination, power and cost (`ss_budget.py`, P10 §7.5.1–7.5.2, §7.8).

## 6. Pre-registration status per row (decision 0043 §3 "Freezing", corrected text, updated to `f694b747`)

**May claim pre-registration (values unseen at `f694b747`):**
- capability and implementation: I0 except the sub-items listed below, I1 except the mapping callables, I5, I7, I16,
  C-res, T1-ion, T12-ion, T13-ion;
- theory identities: I4, I4b;
- stopping comparisons: I3, I3b (i), I3c;
- transport: I6, I8, I8b, I9, I10, I12, I13, I14, I15;
- reference rows: V1, V1-held, V2 (Case C, report-only), V4, and the V5-SS statistic (outcome unseen; its budget,
  bias, discrimination, power and cost are SEEN planning values);
- R2 (not yet evaluated against an ion-enabled build).

**Carried as disclosed SEEN, not pre-registered:**
- I11-rel (gate, 2 %, tolerance unchanged), I11-abs and I11-ext (report-only), I11-map (identities ran at 070b07f0;
  the src implementation tests ran at 683480b/f694b747);
- I2 (regression on seen table values), I3b (ii) (the P3 curve; the three LS variants are compared with each other
  in `test_ion_em_ls_option.py`, not with the I3 reference);
- the P1/P1b reference tables;
- **new at `f694b747` (planner, 2026-11-01):** the implementation tests of 683480b ran before the freeze, so the
  following sub-items are SEEN implementation checks, not pre-registered: in I0, the `fdm_domain_status` fields, the
  `fdm_extrapolated` flag at the standard/end cases, the fail-closed mapping-edge rule (4000 MeV) and the
  `ls_requested`/`ls_effective` persistence; in I1, the twin-vs-Warp CPU f64/f32 parity and the completeness check of
  the exact-Ω mapping callables. The remaining I0/I1 sub-items (species, other fail-closed rules, the
  charge-argument functions) are unseen.

## 7. Acceptance tests of the freeze commit (to be run by the orchestrator, recorded in the freeze review)

1. Every SHA of §1 is an ancestor of the freeze commit.
2. `git status --porcelain` is empty; `git diff --name-only f694b747 <freeze>` obeys the rule of §1.
3. `sha256sum` of every path in §2 equals the recorded value (53 repository paths; count below).
4. FC-1 and FC-2 are disposed by committed records or dated amendments, named in §0 by dated notes.
5. The logged seed-family test passes at the freeze commit
   (`tests/ionmc/test_v3_008_seed_families.py`).
6. No file under `validation/generated/`, no LV/HR record and no scored REF run for V3-008 exists at the freeze commit.

**Count.** Repository files bound by sha256 in §2: 6 (2.1) + 7 (2.2) + 5 (2.3) + 14 (2.4) + 21 (2.5: 20 case files and
the README) = **53**. External sources with full sha256 in §3: 29 acquired + 5 scratch-referenced = **34**. Total
full sha256 values in this record: **87**.

**Dated note 2026-11-01 (FC-2):** repository files bound in §2: 53 + 1 (`p2b-follow-ups.md`) = **54**; external
sources with full sha256 in §3: 32 acquired (29 + [BB], [BB-h], [NMB]) + 4 scratch-referenced = **36**; total full
sha256 values in this record: **90**. Test 3 above checks the 54 repository paths. Any record or amendment that
disposes FC-1 adds to these counts by its own dated note.

**Dated note (2026-11-02, orchestrator) — counts at the freeze commit:** repository files bound by sha256: 54 (FC-2 note) + `p1c_attribution.py` + `test_ion_em_p1c.py` + `p1c-aggregates.json` + the two intervention record files = **59** (every path recomputed by the orchestrator's script at this tree; `p1c-completion.md` and `p1c-aggregates.json` carry their final hashes after the JSON regeneration). External sources with full sha256: 36 + 8 (FC-1 inputs, §3 dated note) = **44**. Total full sha256 values: **103**. Item 3 of §7 reads "59 repository paths" accordingly.

**Dated note (2026-11-07, planner; Codex REVIEW-ad85dff1 bookkeeping finding (a)) — count correction.** The orchestrator's
note of 2026-11-02 above omitted `validation/plans/v3-008/p1c-completion.md` (added to §2.4 by the FC-1 dated note).
Repository files bound by sha256 at the revision-1 freeze commit: 53 + `p2b-follow-ups.md` + `p1c-completion.md` +
`p1c_attribution.py` + `test_ion_em_p1c.py` + `p1c-aggregates.json` + the two intervention record files = **60**.
External sources with full sha256: **44** (unchanged). Total full sha256 values: **104**. **Item 3 of §7 reads "60
repository paths"** for the revision-1 freeze commit (superseding "53" and the 2026-11-02 "59"). Nothing else changes.

## Revision 2 additions (2026-11-07, planner, at `d082960169b00960f7a192ce6e6e3a44edf65581`; dated note, earlier text kept)

**Rule.** Plan revision 2 (`validation/plans/v3-008-acceptance.md`, status block and Amendments 9–15) and decision 0043
revision 2 are **the frozen revision** once committed by a reviewed commit (the "revision-2 freeze commit"). Revision
1's freeze (the reviewed freeze commit d0829601 and its predecessors) is superseded **before any result**: at
d0829601 no V3-008 ionmc ion-transport run, LV/HR qualification step or scored ion TOPAS run exists. Everything
bound in §1–§4 above stays bound (hashes unchanged except the paths rebound in R2.2 below: `validation/plans/v3-008/p10-single-scattering-availability.md`, `validation/plans/v3-008/p11-exact-omega-mapping.md` and `validation/reference_cases/README.md`, and, from round R2-TOPAS-4, the intervention `README.md` (status-update section appended; its revision-1 entry below is superseded history, current entry in the §2 dated note of R2-TOPAS-4); the intervention `request.json` additionally carries its own current entry in the §2 dated note of 2026-11-07); the revision-1 history is
kept.

**R2.1 Commit lineage.** Required ancestor of the revision-2 freeze commit, in addition to §1:
`d082960169b00960f7a192ce6e6e3a44edf65581` (revision-1 freeze commit). Clean tree at commit and review. `git diff
--name-only d0829601 <revision-2 freeze>` lists only: the plan, decision 0043, `review.md`, this file,
`validation/plans/v3-008/r2-topas-target.md`, `validation/plans/v3-008/p10-single-scattering-availability.md`,
`validation/plans/v3-008/p11-exact-omega-mapping.md`, `reference_requests/operator-directives/2026-11-07-topas-validation-target.md`,
the intervention resolution files under `reference_requests/interventions/IR-20261102-010435-B59B2D/` (of which
`response.json` is bound, see R2.5), `validation/scripts/ion_em/gen_r2_bundles.py`, `tests/ionmc/test_r2_bundles.py`,
`validation/reference_cases/README.md`, and the `case.json`/`input.txt` files of the 205 directories
`validation/reference_cases/topas/*-r2-*` (410 files, of which 26 files, the 13 smoke bundles, changed in
revision 2 round R2-TOPAS-2 and restored to their b8c9836f content in round R2-TOPAS-3). Any other path is a STOP for the revision-2 freeze. **Dated note (2026-11-07, Codex
REVIEW-7a7186a2 finding 1):** this whitelist is exactly the set of paths changed relative to d0829601 that are not
the intervention README/request files already bound at revision 1; the three paths added here
(`gen_r2_bundles.py`, `test_r2_bundles.py`, `reference_cases/README.md`) were missing from the earlier text. **Dated note (R2-TOPAS-4, Codex REVIEW-e637178c):** the intervention `README.md` is now also changed relative to d0829601 (status-update section appended, original lines verbatim) and is bound at its current hash; the 26 smoke-bundle files changed again in round R2-TOPAS-4 (accurate post-seal archival wording; the 192 evidence bundles stay byte-identical).

**R2.2 Files added or rebound (sha256 after pre-commit, 2026-11-07).**

| File | sha256 | Note |
|---|---|---|
| `reference_requests/operator-directives/2026-11-07-topas-validation-target.md` | (single entry: §2 dated note of 2026-11-07) | the operator's directive, verbatim (new) |
| `validation/plans/v3-008/r2-topas-target.md` | (single entry: §2 dated note of 2026-11-07) | revision-2 design record with appendix A (new) |
| `validation/plans/v3-008/p10-single-scattering-availability.md` | `a7f539a0b71a39939aed2644a9d7c2969133251ceb1676a3e48572a9292f7fda` | rebound (was bea86277…; dated note appended) |
| `validation/plans/v3-008/p11-exact-omega-mapping.md` | `b40530f48391d6eca98656d5bdb06a5e4272bf22ee0a469e08357ebcd457e330` | rebound (was 6116648a…; dated note appended) |
| `validation/reference_cases/README.md` | `e16cc9d3a2eecf2717e405ae2c4fc53c83b884f37af64c758d6ac6bddb4a482e` | rebound (was `2c160fb5…`, kept in §2.5 as superseded revision-1 history; r2 layout entry added) |
| 205 bundles `validation/reference_cases/topas/*-r2-*` (`case.json` + `input.txt`, 410 files) | bundle-set digest `a8a1c27d0177d38e763344248753f0cb9792bbf2277158622b41bcbf9d8f423a` | per-file sha256 in `r2-topas-target.md` appendix A; digest = sha256 of the appendix listing block |

Each of the 205 `case.json` files passes `infrastructure/experiment_v3/reference.py::validate_case`; the seed equals
the native `i:Ts/Seed` line; `histories` equals `NumberOfHistoriesInRun` (checked 2026-11-07). 192 evidence parts
(seeds 20619001–20620512 by `20619000 + 100·part + case`, cases 1–12, parts 0–15) and 13 execution-only smokes
(rehearsal seeds 20600411–20600423).

**R2.3 Seeds.** Unchanged formulas and guard test. TOPAS cases 1–12 re-assigned (none ever ran); part blocks 0–3 /
4–7 / 8–11 / 12–15; parts 16–19 reserved for the FULL sizing fallback; V5-SS cases 22, 25–31 unchanged; 23, 24, 32,
33 retired; 13–21 and 34–40 unassigned. ionmc TD/TL samples: lv8 r_index 30–41. **No V3-008 qualification seed has
been consumed at d0829601.**

**R2.4 Pre-registration status (revision 2).** May claim pre-registration (values unseen): TD-EM, TD-FULL, TL-EM,
TL-FULL, TD-H, TL-H, TF-DIAG, in addition to the §6 list (V1, V1-held, V2 and V4 are superseded, never run). Disclosed
SEEN inputs to the revision-2 tolerances: decision 0038 stopping comparisons (0.18 % / 0.10 %), the V3-003 proton
EM-vs-full R80 shift (−0.2 mm, `validation/results/reference/README.md`), the P2 §6.3 Urban scratch estimate, the P11b
domain status. The planning CSDA ranges of `r2-topas-target.md` §3 are project-table construction values (box and bin
sizing), not reference values.

**R2.5 Counts (counting rule: one count per distinct repository path that has a current full sha256 entry in this
record; a rebound path counts once, at its current entry; superseded hashes kept as history do not count; the
bundle-set digest is not counted).** Repository paths bound at revision 1: 60 (README, p10 and p11 included).
Added by revision 2: the directive (1), `r2-topas-target.md` (1), `validation/scripts/ion_em/gen_r2_bundles.py` (1),
`tests/ionmc/test_r2_bundles.py` (1), the intervention `response.json` (1), and the 410 bundle files (410) = 415.
Rebound, not added: p10, p11, `validation/reference_cases/README.md` and (R2-TOPAS-4) the intervention `README.md` (its revision-1 entry in §2 is superseded history and counts once, at its current entry). Total repository paths: 60 + 415 =
**475**. External sources with full sha256: **44** (no new acquisition). Total full sha256 values: **519**, plus the
bundle-set digest. The 26 smoke-bundle files changed in round R2-TOPAS-2 were restored to their b8c9836f content in round R2-TOPAS-3, and changed again in round R2-TOPAS-4 (accurate post-seal archival wording; the count of paths is
unchanged; their hashes and the digest are the values of `r2-topas-target.md` appendix A, bundle-set digest `a8a1c27d…`). Each of gen_r2_bundles.py,
test_r2_bundles.py, `response.json` and the directive has exactly one current hash entry, in the 2026-11-07 dated
note of §2 (the R2.2 table rows of the directive and `r2-topas-target.md` carry no hash of their own).

**R2.6 Acceptance tests of the revision-2 freeze commit** (orchestrator, recorded in its review): (1) d0829601 and
every §1 commit are ancestors; (2) clean tree; the R2.1 diff rule holds; (3) `sha256sum` of every one of the **475** repository paths that has a current entry in this record (§2 including its
dated notes, with p10/p11 at their current values and `p1c-completion.md` at a1a24cd3…, `validation/reference_cases/README.md`
at its R2.2 rebound value, and the 410 bundle files at their appendix-A values) equals the recorded value at the
freeze commit; the appendix-A block hashes to the bundle-set digest `a8a1c27d…` (also printed by
`gen_r2_bundles.py --check`); every path changed relative to d0829601 is in the R2.1 whitelist and every
whitelisted path is bound by a current hash entry in this record, with exactly these exceptions: the plan, decision 0043, `review.md` and this record itself (bound by the freeze commit SHA; a file cannot contain its own hash), `response.json` and the intervention `README.md` (rebound in R2-TOPAS-4 to its current hash after the appended status update) are bound and counted (R2.5), not excepted; (4) the logged
seed-family test passes; (5) `validate_case` passes on all 205 bundles; (6) no file under `validation/generated/`, no
LV/HR record and no scored V3-008 ion REF run exists.

**Dated note (2026-11-08, Amendment):** mapping evaluation = tabulation (decision 0043 section 3, "Evaluation", the allowed alternative), chosen before any revision-2 ionmc sample or TOPAS evidence run; reason: measured cost. Table: `pv_eq/pv` per (material, ion) on a uniform 128 x 128 grid in (ln pv, ln(p1v1/pv)) (a shear of the decision's (ln pv, ln p1v1) plane), nodes computed by the unchanged root solve `pv_equivalent` at validation (only the nodes of cells meeting the band of the run's initial states pv(E0 - 6 sigma_E) ... pv(E0 + 6 sigma_E); every other node invalid, a state there fails closed as `mcs_mapping_failure`), bilinear interpolation (`exact_omega.pv_equivalent_table`, shared by the python twin and the Warp kernel); option `PhysicsOptions.mcs_mapping_eval` = "table" (default for z >= 2) or "root" (the per-step root solve, retained); the summary records the evaluation, table sizes, axes, domain and the sha256 of the table rows. Bound test `tests/ionmc/test_ion_mcs_table.py::test_table_variance_error_bound` (He, C, O in water and skeletal muscle; >= 2000 random states plus edge states per ion and material, upper and lower domain end, thin-slab edge included; bound |var_table/var_root - 1| <= 1e-3): measured maximum 2.05e-5 (C in skeletal muscle; He 4.4e-6 / 5.5e-6, C 1.83e-5 / 2.05e-5, O 1.84e-5 / 2.04e-5 in water / muscle); transport-level counterpart `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states` max 6.8e-6. Measured MCS cost per step (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 0.16 us (root 39.85 us), C 150 MeV/u 0.17 us (root 39.57 us). The model (decision 0043 section 3) is unchanged; the root solve remains available as an option. The rebinding of the changed bound files and the binding of the new implementation paths are recorded in the §2 dated note "2026-11-08, post-freeze implementation commits". (2026-11-08, Codex REVIEW-d3a17ebf) The bound is in addition certified deterministically at build for every produced table (`exact_omega.certify_pv_eq_table`: the root solve at 5 interior points per solved cell, the centre and +-1/4 of the cell in both axes; per-step variance of the first and later steps; a table above TABLE_VARIANCE_BOUND = 1e-3 fails closed; the maximum, cells and points are recorded in the summary), measured maxima at the upper domain end (128 x 128, 628 cells, 3140 points per table) He / C / O: water 4.11e-6 / 1.97e-5 / 1.98e-5, skeletal muscle 4.02e-6 / 2.20e-5 / 2.20e-5, over all eleven built-in materials and a user-defined material at most 4.44e-6 / 2.87e-5 / 3.21e-5 (`test_ion_mcs_table.py::test_certification_every_builtin_and_user_material`). (2026-11-08, Codex REVIEW-9529fd35) The five-sample certification is replaced by a domain-wide analytic bound per solved cell (`exact_omega.certify_pv_eq_table`): the tensor-product bilinear error bound hx^2/8 sup|q_xx| + hy^2/8 sup|q_yy| of q = pv_eq/pv, second derivatives from second differences of the root solve on a 2x-refined grid (3 x 3 stencil per cell, maximum over the cell and its solved neighbours, which covers cells where a second derivative changes sign), times safety factor 2 (q is smooth in ln pv and ln r, P11; the factor covers the variation of the second derivative within a cell; empirically, over the 70 tables below the sampled cross-check passes everywhere with factor 1 and fails in every table with 0.75, so factor 2 leaves a per-cell margin of at least 2), converted exactly into the relative per-step variance error (only pv_eq differs between the paths and f_dM is linear in log10 pv_eq at fixed memory term l1: var_table/var_root - 1 = (c2 + c3 l1) dl2 / f, maximized at the corners of each cell's (l1, l2) rectangle) for later steps and for the first step over every physical end point pv_end in (0, min(pv_mid, (1 - 1e-4) p1v1)] (the first-step ratio is endpoint-dependent, c3 != 0); 5 x 5 interior samples per cell must stay below each cell's analytic bound, else the build fails closed, as it does when the bound exceeds 1e-3; scope: states of a particle that has lost less than 1e-4 of its initial pv (TABLE_MIN_LOSS) are outside the relative contract, because the memory term drives the f_dM fit to its zero there (e.g. He at 50 MeV/u in aluminium, copper, lead, compact bone and the user material, where the fit is clamped to 0 for both paths); measured analytic bound / sampled maximum, upper and lower domain end: water He 1.10e-5 / 4.1e-6 and 1.36e-5 / 3.2e-6, C 4.38e-5 / 1.97e-5 and 1.47e-5 / 7.1e-6, O 4.37e-5 / 1.98e-5 and 1.94e-5 / 9.3e-6; skeletal muscle He 1.09e-5 / 4.0e-6 and 1.32e-5 / 3.1e-6, C 4.90e-5 / 2.20e-5 and 1.46e-5 / 7.0e-6, O 4.91e-5 / 2.20e-5 and 1.94e-5 / 9.3e-6; largest analytic bound over the eleven built-in materials and a user material (70 tables) He 1.36e-5 (water, lower), C 6.37e-5, O 7.20e-5 (beryllium, upper), at most 7.2 % of the 1e-3 threshold; oxygen at 100 MeV/u in copper and lead cannot be tabulated (a node has no root; fails closed before certification, root option remains); build plus certification 0.5-1.6 s per (material, ion) table. (2026-11-08, Codex REVIEW-ee53867c) The table is demoted to an accelerator and the two preceding certification sentences are superseded: mapping evaluation = exact root solve on every step in both modes (the same Illinois iteration, tolerance |residual| <= 1e-13 or step/bracket <= 1e-14 in ln pv_eq, float64, and status codes as the root path); with mcs_mapping_eval = "table" the bilinear pv_eq/pv value only supplies the starting bracket (one end at the table value, the other max(10 x the table's build-time interpolation estimate, 1e-5) away on the side of the sign change, widened once x10, else the full proton-path bracket, counted in the diagnostic counter mcs_table_bracket_miss, which never invalidates a result); no interpolation bound is claimed (the build-time estimate is a performance diagnostic) and the TABLE_MIN_LOSS scope is withdrawn; the decision-0043 test of the table evaluation is the per-step root guarantee (`test_ion_mcs_table.py::test_accelerated_root_equals_root`: He, C, O at both domain ends in the eleven built-in materials and a user material, >= 2000 random plus edge states each, later and first steps with both end-point extremes: equal status, pv_eq within 1e-8 of the root path (observed <= 1.43e-9), variance within 1e-8 times the f_dM sensitivity (observed <= 5.3e-10), 0 bracket misses in 427 704 evaluations; `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states`: residual <= 1e-10, observed <= 9.4e-14); oxygen at 100 MeV/u in copper and lead still cannot be tabulated (build fails closed, root option); measured MCS cost (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 21.00 us/step (root 40.04), C 150 MeV/u 21.65 us/step (root 39.63), 4.58 / 4.84 ms per history against 8.66 / 8.79; Illinois iterations per step on random table-domain states: accelerated mean 3.8-5.7 (max 26), root 7.8-8.7 (max 33); build 0.11-0.16 s (water) and 0.21-0.31 s (skeletal muscle) per (material, ion) table. (2026-11-09, Codex REVIEW-16c63514, after the sample seal commit 84b8b2cd) Mapping domain rule, decision D1: the exact-Omega mapping is defined where the proton-equivalent state stays inside the proton b-slab range; with multiple scattering on, `validate` root-solves the table nodes of the run's reachable band of initial states in BOTH evaluations (`exact_omega.pv_eq_table_coverage`, cached node solves shared with the table build, so the table rows are bit-identical) and fails closed (`UnsupportedCombinationError` naming the material, the uncovered cells, the lowest pv and largest p1v1/pv and the root status) when any node has no root, so every validated run has table coverage 1.0 (recorded, `mcs_mapping.table.table_coverage`); previously "root" mode validated such a run and every history then failed as `mcs_mapping_failure` (measured: oxygen16 100 MeV/u, 20 histories, python backend: copper 15 / 14 failures and lead 46 / 45 in root / table mode, identical counters), and "table" mode failed closed at the build. A table cell with a failed node is otherwise left invalid and a state there uses the full proton-path bracket (the same exact root; such a step has no starting guess and is NOT counted in `mcs_table_bracket_miss`, whose meaning is now: a tight bracket was tried and had no sign change; the sealed water samples count 0 misses under both readings). With `multiple_scattering=False` no mapping table is built and no coverage check is made (summary `mcs_mapping.evaluation = "root"`, `unused = true`). Measured restriction (`test_ion_mcs_table.py::test_mapping_domain_restrictions_are_measured`, eleven built-in materials x three ions x both domain ends, 66 cases): exactly oxygen16 at the lower domain end (100 MeV/u) in copper (2 of 379 band cells uncovered) and lead (5 of 379), failing nodes at pv 62.6541 MeV with p1v1/pv up to 48.6456, root status 2; a monoenergetic spot scan (not a test) found the failure in copper at 100-175 MeV/u (200 and above validate) and in lead at 100-350 MeV/u (405 validates). The twelve sealed sample configurations (water) are unaffected: table rows and results of full-coverage configurations are bit-identical before and after this change (python backend, 50 histories, seed 11: oxygen16 150 MeV/u dose c0ba46e5…, carbon12 80 MeV/u dose 94b369b9…, table rows e5a79539… / a3575b5d…). Extending the proton path below its current lower edge is a model change outside revision 2 (follow-up). (2026-11-09, Codex REVIEW-f3361441) With `multiple_scattering=False` an ion source now performs NO mapping work at validation: no f_dM domain check (`mcs_mapping.fdm_domain = "not evaluated (multiple_scattering=False)"`, `fdm_extrapolated` false, proton-equivalent energies null), no proton path, no packing, no coverage probe, no table and no root solve; the kernels receive zero placeholder arrays that are read only inside the MCS branch (`ctl.mcs == 1`), verified on the python and warp-cpu backends (oxygen16 100 MeV/u in lead, MCS off: valid result, every counter 0); regression test `test_ion_mcs_table.py::test_mcs_off_runs_no_mapping_machinery` patches every builder (`fdm_domain_status`, `proton_path`, `pack_*`, `pv_eq_table_coverage`, `build_pv_eq_table_and_estimate`, `root_pv_equivalent_batch`, `ExactOmegaMapping.build`) to raise and validates MCS-off ion configurations in both evaluations (the same patches trip with MCS on). MCS-on behaviour, summaries and the full-coverage digests above are unchanged.
