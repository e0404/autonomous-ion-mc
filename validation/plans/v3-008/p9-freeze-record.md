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

**Dated note (2026-11-01, orchestrator) — FC-1 DISPOSED.** Record `validation/plans/v3-008/p1c-completion.md` (sha256 e3cd32585140fc36a45c7059ec21fa3e782144f7ad5557bc953856c1f1356e53; computed in the scratchpad from the committed reference `moliere_hanson.py` 5a99b3b8… unchanged, variants as runtime monkeypatches only). Outcome: **Pb CLEARED by attribution** — the P1b Pb failure (+1.10/+1.32 % at x/R1 ≥ 0.9) is reproduced by the range-table input (the ionmc Bethe Pb table is +0.85 % longer than Gottschalk's ICRU 49 R1 and +2.6–3.1 % in residual range at 20–41 MeV, amplified by d ln θ/d ln x = 2.3–5.7 near the end of range); on an ICRU 49 node path the committed reference reproduces Gottschalk's Pb θ_Hanson to ≤ 0.044 % (≤ 0.024 % on a PSTAR path), controls Be/Al/Cu ≤ 0.105 %; the Fano form (Δb 0.05–0.06 for Pb), the compound rule (n/a) and the quadrature (≤ 1.2e-6) do not explain the failure. The committed test (ionmc Bethe ranges) still fails at x/R1 ≥ 0.9, so its xfail(strict=True) stays as committed (no bound file changes before the freeze); P1 is CLOSED with the attribution and the proton validation statement reads 'passes Be/Al/Cu; Pb cleared by attribution to the range-table input'. Water/compound check: Gottschalk 1993 (NIM B 74, 467) is paywalled and lists no water target; arXiv:1610.01279 tabulates measured θ0 and θ_Hanson for polystyrene, Lexan, Lucite, Nylon, Teflon and C (no uncertainties) — the committed reference agrees with Gottschalk's θ_Hanson (LOOKUP) within −0.57 to −0.04 % (polystyrene −0.16/−0.20 %, Lexan C16H14O3 −0.15/−0.10 %), which closes OPEN item 7 for H/C/N/O/F compounds and discriminates against the electron- and mass-weighted Fano averages (+4.1 % / +1.8 %); **water remains not validated directly** (the only water θ_Hanson value available has two digits; used as a check on b: ours low by 0.040 ± 0.002, which moves θ by ≤ 0.16 % at I11-like B). Side finding (outside the V3-008 domain, recorded for high-Z use): ionmc Bethe Pb ranges are 0.8–3 % longer than ICRU 49 at 20–160 MeV. TO-VERIFY: Gottschalk 1993's measurement uncertainties; the exact LOOKUP compositions. The P1 proton tests skip unless IONMC_CACHE_DIR names the objects store (test default cache path; recorded, not changed).

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
| `tests/ionmc/test_ion_em_moliere_hanson.py` | `789fd4b1714f18a53fdb3db840228209c5b09a910c1651b1d67b646d93b87e0c` | P1 |
| `tests/ionmc/test_ion_em_exact_omega.py` | `12cbdba28530db9561ac469acba5d45f7533690f5587dd74d44b2143b542baee` | P11 (I11-map/I11-rel, reference script) |
| `tests/ionmc/test_ion_em_lindhard_sorensen.py` | `d35ca4e34610e153a6aa31f672a180cab4092572c8c839dc4087e606024dc6ef` | P3 |
| `tests/ionmc/test_ion_em_ss_budget.py` | `ca49735b5b00831aeaabd59d63faecfc23b98fe0d4247024ca963d7e9a77f097` | P9 (iv) |
| `tests/ionmc/test_ion_em_ls_option.py` | `c67966c8d2c3d9564032e33e666321dd11bd8d2386780afa0cb74850455ca84d` | P9 (iii) |
| `tests/ionmc/test_ion_em_mapping_src.py` | `e11b7b655044469efadb1d42d3d7a94d5b51c26059aff3da7bdf6f735294c93b` | I11-map (src), I0, I1 |
| `tests/ionmc/test_v3_008_seed_families.py` | `5a11ca5128709fe23618a8179dab0288c2d01fbe0431f5b328ca45d5107e9114` | P7, §C |

### 2.3 Implementation already on the branch (src)

| File | sha256 | Bound by |
|---|---|---|
| `src/ionmc/physics/lindhard_sorensen.py` | `800d1f0a13524af5bad3946032f894479f606ba9b1482bf3ea190eeb9b4d3a6a` | P9 (iii) |
| `src/ionmc/physics/exact_omega.py` | `e8427d23f71e7b524b434237d3343a7783ff375822da71ff38a67fac1c9d541e` | I0/I11-map notes |
| `src/ionmc/physics/stopping.py` | `951f29acaae262255e1f0b6bc831d11118fe55025bf0a8cfb93c4cefad8db828` | P9 (iii) |
| `src/ionmc/config.py` | `2ead8d09a31dd3eb16ed97dba12f90d8a4e937efc56eabc0a8d718f2f5ee3f90` | I0 note |
| `pyproject.toml` | `99ff83af258ccead7d0e81516abdf020730720aa3596a64d5292493b250dfc9f` | lint exceptions for the two mapping files (added by the planner; the only other changed file outside the plan set) |

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
| `validation/plans/v3-008/p10-single-scattering-availability.md` | `bea86277604b20c2d396c12d0c562d5f9e21acd0fa94ba0c560cf5ee0ce02881` |
| `validation/plans/v3-008/p11-exact-omega-mapping.md` | `6116648ae7ef2de61aca353616213e7709f9eeb012af3ab8aaf3d3ae55d81a07` |
| `validation/plans/v3-008/research/research-models.md` | `d194a516967f290e9aa65f693c3c414282c7d45cde5c1226ab08464a46b396cd` |
| `validation/plans/v3-008/research/research-data-sources.md` | `2bf40c44d73dd1d5ce6afe6d05d7891128cb3b4468af2afbf196f4b99afa0732` |

**Dated note 2026-11-01 (FC-2 record, bound by the §0 rule):**

| File | sha256 |
|---|---|
| `validation/plans/v3-008/p2b-follow-ups.md` | `9c70eb1b9d49f0375af7976232817dbda925576eeb4b64e31aee61b33dfb5457` |

Not bound by hash, because the freeze commit itself changes them: the plan, decision 0043, `review.md` and this
file. They are bound by the freeze commit SHA.

**Dated note (2026-11-01, FC-1):** `validation/plans/v3-008/p1c-completion.md` `e3cd32585140fc36a45c7059ec21fa3e782144f7ad5557bc953856c1f1356e53` (added by the freeze commit; §7 counts: +1 repository file).

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

`validation/reference_cases/README.md`: `2c160fb585dc571d26753aa4a1a43e44279ae4f41dfbeba480f96452a07b491e`.

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
