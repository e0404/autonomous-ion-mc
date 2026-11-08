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
| `src/ionmc/physics/exact_omega.py` | `e43eecc309d19d3b0374ed9217a40480e072e284c39794505639a035f321d104` | I0/I11-map notes |
| `src/ionmc/physics/stopping.py` | `89db1de28e5cd1f03efc1622669b58dd766fc88d5d0a91405af190e96094753e` | P9 (iii) |
| `src/ionmc/config.py` | `34be11694b94dc0df3a5c47dba58fe5e8bbe26e13b03e5cc8374beacbef82d77` | I0 note |
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
| `validation/scripts/ion_em/gen_r2_bundles.py` | `fe122218804cbe74f2f6966fbcc67a2d66b1a9e85ec7d2aab57cee7350c8f60e` |
| `tests/ionmc/test_r2_bundles.py` | `54af24fa0baf13d43f4ed25c96989344e09c01354efdda0dc34ba83a5ab60fdb` |
| `reference_requests/operator-directives/2026-11-07-topas-validation-target.md` | `ca079e9f0af990cdff9d8efde60d38ac55e2d0fc0fb9ef22355a81f4a2ab0b26` |
| `reference_requests/interventions/IR-20261102-010435-B59B2D/response.json` | `dc6a610cb67dee50fa0bb96327026eea0b81bb9506e7b898c2e9e44322f69bd0` |
| `validation/plans/v3-008/r2-topas-target.md` | `346ad9518cf14e7090fa5c14548fb8a0a0197489442341629708bc68b048ae11` |


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
| `src/ionmc/physics/ion_mcs.py` | `36f4f3a7246d8e917a1ad74a80d706568c14516b2cd3999d7280b4dc6202ccf0` |
| `src/ionmc/simulation.py` | `2958dbd3069604a25147aafe22a02a988b3ac923156d9e3c5c43ab3fc7b7d814` |
| `src/ionmc/species.py` | `19df88c722ad93ca65e636b870a4b56b1a80df2e215bee39a7623f441daa2322` |
| `src/ionmc/transport/kernels.py` | `db4ea1f096e4919303b5e1d1413f25d78dd7ced149c331983016172b0b099266` |
| `src/ionmc/transport/reference.py` | `f3e3dec9ed5c7bffed827213783af19eee2bbe833f47120550c2c88a0a6b9bb3` |
| `src/ionmc/transport/run.py` | `bbfae2a4c1c2168397ab898204ede3bc2d7c92e1968a08b8156a03326217436b` |
| `src/ionmc/transport/tables.py` | `45542d46fb208ba1ca9ccdcb9e8baea41e47f86a7429b46f8566ad12b1b2619d` |
| `src/ionmc/transport/tally.py` | `b941a46bda49c12e2cda1e27978528b85c08120fc56a7d4c7e172f32953c87d0` |
| `src/ionmc/transport/warp_driver.py` | `cb32245e48d6bffd03092d09fb7421d8b11a8191158a311fd0b9a92e568aeb86` |
| `validation/scripts/ion_em/r2_evaluate.py` | `be1eadcc504968d75fce09b0f7555b80ccb23c2ece18477b08159091f226ec14` |
| `validation/scripts/transport/steps_v8.py` | `e8b73006921b479797e762c1319adf7e1e4edb95594ffac4109827b4ee185a7d` |
| `tests/ionmc/test_ion_mcs_table.py` | `62d495bc4bba5557589e8670cd8f58ced5fc40cf194320fcacfe49b88a4dbb18` |
| `tests/ionmc/test_v3_008_rows_ci.py` | `efce49724d4eef7060bc755e8b4c2ecb60d76141d61993120be0b5577336987d` |
| `tests/ionmc/test_r2_evaluate.py` | `5b996ebb62304e4720bc631ec4bae4b41c7a4bf13ddabdfd1ccc7ea9cf18f222` |
| `tests/ionmc/test_v3_008_steps_v8.py` | `ab4922ea8842751e22c9d1882e9e1dec0ae88a25b55174c99586c4802868b731` |

(c) R2.5 counts, superseded by this note (the R2.5 text is kept as written): repository paths 475 + 15 = **490**; total full sha256 values 519 + 15 = **534**, plus the bundle-set digest (unchanged, `a8a1c27d…`). The count in R2.6 item (3) (475) is superseded accordingly by 490 for the acceptance tests of the seal.

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-7, Codex REVIEW-d3a17ebf: per-table certification of the tabulation bound; attest_verify bootstrap migration-only; manifest sentinel) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` 1b67ae5f… → fbeaaf18…; `src/ionmc/config.py` f8b1c447… → 0b9b34a3…; `validation/plans/v3-008/r2-topas-target.md` 474480e7… → 7aa4d9ee…; `src/ionmc/simulation.py` 783fa561… → cfb8f52d…; `tests/ionmc/test_ion_mcs_table.py` ce87db6f… → 95b02229…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts of the 2026-11-08 post-freeze note (490 paths / 534 values) are unchanged: no path added or removed.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-8, Codex REVIEW-9529fd35: analytic domain-wide certification of the tabulation bound with the relative-contract scope TABLE_MIN_LOSS = 1e-4) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` fbeaaf18… → 53bbc7f9…; `src/ionmc/config.py` 0b9b34a3… → 39bddffb…; `validation/plans/v3-008/r2-topas-target.md` 7aa4d9ee… → a90dc6df…; `src/ionmc/simulation.py` cfb8f52d… → 1d4d701b…; `tests/ionmc/test_ion_mcs_table.py` 95b02229… → 03e1925f…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts (490 paths / 534 values) are unchanged.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-9, Codex REVIEW-ee53867c: the table demoted to an accelerator of the exact root solve on every step; no interpolation-bound claim) — rebound entries (old → new): `tests/ionmc/test_ion_em_mapping_src.py` d81e0915… → 865d3143…; `src/ionmc/physics/exact_omega.py` 53bbc7f9… → 5f4aa165…; `src/ionmc/config.py` 39bddffb… → 34be1169…; `validation/plans/v3-008/r2-topas-target.md` a90dc6df… → 4722d837…; `src/ionmc/physics/ion_mcs.py` 709c61b1… → 36f4f3a7…; `src/ionmc/simulation.py` 1d4d701b… → b855273a…; `src/ionmc/transport/kernels.py` 73bcaf87… → db4ea1f0…; `src/ionmc/transport/reference.py` 103672f5… → f3e3dec9…; `src/ionmc/transport/tally.py` debbe3ed… → 7a746c99…; `tests/ionmc/test_ion_mcs_table.py` 03e1925f… → 62d495bc…; `tests/ionmc/test_v3_008_rows_ci.py` f40b7006… → efce4972…. All other entries unchanged (recomputed against the tree). Freeze commit 0245158b, the R2.1 rule and the counts (490 paths / 534 values) are unchanged.**

**Dated note (2026-11-08, post-freeze implementation commit R2-TOPAS-9b, Codex REVIEW-4bf4f5bc minor: the `mcs_table_bracket_miss` descriptions now state that the counter counts every miss of the initial tight bracket, widened or full; documentation only) — rebound entries (old → new): `src/ionmc/physics/exact_omega.py` 5f4aa165… → e43eecc3…; `src/ionmc/simulation.py` b855273a… → 2958dbd3…; `src/ionmc/transport/tally.py` 7a746c99… → b941a46b…. All other entries unchanged (recomputed against the tree). Counts unchanged (490 paths / 534 values).**

**Dated note (2026-11-08, ionmc sample chunk plan appended to r2-topas-target.md section 8; records only) — rebound entries (old → new): `validation/plans/v3-008/r2-topas-target.md` 4722d837… → 346ad951…. Counts unchanged (490 paths / 534 values).**

**Dated note (2026-11-08, steps_v8 writes the sample npz with np.savez_compressed (lossless; the lateral grids are mostly zero, ~385 MB → tens of MB per chunk; chosen before any sealed sample, the single chunk run so far is repeated) — rebound entries (old → new): `validation/scripts/transport/steps_v8.py` 8d2617cd… → e8b73006…. Counts unchanged (490 paths / 534 values).**
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

**Dated note (2026-11-08, Amendment):** mapping evaluation = tabulation (decision 0043 section 3, "Evaluation", the allowed alternative), chosen before any revision-2 ionmc sample or TOPAS evidence run; reason: measured cost. Table: `pv_eq/pv` per (material, ion) on a uniform 128 x 128 grid in (ln pv, ln(p1v1/pv)) (a shear of the decision's (ln pv, ln p1v1) plane), nodes computed by the unchanged root solve `pv_equivalent` at validation (only the nodes of cells meeting the band of the run's initial states pv(E0 - 6 sigma_E) ... pv(E0 + 6 sigma_E); every other node invalid, a state there fails closed as `mcs_mapping_failure`), bilinear interpolation (`exact_omega.pv_equivalent_table`, shared by the python twin and the Warp kernel); option `PhysicsOptions.mcs_mapping_eval` = "table" (default for z >= 2) or "root" (the per-step root solve, retained); the summary records the evaluation, table sizes, axes, domain and the sha256 of the table rows. Bound test `tests/ionmc/test_ion_mcs_table.py::test_table_variance_error_bound` (He, C, O in water and skeletal muscle; >= 2000 random states plus edge states per ion and material, upper and lower domain end, thin-slab edge included; bound |var_table/var_root - 1| <= 1e-3): measured maximum 2.05e-5 (C in skeletal muscle; He 4.4e-6 / 5.5e-6, C 1.83e-5 / 2.05e-5, O 1.84e-5 / 2.04e-5 in water / muscle); transport-level counterpart `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states` max 6.8e-6. Measured MCS cost per step (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 0.16 us (root 39.85 us), C 150 MeV/u 0.17 us (root 39.57 us). The model (decision 0043 section 3) is unchanged; the root solve remains available as an option. The rebinding of the changed bound files and the binding of the new implementation paths are recorded in the §2 dated note "2026-11-08, post-freeze implementation commits". (2026-11-08, Codex REVIEW-d3a17ebf) The bound is in addition certified deterministically at build for every produced table (`exact_omega.certify_pv_eq_table`: the root solve at 5 interior points per solved cell, the centre and +-1/4 of the cell in both axes; per-step variance of the first and later steps; a table above TABLE_VARIANCE_BOUND = 1e-3 fails closed; the maximum, cells and points are recorded in the summary), measured maxima at the upper domain end (128 x 128, 628 cells, 3140 points per table) He / C / O: water 4.11e-6 / 1.97e-5 / 1.98e-5, skeletal muscle 4.02e-6 / 2.20e-5 / 2.20e-5, over all eleven built-in materials and a user-defined material at most 4.44e-6 / 2.87e-5 / 3.21e-5 (`test_ion_mcs_table.py::test_certification_every_builtin_and_user_material`). (2026-11-08, Codex REVIEW-9529fd35) The five-sample certification is replaced by a domain-wide analytic bound per solved cell (`exact_omega.certify_pv_eq_table`): the tensor-product bilinear error bound hx^2/8 sup|q_xx| + hy^2/8 sup|q_yy| of q = pv_eq/pv, second derivatives from second differences of the root solve on a 2x-refined grid (3 x 3 stencil per cell, maximum over the cell and its solved neighbours, which covers cells where a second derivative changes sign), times safety factor 2 (q is smooth in ln pv and ln r, P11; the factor covers the variation of the second derivative within a cell; empirically, over the 70 tables below the sampled cross-check passes everywhere with factor 1 and fails in every table with 0.75, so factor 2 leaves a per-cell margin of at least 2), converted exactly into the relative per-step variance error (only pv_eq differs between the paths and f_dM is linear in log10 pv_eq at fixed memory term l1: var_table/var_root - 1 = (c2 + c3 l1) dl2 / f, maximized at the corners of each cell's (l1, l2) rectangle) for later steps and for the first step over every physical end point pv_end in (0, min(pv_mid, (1 - 1e-4) p1v1)] (the first-step ratio is endpoint-dependent, c3 != 0); 5 x 5 interior samples per cell must stay below each cell's analytic bound, else the build fails closed, as it does when the bound exceeds 1e-3; scope: states of a particle that has lost less than 1e-4 of its initial pv (TABLE_MIN_LOSS) are outside the relative contract, because the memory term drives the f_dM fit to its zero there (e.g. He at 50 MeV/u in aluminium, copper, lead, compact bone and the user material, where the fit is clamped to 0 for both paths); measured analytic bound / sampled maximum, upper and lower domain end: water He 1.10e-5 / 4.1e-6 and 1.36e-5 / 3.2e-6, C 4.38e-5 / 1.97e-5 and 1.47e-5 / 7.1e-6, O 4.37e-5 / 1.98e-5 and 1.94e-5 / 9.3e-6; skeletal muscle He 1.09e-5 / 4.0e-6 and 1.32e-5 / 3.1e-6, C 4.90e-5 / 2.20e-5 and 1.46e-5 / 7.0e-6, O 4.91e-5 / 2.20e-5 and 1.94e-5 / 9.3e-6; largest analytic bound over the eleven built-in materials and a user material (70 tables) He 1.36e-5 (water, lower), C 6.37e-5, O 7.20e-5 (beryllium, upper), at most 7.2 % of the 1e-3 threshold; oxygen at 100 MeV/u in copper and lead cannot be tabulated (a node has no root; fails closed before certification, root option remains); build plus certification 0.5-1.6 s per (material, ion) table. (2026-11-08, Codex REVIEW-ee53867c) The table is demoted to an accelerator and the two preceding certification sentences are superseded: mapping evaluation = exact root solve on every step in both modes (the same Illinois iteration, tolerance |residual| <= 1e-13 or step/bracket <= 1e-14 in ln pv_eq, float64, and status codes as the root path); with mcs_mapping_eval = "table" the bilinear pv_eq/pv value only supplies the starting bracket (one end at the table value, the other max(10 x the table's build-time interpolation estimate, 1e-5) away on the side of the sign change, widened once x10, else the full proton-path bracket, counted in the diagnostic counter mcs_table_bracket_miss, which never invalidates a result); no interpolation bound is claimed (the build-time estimate is a performance diagnostic) and the TABLE_MIN_LOSS scope is withdrawn; the decision-0043 test of the table evaluation is the per-step root guarantee (`test_ion_mcs_table.py::test_accelerated_root_equals_root`: He, C, O at both domain ends in the eleven built-in materials and a user material, >= 2000 random plus edge states each, later and first steps with both end-point extremes: equal status, pv_eq within 1e-8 of the root path (observed <= 1.43e-9), variance within 1e-8 times the f_dM sensitivity (observed <= 5.3e-10), 0 bracket misses in 427 704 evaluations; `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states`: residual <= 1e-10, observed <= 9.4e-14); oxygen at 100 MeV/u in copper and lead still cannot be tabulated (build fails closed, root option); measured MCS cost (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 21.00 us/step (root 40.04), C 150 MeV/u 21.65 us/step (root 39.63), 4.58 / 4.84 ms per history against 8.66 / 8.79; Illinois iterations per step on random table-domain states: accelerated mean 3.8-5.7 (max 26), root 7.8-8.7 (max 33); build 0.11-0.16 s (water) and 0.21-0.31 s (skeletal muscle) per (material, ion) table.
