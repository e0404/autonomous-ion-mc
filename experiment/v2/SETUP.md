# Preparing and running experiment v2

This is infrastructure setup, not IonMC v2 scientific implementation. Historical
root protocol/requirements and kickoff-v1 remain unchanged; v1 is reconstructable
from `experiment/v1/provenance.json`. The setup branch keeps the completed v1
package. The active v2 checkout is now an independent clone with restricted
v2-only history access; see [the condition 2.1 amendment](HISTORY-ISOLATION.md).
The old fresh-start generator below records how the initial snapshot was made;
it does not produce the final isolated runtime checkout.

## What changes from v1

| V1 observation | V2 condition |
| --- | --- |
| Empty release selection could claim readiness | Nonempty frozen membership, exact-SHA evidence, complete MUST coverage and full-scope gate |
| Planning flags silently omitted physics | Public capability contract; requested/effective configuration and discriminating tests |
| Thin steps, grid shifts and source translations exposed defects | Required adversarial numerical falsification; explicit domains and statistical scaling |
| Related/proton-derived evidence was used for ion claims | Evidence taxonomy, lineage/use roles and independently justified ion-specific evaluation |
| References depended on ad hoc external execution | Autonomous native-input TOPAS, MCsquare and FRED controlled services |
| Fragmentation, biology, provenance and user workflow incomplete | Central heavy-ion scope, external lookup support, persisted metadata and clean-install acceptance |
| Performance emphasized stripped kernels | Fixed planning workload envelopes, hardware/resource metrics and pre-optimization target freeze |
| No recorded interventions; approval prompts could stall | Operational evidence-blocker rule, distinct notification hooks and delivery preflight |
| Direct CI CLI / outside-path typo / early milestone stops | Restricted CLI, portable checkout paths, controlled CI, continuous-execution prompt |

These are new experimental boundaries, not retroactive claims that v1 met them.
The audit report's SHA256 and archive SHA256 are recorded, without copying its
large raw reference corpus into a prescribed v2 validation matrix. Its exact
numerical thresholds/cases are proposals, not imposed scientific answers.

## Reference installation and smoke execution

Run these maintainer commands from the setup checkout (Python 3.12+, bubblewrap).
No engine downloads, builds, license acceptance or runtime tests execute on import.
Registration snapshots the hashes of every runtime file/library/data file and
in-tree symlink. A later change fails the runner until explicitly reprovisioned.
Provisioning is maintainer-only and is not exposed to the autonomous agent.

```bash
python3 -m infrastructure.experiment_v2.provision topas \
  --root /home/wahln/topas/install --geant4-data /home/wahln/geant4/G4DATA \
  --version 'OpenTOPAS 4.3 / Geant4 11.4.2' \
  --source https://github.com/OpenTOPAS/OpenTOPAS
python3 -m infrastructure.experiment_v2.provision mcsquare \
  --root /path/to/dedicated/MCsquare --version '<exact source revision/build>' \
  --source https://gitlab.com/openmcsquare/MCsquare
python3 -m infrastructure.experiment_v2.provision fred \
  --root /path/to/licensed/fred_3.76.0_Linux_gcc_14.2.1 \
  --version 'FRED 3.76.0 gcc14.2.1' --source https://www.fred-mc.org
```

TOPAS uses `<runtime>/bin/topas input.txt` and its registered Geant4 data. MCsquare
uses an explicitly selected `MCsquare_linux config.txt`, with its Materials/BDL/
scanner data hashed; no CPU-dispatch launcher can silently choose another binary.
The setup workstation can use the public static binary from the developer mirror
at `211eefe6eaf2b8572d196d17f546f35ffb0ae0cf`. Its native banner and binary hash,
not the repository date alone, identify the executable. Make it executable before
registration. Upstream source at `85bf2911ddadbb40bd41b1a18b4b247536bb33ab`
uses Intel classic compiler conventions; this setup does not claim a rebuilt
modern-compiler MCsquare. FRED uses its licensed native Python launcher (`-f`);
it may invoke engine subprocesses/preparsers inside the sandbox. The host service
itself never shell-interpolates argv. FRED defaults here to CPU with explicit
resource configuration; OpenCL/GPU provisioning is separate and unqualified.

Public interface documentation:
[OpenTOPAS installation](https://opentopas.readthedocs.io/en/latest/getting-started/install.html),
[MCsquare source and README](https://gitlab.com/openmcsquare/MCsquare),
[FRED command options](https://www.fred-mc.org/Manual_3.76/Control/Command%20Line%20Options.html).
These interfaces leave native inputs under autonomous scientific control.

For a clean committed checkout containing the smoke bundles:

```bash
python3 -m infrastructure.experiment_v2.reference --engine topas \
  --case experiment/v2/references/topas
python3 -m infrastructure.experiment_v2.reference --engine mcsquare \
  --case experiment/v2/references/mcsquare
python3 -m infrastructure.experiment_v2.reference --engine fred \
  --case experiment/v2/references/fred
```

Cases include metadata and expected native outputs. The tiny synthetic CT and
single-spot plan are execution fixtures, not a validation corpus or commissioned
beam model. Runtime fields are **declared**, not independently parsed/verified
physics. Empty/missing outputs fail even if the native process returns zero.
Inspect native banners, actual histories, settings, units and output content
before scientific comparison. Smoke success is not physics qualification.

Default state is `~/.local/share/ionmc-experiment/experiment-v2`, with `engines.json`, raw
reference attempts and content-addressed data; it is inaccessible for agent-shell
writes in the v2 settings. MCP exposes bounded raw artifact reads and copies
verified datasets into ignored task caches. A configuration ID hashes exact code,
inputs, registered runtime and execution settings; each retry has its own attempt
ID. Full stdout/stderr and raw outputs remain archived. Limits: 2 hours wall time,
512 MiB per output/log file, 16 kB per inline stream. There is no aggregate disk
quota; the operator should allocate/monitor storage. System libraries come from
read-only host OS mounts, so reproduction requires the recorded host environment
in addition to registered runtime hashes. This is not a portable container image.
The native code has no network, host home or credentials; only runtime trees and
the disposable run directory are mounted. WSL GPU access is opt-in. Execution
limits and isolation are infrastructure controls, not a hostile-code proof.

## Data and release evidence interfaces

`v2-research` supplies engine discovery/run/artifact reads, public HTTPS GET or
form POST acquisition, role assignment, verified cache materialization, and
scientific event recording. No complete dataset registry is preselected. Downloads
have source/license/citation/request/hash metadata; transforms should record parent
hashes and exact code in the agent's evidence records. Private addresses, embedded
credentials, credential-bearing proxies and non-HTTPS destinations are rejected. Authentication-free
service proxy settings are used when required for public egress. Public
DNS is checked on each redirect; this application check is not a DNS-rebinding
firewall. Use trusted public scientific sources.

The agent creates `validation/release-plan.json` and commits it. The schema is
exercised by `tests/infrastructure/v2/test_release.py`: each suite has id, category,
requirements, backends, acceptance/rationale and, for physics, species and required
independence. Categories include physics/numerical/capability/provenance/performance/
workflow. `performance_targets` maps each workload to numeric metric bounds
(`direction: min|max`, `value`), with hardware and calibration artifacts.
No passing example is supplied as real scientific evidence.

Reports contain code_sha, dirty=false, scope=full, filters=[], plan_sha256,
requirements (every ID satisfied), suite results with nonzero counts/no skips,
backends, hashed artifact paths, independent evidence lineage/review and measured
performance against each target. The checker verifies artifact bytes; scientific
truth and sufficiency still require review. Reports/artifacts belong outside Git
(or ignored output paths) so exact-SHA evidence does not require changing the SHA.

```bash
python3 -m infrastructure.experiment_v2.release --freeze-sha <plan-commit> \
  --report /path/to/report.json --artifacts /path/to/artifacts \
  --output /path/to/qualification.json
```

The `v2-release` MCP re-evaluates evidence and performs the protected release PR,
exact-head CI check, merge, immutable tag and GitHub release. Its state-relative
report paths use the protected v2 state directory. Stage release evidence there
through the evidence-staging tool. Promotion is retriable and never accepts a
standalone readiness boolean. Numeric targets are intentionally unratified until
the bounded calibration phase; this blocks scientific release, not infrastructure
preparation.

## Notifications and runtime policy

The existing `IONMC_NTFY_URL` is supplied privately by the operator's launcher;
never commit or print it. Both permission hooks forward generic event IDs using
the intervention notification transport, record delivery status, and return no
approval. Actual approval remains explicit in Claude's UI. Structured scientific
interventions retain their separate request/response records. Hook details are in
the [official Claude hook reference](https://code.claude.com/docs/en/hooks).
Delivery and hook availability must be checked on the installed Claude version.
No endpoint in the launch environment means no claimed out-of-band delivery.

Direct `gh` permission is removed and a PreToolUse guard rejects ordinary direct
or wrapped gh commands. This prevents accidental policy bypass, not a determined
agent reconstructing GitHub requests in arbitrary Python. Credentials remain
outside its shell sandbox; trusted integration services have the required access.
No special rule is added for the historical `/home/wahlm` typo. The user-level
blockReadsOutsideWorkingDirectories observation is retained as an incident;
v2 uses portable declared worktree paths and does not modify user configuration.

## Original fresh-start preparation (operator provenance)

From a clean, committed setup checkout:

```bash
python3 -m infrastructure.experiment_v2.bootstrap \
  --destination /home/wahln/aiprojects/ion-mc-v2
```

This creates `v2/develop`, `v2/main` and annotated `experiment-v2-start` without
rewriting v1. It fails rather than overwriting existing branches/worktrees/tags.
The committed condition manifest identifies the setup SHA and baseline. Maintainer
publication must push these new refs and set the same protective requirements as
v1 (PRs, no force/deletion, required pre-commit/tests/docs, linear v2/develop).
Do not start unattended until protection, exact-head lightweight CI, reference
smokes and notification delivery pass. Task pushes name origin and branch explicitly; they do not write Git tracking
configuration (the setup environment has a pre-existing repository config lock).
The current setup test record distinguishes
local checks from remote CI and delivery checks.

For an independent replacement checkout, follow HISTORY-ISOLATION.md rather than
using the shared-worktree bootstrap. With the operator notification variable
supplied by the existing launcher, start a new session from the v2 clone:

```bash
cd /home/wahln/aiprojects/ion-mc-v2
~/bin/claude-ionmc --preflight
~/bin/claude-ionmc "Read experiment/prompts/kickoff-v2-isolated.md and carry out that autonomous run."
```

The bootstrap and ready conditions are preserved separately as
`experiment-v2-start` and `experiment-v2-ready`; use current `v2/develop` and the
ready tag for the original setup state. Condition 2.1 is a subsequent versioned
amendment on v2/develop; its clone intentionally has no inherited tags.
No historical tag is moved.

Run the Claude command only after preflight succeeds. Calibration remains the
first scientific-run phase. This preparation task never launches it.

## Launcher and experiment storage

The operator launcher `~/bin/claude-ionmc` sets IONMC_EXPERIMENT_ID first.
Active caches are under `~/.cache/ionmc-experiment/<experiment-id>` and state under
`~/.local/share/ionmc-experiment/<experiment-id>`. For this condition the ID is
`experiment-v2`. The repository's Claude settings use the same literal ID;
changing an experimental condition requires updating those settings too.

| Purpose | Location relative to the experiment root |
| --- | --- |
| General development caches | cache root (XDG_CACHE_HOME) |
| Host environment and CUDA/Warp caches | cache root / host-runner |
| Engine registry, reference runs/events, datasets and runtimes | state root |
| Hook and orchestration telemetry | state root / telemetry |
| Codex raw output | state root / raw/codex |
| Exact-SHA validation records | state root / validation |
| Host execution records | state root / host-runs |
| Launcher temporary files | state root / tmp |

IONMC_V2_STATE remains an optional state override. Existing logging variable names
remain supported. Registered runtime data and historical execution evidence are
preserved during migration; do not clear state directories to start a run.
The host virtual environment is recreated at its new path with the same installed
package versions because moving a venv can leave invalid absolute shebangs.
The old host cache remains archived. Legacy state paths have compatibility
symlinks for earlier run records. New execution writes to the experiment-ID roots.

The launcher must enter `ion-mc-v2`, invoke the Python module
`infrastructure.experiment_v2.preflight` (the module name is not the experiment
ID), and end with exactly one `exec claude "$@"`. Its `--preflight` mode runs
smokes and notification delivery using the same exported environment as Claude.
Keep the operator's notification URL private. The launcher itself is maintained
outside Git; the repository contains no copy of its secret-bearing contents.
