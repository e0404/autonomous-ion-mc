"""Write a Markdown summary of an archived run_all.sh results directory.

Every value is copied or computed from the raw step files; nothing is typed by
hand. The expected set of step files is derived here from the fixed suite
definition and the validated ``repeats`` value recorded in ``environment.txt``,
never from the archive's own manifest; the manifest and the files must both
match that set exactly. Steps with a non-zero exit code are reported as
failed/timed out and excluded from statistics. For repeated steps
(``NN-rK-name``) the median and range of the metrics found in the output are
computed.
"""

import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

# Fixed suite definition (must match run_all.sh; a regression test checks that).
REPEATED_STEPS = [
    "01-compile-cold-backward-on",
    "02-compile-cold-backward-off",
    "03-toy-f32-1thread",
    "04-toy-f64-1thread",
    "05-toy-f32-4threads",
    "06-toy-f32-8threads",
    "07-toy-f32-16threads",
    "13-philox-kat-and-cost",
]
SINGLE_STEPS = [
    "08-precision",
    "09-rng-seed-dupes",
    "10-rng-overlap-1e5x1000",
    "11-rng-overlap-1e6x1000",
    "12-rng-overlap-1e6x2000",
]
NOISE = re.compile(r"^(Warp CUDA (error|warning)|Warp DeprecationWarning|#|$)")
METRICS = {
    "hist/s": re.compile(r"hist/s=([0-9.]+)"),
    "cold compile s": re.compile(r"cold compile ([0-9.]+)s"),
    "philox ns/uniform": re.compile(r"philox: ([0-9.]+) ns/uniform"),
    "warp-pcg ns/uniform": re.compile(r"warp-pcg: ([0-9.]+) ns/uniform"),
}


def expected_names(repeats: int) -> set[str]:
    names = set(SINGLE_STEPS)
    for step in REPEATED_STEPS:
        number, rest = step.split("-", 1)
        names.update(f"{number}-r{k}-{rest}" for k in range(1, repeats + 1))
    return names


def read_environment(root: Path) -> tuple[str, int]:
    values = dict(line.split("=", 1) for line in (root / "environment.txt").read_text().splitlines() if "=" in line)
    sha = values.get("git_sha", "")
    if not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise SystemExit("environment.txt: missing or malformed git_sha")
    repeats = values.get("repeats", "")
    if not re.fullmatch(r"[1-9][0-9]*", repeats):
        raise SystemExit("environment.txt: missing or non-positive repeats")
    return sha, int(repeats)


def verify_archive(root: Path) -> str:
    sha, repeats = read_environment(root)
    expected = expected_names(repeats)
    manifest_lines = [l for l in (root / "manifest.txt").read_text().splitlines() if l.strip()]
    if not manifest_lines or len(set(manifest_lines)) != len(manifest_lines):
        raise SystemExit("manifest.txt is empty or contains duplicates")
    if set(manifest_lines) != expected:
        raise SystemExit(
            f"manifest.txt does not match the suite for repeats={repeats}: "
            f"missing={sorted(expected - set(manifest_lines))} extra={sorted(set(manifest_lines) - expected)}"
        )
    present = {p.stem for p in root.glob("[0-9][0-9]-*.txt")}
    if present != expected:
        raise SystemExit(
            f"archive files do not match the suite for repeats={repeats}: "
            f"missing={sorted(expected - present)} extra={sorted(present - expected)}"
        )
    for path in root.glob("[0-9][0-9]-*.txt"):
        header = path.read_text().splitlines()[:4]
        if f"# git_sha: {sha}" not in header:
            raise SystemExit(f"{path.name}: step provenance does not match environment.txt (SHA {sha})")
    return sha


def parse(path: Path) -> tuple[str, str, list[str]]:
    lines = path.read_text().splitlines()
    command = next((l[len("# command: ") :] for l in lines if l.startswith("# command: ")), "?")
    exits = [l[len("# exit=") :] for l in lines if l.startswith("# exit=")]
    code = exits[-1] if exits else "missing"
    return command, code, [l for l in lines if not NOISE.match(l)]


def main(directory: str) -> None:
    root = Path(directory)
    verify_archive(root)
    print(f"# Archived measurement run `{root.name}`\n")
    print("## Environment\n\n```")
    print((root / "environment.txt").read_text().rstrip())
    print("```\n\n## Steps\n")
    print("| Step | Command | Exit | Result lines (verbatim, warnings removed) |")
    print("|---|---|---|---|")
    groups: dict[str, list[tuple[str, list[str]]]] = defaultdict(list)
    for path in sorted(root.glob("[0-9][0-9]-*.txt")):
        command, code, result = parse(path)
        status = code if code == "0" else f"**{code} (failed or timed out)**"
        cell = "<br>".join(l.replace("|", "\\|") for l in result) or "(no output)"
        print(f"| {path.stem} | `{command}` | {status} | {cell} |")
        m = re.match(r"(\d\d)-r\d+-(.*)", path.stem)
        if m and code == "0":
            groups[f"{m.group(1)}-{m.group(2)}"].append((command, result))
    print("\n## Statistics of repeated steps (successful repeats only)\n")
    print("| Step | Metric | Repeats | Median | Min | Max |")
    print("|---|---|---|---|---|---|")
    for group, runs in sorted(groups.items()):
        for metric, pattern in METRICS.items():
            values = [float(pattern.search(l).group(1)) for _, result in runs for l in result if pattern.search(l)]
            if values:
                print(
                    f"| {group} | {metric} | {len(values)} | {statistics.median(values):g} | "
                    f"{min(values):g} | {max(values):g} |"
                )


if __name__ == "__main__":
    main(sys.argv[1])
