"""Write a Markdown summary of an archived run_all.sh results directory.

Every value is copied or computed from the raw step files; nothing is typed by
hand. Steps with a non-zero exit code are reported as failed/timed out and
excluded from statistics. For repeated steps (``NN-rK-name``) the median and
range of the metrics found in the output are computed.
"""

import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

NOISE = re.compile(r"^(Warp CUDA (error|warning)|Warp DeprecationWarning|#|$)")
METRICS = {
    "hist/s": re.compile(r"hist/s=([0-9.]+)"),
    "cold compile s": re.compile(r"cold compile ([0-9.]+)s"),
    "philox ns/uniform": re.compile(r"philox: ([0-9.]+) ns/uniform"),
    "warp-pcg ns/uniform": re.compile(r"warp-pcg: ([0-9.]+) ns/uniform"),
}


def parse(path: Path) -> tuple[str, str, list[str]]:
    lines = path.read_text().splitlines()
    command = next((l[len("# command: ") :] for l in lines if l.startswith("# command: ")), "?")
    exits = [l[len("# exit=") :] for l in lines if l.startswith("# exit=")]
    code = exits[-1] if exits else "missing"
    return command, code, [l for l in lines if not NOISE.match(l)]


def main(directory: str) -> None:
    root = Path(directory)
    env = (root / "environment.txt").read_text()
    sha = next(l.split("=", 1)[1] for l in env.splitlines() if l.startswith("git_sha="))
    expected = set((root / "manifest.txt").read_text().split())
    present = {p.stem for p in root.glob("[0-9][0-9]-*.txt")}
    if present != expected:
        missing, extra = sorted(expected - present), sorted(present - expected)
        raise SystemExit(f"archive does not match manifest.txt: missing={missing} extra={extra}")
    for path in root.glob("[0-9][0-9]-*.txt"):
        header = path.read_text().splitlines()[:4]
        if f"# git_sha: {sha}" not in header:
            raise SystemExit(f"{path.name}: step provenance does not match environment.txt (SHA {sha})")
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
