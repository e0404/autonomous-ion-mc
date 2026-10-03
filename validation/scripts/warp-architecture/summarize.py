"""Write a Markdown summary of an archived run_all.sh results directory.

Every value is copied from the raw step files; nothing is typed by hand.
A step whose exit code is not 0 is reported as failed/timed out.
"""

import re
import sys
from pathlib import Path

NOISE = re.compile(r"^(Warp CUDA (error|warning)|Warp DeprecationWarning|#|$)")


def main(directory: str) -> None:
    root = Path(directory)
    print(f"# Archived measurement run `{root.name}`\n")
    print("## Environment\n\n```")
    print((root / "environment.txt").read_text().rstrip())
    print("```\n\n## Steps\n")
    print("| Step | Command | Exit | Result lines (verbatim, warnings removed) |")
    print("|---|---|---|---|")
    for path in sorted(root.glob("[0-9][0-9]-*.txt")):
        lines = path.read_text().splitlines()
        command = next((l[len("# command: ") :] for l in lines if l.startswith("# command: ")), "?")
        exits = [l[len("# exit=") :] for l in lines if l.startswith("# exit=")]
        code = exits[-1] if exits else "missing"
        result = [l for l in lines if not NOISE.match(l)]
        status = code if code == "0" else f"**{code} (failed or timed out)**"
        cell = "<br>".join(l.replace("|", "\\|") for l in result) or "(no output)"
        print(f"| {path.stem} | `{command}` | {status} | {cell} |")


if __name__ == "__main__":
    main(sys.argv[1])
