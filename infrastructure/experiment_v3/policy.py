"""Defense-in-depth guard against direct GitHub CLI CI/integration access."""

import json
import re
import sys


def direct_github_command(command):
    # Covers ordinary qualified paths, env wrappers and shell sequences. This is
    # an accidental-use guard, NOT a security boundary against arbitrary Python.
    return bool(re.search(r'(?<![\w-])(?:[^\s;|&"\x27]*/)?gh(?=[\s"\x27]|$)', command))


def main():
    payload = json.load(sys.stdin)
    if payload.get("tool_name") == "Bash" and direct_github_command(
        payload.get("tool_input", {}).get("command", "")
    ):
        print(
            "Use controlled CI/integration/release MCP tools; "
            "direct gh CLI is disabled.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
