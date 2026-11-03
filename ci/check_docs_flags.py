#!/usr/bin/env python3
"""Validate that every CLI flag used in the docs actually exists.

The user-facing documentation is the first thing a new user runs, so a flag that
does not exist (``-b 10`` where only ``--batch-size`` is declared) or a flag
spelled with an underscore (``--custom_db`` vs ``--custom-db``) is a hard stop
for them. This script parses the click option declarations out of
``src/baktfold/__init__.py`` and checks every ``baktfold ...`` command line in
``docs/`` and ``README.md`` against them. It also catches reversed markdown
links, i.e. ``(text)[url]`` instead of ``[text](url)``.

Exits non-zero on any finding so CI fails loudly.

Usage:
    python ci/check_docs_flags.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Command lines that invoke another tool entirely — their flags are not ours.
FOREIGN_TOOLS = re.compile(r"\b(conda|mamba|pip|docker|singularity|apptainer|foldseek|wget|curl)\b")

CMD = re.compile(r"(^|[`$\s])baktfold\s+[a-z-]+")
FLAG = re.compile(r"(?<![\w-])(--?[A-Za-z][\w-]*)")
REVERSED_LINK = re.compile(r"\(([^()\n]{1,60})\)\[(https?://[^\]\n]+)\]")


def declared_flags() -> set[str]:
    src = (REPO / "src" / "baktfold" / "__init__.py").read_text()
    flags = set(re.findall(r'"(--[a-z0-9-]+)"', src))
    flags |= set(re.findall(r'"(-[a-zA-Z])"', src))
    # click supplies these on every command
    flags |= {"--help", "-h", "--version", "-V"}
    return flags


def doc_files() -> list[Path]:
    files = sorted((REPO / "docs").rglob("*.md"))
    readme = REPO / "README.md"
    if readme.exists():
        files.append(readme)
    return files


def main() -> int:
    declared = declared_flags()
    problems: list[str] = []

    for path in doc_files():
        rel = path.relative_to(REPO)
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            for match in REVERSED_LINK.finditer(line):
                problems.append(
                    f"{rel}:{lineno}: reversed markdown link "
                    f"({match.group(1)})[{match.group(2)}] — should be [text](url)"
                )

            if not CMD.search(line) or FOREIGN_TOOLS.search(line):
                continue

            for flag in FLAG.findall(line):
                if not (flag.startswith("--") or len(flag) == 2):
                    continue
                if flag not in declared:
                    problems.append(
                        f"{rel}:{lineno}: undeclared baktfold flag {flag!r}\n"
                        f"    {line.strip()[:100]}"
                    )

    if problems:
        print(f"{len(problems)} documentation problem(s) found:\n", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        print(
            "\nEvery flag in a documented `baktfold ...` command must be declared "
            "in src/baktfold/__init__.py.",
            file=sys.stderr,
        )
        return 1

    print(f"OK — {len(declared)} declared flags; docs use none that are undeclared.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
