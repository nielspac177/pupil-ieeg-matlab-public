#!/usr/bin/env python3
"""
Static checks that do not need a MATLAB licence.

CI cannot run the pipeline — MATLAB is licensed software and the input data are
not public — so it checks the things that break silently and are cheap to
verify without executing anything:

  1. every `phg.<name>(` referenced anywhere in the project resolves to a file
     in src/+phg/, which catches a renamed or deleted function that no test
     happens to exercise;
  2. every function file in src/+phg/ declares the function name its filename
     promises, which MATLAB would otherwise resolve to something surprising;
  3. .gitignore still forbids the raw-recording extensions, so a future edit
     cannot quietly open the door to committing patient data;
  4. no tracked file contains an absolute path under /Users or /home.

Exit code 0 if everything passes, 1 otherwise, with every failure printed.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "src" / "+phg"

RAW_EXTENSIONS = [
    "*.ns[1-9]", "*.nev", "*.nf3", "*.pldata", "*.edf", "*.eeg",
    "*.fdt", "*.set", "*.nwb", "*.mefd/", "*.nii", "*.nii.gz",
]

failures: list[str] = []


def tracked_files() -> list[Path]:
    try:
        output = subprocess.run(
            ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True,
            check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return sorted(ROOT.rglob("*.m"))
    return [ROOT / line for line in output.splitlines() if line]


def check_package_references() -> None:
    available = {path.stem for path in PACKAGE.glob("*.m")}
    pattern = re.compile(r"\bphg\.([A-Za-z]\w*)\s*\(")
    for path in tracked_files():
        if path.suffix != ".m" or not path.exists():
            continue
        text = path.read_text(errors="replace")
        for line_number, line in enumerate(text.splitlines(), 1):
            if line.lstrip().startswith("%"):
                continue
            for name in pattern.findall(line):
                if name not in available:
                    failures.append(
                        f"{path.relative_to(ROOT)}:{line_number}: "
                        f"phg.{name} has no file in src/+phg/")


def check_function_names() -> None:
    for path in sorted(PACKAGE.glob("*.m")):
        first = ""
        for line in path.read_text(errors="replace").splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("%"):
                first = stripped
                break
        match = re.match(r"function\s+(?:.*?=\s*)?(\w+)\s*\(", first)
        if not match:
            failures.append(f"{path.name}: first statement is not a function definition")
        elif match.group(1) != path.stem:
            failures.append(
                f"{path.name}: declares function {match.group(1)}, "
                f"filename promises {path.stem}")


def check_gitignore() -> None:
    ignore = ROOT / ".gitignore"
    if not ignore.exists():
        failures.append(".gitignore is missing")
        return
    text = ignore.read_text()
    for pattern in RAW_EXTENSIONS:
        if pattern not in text:
            failures.append(f".gitignore no longer excludes {pattern}")


def check_absolute_paths() -> None:
    pattern = re.compile(r"(/Users/|/home/)[A-Za-z0-9._-]+/")
    for path in tracked_files():
        if not path.exists() or path.suffix in {".png", ".pdf", ".docx", ".tiff"}:
            continue
        try:
            text = path.read_text(errors="replace")
        except (OSError, UnicodeDecodeError):
            continue
        for line_number, line in enumerate(text.splitlines(), 1):
            if "example" in path.name.lower() or "/absolute/path/to/" in line:
                continue
            if pattern.search(line):
                failures.append(
                    f"{path.relative_to(ROOT)}:{line_number}: "
                    "absolute home path in a tracked file")


def main() -> int:
    check_package_references()
    check_function_names()
    check_gitignore()
    check_absolute_paths()
    if failures:
        print(f"{len(failures)} problem(s):")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("structure checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
