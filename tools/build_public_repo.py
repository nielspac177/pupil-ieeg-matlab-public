"""
Assemble the public, code-only mirror of this repository.

The public repository exists so that the analysis can be read, criticised and
re-run against another dataset. It therefore contains the code, the methods,
and the one figure that is drawn from synthetic traces — and nothing derived
from participants: no result tables, no rendered electrode positions, no
manuscript.

METHODS.md is exported from the built manuscript rather than written here, so
the public methods cannot drift from the methods that were actually submitted.

Run:
    python tools/build_public_repo.py [destination]

The destination defaults to ../pupil-ieeg-public. Anything already there is
replaced, except its .git directory, so the mirror can be re-synced and
committed as a normal update.
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

from docx import Document

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DESTINATION = REPO.parent / "pupil-ieeg-public"
PUBLIC_URL = "https://github.com/nielspac177/pupil-ieeg-matlab-public"

# Copied verbatim. Everything not listed is excluded by construction rather
# than by a filter, so a new results directory cannot leak in by accident.
COPY_TREES = ["src", "tests", "config", "tools", "examples", ".github"]
COPY_FILES = ["run_all.m", "run_tests.m", "LICENSE",
              ".gitignore", "docs/analysis_plan.md", "docs/raw_data_contract.md",
              "docs/supplied_data_inventory.md"]

METHODS_FIGURE = REPO / "results" / "figures" / "Fig0_methods_schematic.png"
MANUSCRIPT = REPO / "manuscript" / "Pupil_MS_v5.docx"

# Sections of the manuscript that constitute the public methods record.
EXPORT_SECTIONS = ["Materials and methods"]
STOP_SECTIONS = ["Results"]


def export_methods() -> str:
    """Pull the Methods section, and the Figure 1 caption, out of the .docx."""
    document = Document(MANUSCRIPT)
    lines: list[str] = []
    capturing = False
    caption = ""

    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        if not text:
            continue
        style = paragraph.style.name
        if text.startswith("Figure 1."):
            caption = text
        if style == "Heading 1":
            if text in EXPORT_SECTIONS:
                capturing = True
                lines.append(f"## {text}")
                continue
            if text in STOP_SECTIONS:
                capturing = False
            if capturing:
                lines.append(f"## {text}")
                continue
        if not capturing:
            continue
        if style == "Heading 2":
            lines.append(f"### {text}")
        elif text.startswith("Figure 1."):
            continue
        else:
            lines.append(text)

    if not lines:
        raise ValueError("no Methods section found in the manuscript")

    header = [
        "# Methods",
        "",
        "Exported verbatim from the submitted manuscript by "
        "`tools/build_public_repo.py`. Edit the manuscript builder in the "
        "private repository, not this file.",
        "",
        "![Measurement and model](figures/Fig0_methods_schematic.png)",
        "",
        f"*{caption}*" if caption else "",
        "",
    ]
    return "\n\n".join(part for part in header + lines if part) + "\n"


PUBLIC_CITATION = """cff-version: 1.2.0
message: "If you use this code, please cite it as below."
title: "pupil-ieeg-matlab: a hurdle-model pipeline for pupil-linked intracranial electrophysiology"
abstract: >-
  MATLAB analysis code for relating spontaneous high-gamma activity to pupil
  diameter in human intracranial recordings. Implements a two-part hurdle model
  of a semicontinuous coupling outcome, with contacts nested in electrode shafts
  nested in patients, FDR-controlled region families, a within-patient paired
  check, and a selection-threshold sensitivity analysis. Code and methods only;
  no participant data.
type: software
license: MIT
repository-code: "{repository}"
keywords:
  - pupillometry
  - intracranial EEG
  - high-gamma
  - arousal
  - hippocampus
  - mixed-effects models
  - reproducibility
authors:
  - family-names: Pacheco-Barrios
    given-names: Niels
    affiliation: "Department of Neurosurgery, Brigham and Women's Hospital, Boston, MA, USA"
"""

PUBLIC_README = """# pupil-ieeg-matlab (public mirror)

MATLAB analysis code for a study of how spontaneous high-gamma activity in the
human brain relates to pupil diameter, from the Rolston Lab.

**This mirror contains code and methods only.** It contains no participant data,
no result tables, and no rendered electrode positions. The recordings are
covered by consent and institutional approvals that do not authorise public
sharing, so `run_all` cannot be executed from this repository alone.

## What you can run here

```matlab
run_demo     % methods schematic + synthetic ripple detector; needs no data
```

`run_demo` draws [`METHODS.md`](METHODS.md)'s figure from synthetic traces and
exercises the ripple detector on a simulated signal. Neither uses patient data.

```matlab
run_tests    % the unit tests; those that load the derived dataset will skip
run_all      % the full pipeline; requires the derived dataset (not public)
```

## What the analysis does

Two properties of the source data determine the whole design, and both are
easy to get wrong:

1. The stored coupling statistic is **not a p-value.** It is the fraction of
   the peri-peak window spanned by the longest contiguous run above a surrogate
   band — a contiguity criterion with no calibrated null distribution.
   Thresholding it *selects* contacts; it does not identify significant ones.
2. The resulting outcome is **semicontinuous**: exactly zero at any contact
   where no suprathreshold run occurred, which in the study cohort was 69% of
   them. A single regression on the pooled variable is misspecified.

The pipeline therefore fits a two-part hurdle model — prevalence, then
direction, then magnitude — with contacts nested in electrode shafts nested in
patients, FDR-controlled region families, a within-patient paired check, and a
sensitivity analysis across selection thresholds. See [`METHODS.md`](METHODS.md)
and [`docs/analysis_plan.md`](docs/analysis_plan.md).

## Reusing it on your own data

`src/+phg/loadDerivedTables.m` defines the contract the pipeline expects of an
input table: one row per contact, with a signed peri-peak response, a region
label, a patient identifier and a contact label from which the electrode shaft
can be parsed. `docs/raw_data_contract.md` states it in full. If you can put
your own recordings into that shape, the models in
`src/+phg/runDerivedAnalyses.m` apply unchanged.

## Requirements

MATLAB R2024b or newer with the Statistics and Machine Learning, Signal
Processing, and Image Processing Toolboxes.
[LeGUI](https://github.com/Rolston-Lab/LeGUI) is an optional GPL-3.0 dependency
for template-surface generation and is not vendored here.

## Citation

See [`CITATION.cff`](CITATION.cff) to cite this code. The study these methods
were developed for is a multi-author manuscript currently under review; this
section will be updated with its reference on publication, and that paper —
not this repository — is the citation for the scientific findings.

## Licence

MIT — see [`LICENSE`](LICENSE).
"""

RUN_DEMO = """%% RUN_DEMO Everything in this repository that runs without participant data.
% The full pipeline (run_all) needs the derived dataset, which is not public.
% This script exercises the parts that do not: the methods schematic, which is
% drawn from synthetic traces, and the ripple detector, on a simulated signal.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'config'));

cfg = default_config(repoRoot);
phg.ensureProjectDirectories(cfg);

fprintf('[PHG] Methods schematic (synthetic traces)\\n');
phg.makeMethodsFigure(cfg);

fprintf('[PHG] Ripple detector demonstration (simulated signal)\\n');
run(fullfile(repoRoot, 'examples', 'demo_ripple_detector.m'));

fprintf('[PHG] Demo complete. Output in %s\\n', cfg.figureDir);
"""


def main() -> int:
    destination = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 \
        else DEFAULT_DESTINATION
    if not METHODS_FIGURE.exists():
        raise FileNotFoundError(
            f"{METHODS_FIGURE} is missing. Run phg.makeMethodsFigure first.")

    preserved_git = destination / ".git"
    if destination.exists():
        for entry in destination.iterdir():
            if entry == preserved_git:
                continue
            shutil.rmtree(entry) if entry.is_dir() else entry.unlink()
    else:
        destination.mkdir(parents=True)

    # local_config.m holds machine-local absolute paths. It is gitignored in
    # the private repository, but copytree does not consult .gitignore, so it
    # is excluded here by name and the result is checked below.
    ignore_names = shutil.ignore_patterns("local_config.m", "patient_task_map.tsv",
                                          "*.asv", ".DS_Store")
    for tree in COPY_TREES:
        source = REPO / tree
        if source.exists():
            shutil.copytree(source, destination / tree, ignore=ignore_names)
    for name in COPY_FILES:
        source = REPO / name
        if source.exists():
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)

    (destination / "figures").mkdir(exist_ok=True)
    shutil.copy2(METHODS_FIGURE, destination / "figures" / METHODS_FIGURE.name)

    (destination / "CITATION.cff").write_text(
        PUBLIC_CITATION.format(repository=PUBLIC_URL))
    (destination / "METHODS.md").write_text(export_methods())
    (destination / "README.md").write_text(PUBLIC_README)
    (destination / "run_demo.m").write_text(RUN_DEMO)

    # A results directory would be created by the first run; keep it empty and
    # ignored so nothing derived is ever tracked here.
    ignore = destination / ".gitignore"
    ignore.write_text(ignore.read_text() + "\nresults/\nmanuscript/\nsubmission/\n")

    leaked = [p for p in destination.rglob("*")
              if p.is_file() and (p.suffix in {".csv", ".docx"}
                                  or p.name == "local_config.m")
              and ".git" not in p.parts and p.name != "raw_manifest.tsv"]
    if leaked:
        raise RuntimeError(f"derived or local files leaked into the public tree: {leaked}")

    # Nothing in a public repository should name a directory on this machine.
    pattern = re.compile(r"(/Users/|/home/)[A-Za-z0-9._-]+/")
    for path in destination.rglob("*"):
        if not path.is_file() or ".git" in path.parts or path.suffix == ".png":
            continue
        try:
            text = path.read_text(errors="replace")
        except (OSError, UnicodeDecodeError):
            continue
        for number, line in enumerate(text.splitlines(), 1):
            if "/absolute/path/to/" in line or "example" in path.name.lower():
                continue
            if pattern.search(line):
                raise RuntimeError(
                    f"absolute home path in public tree: "
                    f"{path.relative_to(destination)}:{number}")

    files = [p for p in destination.rglob("*") if p.is_file() and ".git" not in p.parts]
    total = sum(p.stat().st_size for p in files)
    print(f"Wrote {destination}: {len(files)} files, {total/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
