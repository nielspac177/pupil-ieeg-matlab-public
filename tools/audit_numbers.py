"""
Audit every number in the manuscript against the result tables.

Two axes, both of which must pass:

  coherence  does the paper agree with itself? Every N must descend from one
             cohort ledger, and a count stated in the abstract, the methods,
             a caption and a table must be the same count.
  fidelity   does the paper agree with the code? Every effect estimate must
             match the CSV the model wrote.

The manuscript is generated from those CSVs, so fidelity should hold by
construction — which is exactly why it is worth testing. Any literal digit in
the builder is a place where that guarantee was bypassed, and this script
finds them.

Exit code 1 on any FAIL, so it can gate a submission step.

Run:
    python tools/audit_numbers.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document

REPO = Path(__file__).resolve().parent.parent
TABLES = REPO / "results" / "tables"
MAIN = REPO / "manuscript" / "Pupil_MS_v5.docx"
SUPPLEMENT = REPO / "manuscript" / "Pupil_MS_v5_supplement.docx"
BUILDER = REPO / "manuscript" / "build_v5_docx.py"

TOLERANCE = {"count": 0.0, "percent": 0.1, "ratio": 0.01, "coefficient": 0.01}

failures: list[str] = []
unmatched: list[str] = []
suspicious: list[str] = []
passes = 0


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name)


def document_text(path: Path) -> str:
    document = Document(path)
    parts = [p.text for p in document.paragraphs]
    for table in document.tables:
        for row in table.rows:
            parts.extend(cell.text for cell in row.cells)
    return "\n".join(parts)


def check(label: str, kind: str, reported, computed) -> None:
    """Record one comparison against its tolerance."""
    global passes
    if reported is None:
        unmatched.append(label)
        return
    if abs(float(reported) - float(computed)) <= TOLERANCE[kind]:
        passes += 1
    else:
        failures.append(
            f"{label}: paper {reported}, computed {computed} "
            f"(tolerance {TOLERANCE[kind]})")


def main() -> int:
    global passes

    if not MAIN.exists():
        print(f"manuscript missing: {MAIN}", file=sys.stderr)
        return 1

    # Freshness: a manuscript older than the tables it quotes is a stale build.
    newest_table = max(p.stat().st_mtime for p in TABLES.glob("*.csv"))
    if MAIN.stat().st_mtime < newest_table:
        failures.append(
            "manuscript is older than the result tables; rebuild before auditing")

    text = document_text(MAIN) + "\n" + document_text(SUPPLEMENT)

    # ---------------------------------------------------------------- ledger
    design = load("analysis_design_summary.csv").set_index("parameter")["value"]
    audit = load("derived_data_audit.csv").set_index("metric")["value"]
    coverage = load("region_coverage_summary.csv")

    ledger = {
        "contacts": int(design["n_contacts"]),
        "excursion contacts": int(design["n_excursion_contacts"]),
        "patients": int(design["n_patients"]),
        "shafts": int(design["n_shafts"]),
        "legacy-selected contacts": int(audit["legacy_selected_channels"]),
        "informative paired patients": int(design["paired_patients_informative"]),
    }
    ledger["zero contacts"] = ledger["contacts"] - ledger["excursion contacts"]

    print("Cohort ledger")
    for name, value in ledger.items():
        occurrences = len(re.findall(rf"\b{value}\b", text))
        print(f"  {name:<30} {value:>5}   appears {occurrences}x in the paper")
        if occurrences == 0:
            unmatched.append(f"ledger value {name}={value} never appears")

    # Coherence: the two halves of the hurdle must add up wherever stated.
    if ledger["excursion contacts"] + ledger["zero contacts"] != ledger["contacts"]:
        failures.append("hurdle split does not sum to the contact total")
    else:
        passes += 1

    stated_percent = re.search(
        rf"{ledger['zero contacts']} of {ledger['contacts']} contacts "
        r"\((\d+)%\)", text)
    if stated_percent:
        check("methods zero-mass percentage", "percent",
              float(stated_percent.group(1)),
              round(100 * ledger["zero contacts"] / ledger["contacts"]))
    else:
        unmatched.append("methods zero-mass percentage sentence not found")

    excursion_percent = re.search(
        rf"{ledger['excursion contacts']} of {ledger['contacts']} contacts "
        r"\((\d+)%\)", text)
    if excursion_percent:
        check("abstract excursion percentage", "percent",
              float(excursion_percent.group(1)),
              round(100 * ledger["excursion contacts"] / ledger["contacts"]))

    # Region coverage totals must reconcile with the ledger.
    if int(coverage["n_contacts"].sum()) != ledger["contacts"]:
        failures.append(
            f"region coverage sums to {int(coverage['n_contacts'].sum())}, "
            f"ledger says {ledger['contacts']}")
    else:
        passes += 1

    # --------------------------------------------------------------- fidelity
    def find_ratio(pattern: str) -> float | None:
        match = re.search(pattern, text)
        return float(match.group(1)) if match else None

    polarity = load("polarity_primary_glme.csv")
    effect = polarity[polarity["term"].str.contains("Hippocampal")].iloc[0]
    check("primary odds ratio", "ratio",
          find_ratio(r"odds ratio for dilation of (\d+\.\d+)"),
          effect["odds_ratio"])
    check("primary CI low", "ratio",
          find_ratio(r"odds ratio for dilation of \d+\.\d+ \(\[(\d+\.\d+)"),
          effect["odds_ratio_ci95_low"])

    # Region estimates are audited against the manuscript TABLE, not the prose.
    # A generated table cell has unambiguous provenance; a sentence does not,
    # and scraping prose for "hippocampus ... 0.70" happily matches a magnitude
    # coefficient three sentences away from the prevalence odds ratio it meant.
    prevalence = load("excursion_prevalence_glme.csv").assign(
        region=lambda d: d["term"].str.replace("Region_", "", regex=False)
    ).set_index("region")
    magnitude = load("signed_magnitude_mixed_model.csv").assign(
        region=lambda d: d["term"].str.replace("Region_", "", regex=False)
    ).set_index("region")

    region_table = None
    for table in Document(MAIN).tables:
        header = [cell.text.strip() for cell in table.rows[0].cells]
        if any("Prevalence OR" in cell for cell in header):
            region_table = table
            break
    if region_table is None:
        unmatched.append("region results table not found in the manuscript")
    else:
        for row in region_table.rows[1:]:
            cells = [cell.text.strip() for cell in row.cells]
            region = cells[0]
            if region not in prevalence.index:
                unmatched.append(f"table row '{region}' has no model row")
                continue
            check(f"Table 2 contacts, {region}", "count",
                  float(cells[1]), prevalence.loc[region, "df"] * 0 +
                  float(coverage.set_index("region").loc[region, "n_contacts"]))
            check(f"Table 2 prevalence OR, {region}", "ratio",
                  float(cells[4].split()[0]), prevalence.loc[region, "odds_ratio"])
            beta = cells[6].split()[0].replace("\u2212", "-").replace("+", "")
            check(f"Table 2 magnitude beta, {region}", "coefficient",
                  float(beta), magnitude.loc[region, "estimate"])

    sensitivity = load("polarity_threshold_sensitivity.csv")
    robust = sensitivity[sensitivity["p_value"] < 0.05]
    bounds = re.search(
        r"leaves the odds ratio between (\d+\.\d+) and (\d+\.\d+)", text)
    # Guard against a future edit reintroducing a shared stem: two analyses
    # reporting an odds-ratio range must not be indistinguishable to a reader
    # or to this audit.
    if len(re.findall(r"the odds ratio (?:between|within)", text)) > 2:
        suspicious.append(
            "more than two sentences report an odds-ratio range with the same "
            "stem; disambiguate them")
    if bounds:
        check("sensitivity lower bound", "ratio", float(bounds.group(1)),
              robust["odds_ratio_hippocampal"].min())
        check("sensitivity upper bound", "ratio", float(bounds.group(2)),
              robust["odds_ratio_hippocampal"].max())
    else:
        unmatched.append("threshold-sensitivity range sentence not found")

    leave_one_out = load("polarity_leave_one_patient_out.csv").dropna(
        subset=["odds_ratio_hippocampal"])
    lopo = re.search(r"holds the odds ratio within (\d+\.\d+) and (\d+\.\d+) "
                     r"across all (\d+) refits", text)
    if lopo:
        check("leave-one-out lower", "ratio", float(lopo.group(1)),
              leave_one_out["odds_ratio_hippocampal"].min())
        check("leave-one-out upper", "ratio", float(lopo.group(2)),
              leave_one_out["odds_ratio_hippocampal"].max())
        check("leave-one-out refit count", "count", float(lopo.group(3)),
              len(leave_one_out))
    else:
        unmatched.append("leave-one-patient-out sentence not found")

    # --------------------------------------------------- hardcoded literals
    # Any digit typed into the builder bypasses the CSV guarantee. Acquisition
    # constants legitimately live there; anything the tables also know must not.
    source = BUILDER.read_text()
    derivable = {str(v) for v in ledger.values()}
    for line_number, line in enumerate(source.splitlines(), 1):
        if (line.lstrip().startswith("#") or "f\"" in line or "f'" in line
                or "http" in line):
            continue
        for literal in re.findall(r'"[^"]*?(\d{2,4})[^"]*?"', line):
            if literal in derivable:
                suspicious.append(
                    f"build_v5_docx.py:{line_number}: literal {literal} is a "
                    f"ledger value and should be interpolated, not typed")

    # -------------------------------------------------- duplicate estimates
    ratios = re.findall(r"odds ratio (\d+\.\d{2,3})", text)
    for value in set(ratios):
        if ratios.count(value) > 2:
            suspicious.append(
                f"odds ratio {value} appears {ratios.count(value)} times — "
                "check these are genuinely distinct estimates")

    # ------------------------------------------------------------- report
    print(f"\nPASS {passes}   FAIL {len(failures)}   "
          f"UNMATCHED {len(unmatched)}   SUSPICIOUS {len(suspicious)}")
    for title, items in [("FAIL", failures), ("UNMATCHED", unmatched),
                         ("SUSPICIOUS", suspicious)]:
        if items:
            print(f"\n{title}")
            for item in items:
                print(f"  - {item}")

    print("\nOverall:", "FAIL" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
