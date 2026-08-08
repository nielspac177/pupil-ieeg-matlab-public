#!/usr/bin/env python3
"""Re-verify the pupil-to-iEEG time anchor for every staged replication run.

The check required by docs/replication_amendment_01.md section 2, run as its own
auditable step rather than only as a side effect of ingestion. Each MEF3 session
is reopened for metadata and records alone -- no signal is decoded -- in its own
subprocess, because pymef segfaults when many sessions are opened in one
interpreter.

For every run it recomputes the anchor from the hardware record stamps,
validates it by asking whether every record lands on a row of events.tsv, and
reports how far the *staged* anchor sits from the recomputed one. That last
column is the point of the tool: a staged anchor that differs by less than one
sample changes nothing, whereas one that differs by seconds means the run was
epoched against the wrong part of the recording and must be re-ingested.

Writes results/tables/replication_time_anchor_audit.csv and patches
anchor_valid in each staged meta.json.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

TOLERANCE_US = 1000.0
MATCH_FRACTION_REQUIRED = 0.95

WORKER = r"""
import csv, json, os, sys
import numpy as np
from pymef.mef_session import MefSession

mefd, events_path = sys.argv[1], sys.argv[2]
with open(events_path, newline="") as handle:
    onsets = np.array([float(r["onset"]) for r in
                       csv.DictReader(handle, delimiter="\t")], dtype=float)
session = MefSession(mefd, None)
record_times = np.array([r["time"] for r in session.read_records()],
                        dtype=np.int64)

common = min(record_times.size, onsets.size)
anchor = float(np.median(record_times[:common] - onsets[:common] * 1e6))

relative = (record_times - anchor) / 1e6
ordered = np.sort(onsets)
position = np.clip(np.searchsorted(ordered, relative), 1, ordered.size - 1)
error = np.minimum(np.abs(relative - ordered[position - 1]),
                   np.abs(relative - ordered[position]))

info = session.read_ts_channel_basic_info()
seeg = [c for c in info if c["channel_description"][0] == b"SEEG"]
signal_anchor = float("nan")
if seeg:
    c = seeg[0]
    signal_anchor = (float(c["end_time"][0])
                     - float(c["nsamp"][0]) / float(c["fsamp"][0]) * 1e6)

print(json.dumps({
    "n_events": int(onsets.size),
    "n_records": int(record_times.size),
    "anchor_uutc": anchor,
    "signal_anchor_uutc": signal_anchor,
    "match_fraction": float(np.mean(error <= 1e-3)),
    "median_error_us": float(np.median(error) * 1e6),
    "max_error_us": float(np.max(error) * 1e6),
}))
"""


def probe(mefd: Path, events: Path) -> dict | None:
    try:
        finished = subprocess.run(
            [sys.executable, "-c", WORKER, str(mefd), str(events)],
            capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return None
    if finished.returncode != 0 or not finished.stdout.strip():
        return None
    return json.loads(finished.stdout.strip().splitlines()[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # No machine-local defaults: this file is published in the public mirror.
    parser.add_argument("--dataset", default=os.environ.get("EBRAINS_DATASET"))
    parser.add_argument("--stage", default=os.environ.get("EBRAINS_STAGE"))
    parser.add_argument("--out", default="results/tables/"
                        "replication_time_anchor_audit.csv")
    args = parser.parse_args()

    if not args.dataset or not args.stage:
        parser.error("--dataset and --stage are required "
                     "(or set EBRAINS_DATASET and EBRAINS_STAGE)")

    dataset, stage = Path(args.dataset), Path(args.stage)
    rows = []
    for stage_dir in sorted(stage.glob("sub-*")):
        meta_path = stage_dir / "meta.json"
        if not meta_path.is_file():
            continue
        meta = json.loads(meta_path.read_text())
        subject, task, run = meta["subject"], meta["task"], meta["run"]
        ieeg_dir = dataset / subject / "ses-001" / "ieeg"
        stem = f"{subject}_ses-001_task-{task}_run-{run}"

        probed = probe(ieeg_dir / f"{stem}_ieeg.mefd",
                       ieeg_dir / f"{stem}_events.tsv")
        if probed is None:
            rows.append({"tag": meta["tag"], "subject": subject, "task": task,
                         "run": run, "n_events": "", "n_records": "",
                         "match_fraction": "", "median_error_us": "",
                         "max_error_us": "", "staged_minus_recomputed_us": "",
                         "staged_offset_samples_at_1khz": "",
                         "anchor_valid": False, "note": "probe_failed"})
            print(f"[anchor] {meta['tag']:28} PROBE FAILED", flush=True)
            continue

        staged = float(meta["anchor_uutc"])
        delta_us = staged - probed["anchor_uutc"]
        valid = bool(probed["match_fraction"] >= MATCH_FRACTION_REQUIRED
                     and probed["median_error_us"] <= TOLERANCE_US
                     and abs(delta_us) <= TOLERANCE_US)

        meta["anchor_valid"] = valid
        meta["anchor_record_match_fraction"] = probed["match_fraction"]
        meta["anchor_staged_minus_recomputed_us"] = delta_us
        meta_path.write_text(json.dumps(meta, indent=2))

        rows.append({
            "tag": meta["tag"], "subject": subject, "task": task, "run": run,
            "n_events": probed["n_events"], "n_records": probed["n_records"],
            "match_fraction": round(probed["match_fraction"], 6),
            "median_error_us": round(probed["median_error_us"], 3),
            "max_error_us": round(probed["max_error_us"], 3),
            "staged_minus_recomputed_us": round(delta_us, 3),
            "staged_offset_samples_at_1khz": round(abs(delta_us) / 1000.0, 6),
            "anchor_valid": valid, "note": "",
        })
        print(f"[anchor] {meta['tag']:28} ev={probed['n_events']:5d} "
              f"rec={probed['n_records']:5d} matched={probed['match_fraction']:.3f} "
              f"staged-recomputed={delta_us:+.1f} us -> "
              f"{'VALID' if valid else 'INVALID'}", flush=True)

    if not rows:
        print("[anchor] no staged runs found")
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    valid = sum(1 for r in rows if r["anchor_valid"])
    print(f"[anchor] {valid}/{len(rows)} runs pass; wrote {out}")
    return 0 if valid == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
