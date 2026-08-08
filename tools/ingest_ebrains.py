#!/usr/bin/env python3
"""Decode the EBRAINS MEF3 replication cohort into MATLAB-readable staging files.

This is an I/O step, not an analysis step. It decodes MEF3, establishes the
shared time base, applies the prespecified inter-trial epoch restriction, and
decimates the wideband signal to a rate that fully resolves the 70-170 Hz
high-gamma band and the 80-120 Hz ripple band. Every measurement, every model
and every statistical claim happens afterwards in MATLAB.

Protocol: docs/replication_plan_ebrains.md, as amended by
docs/replication_amendment_01.md. Nothing here may be changed on the basis of a
result; the exclusions and the time-anchor validation are fixed in those
documents.

Usage:
    export EBRAINS_DATASET=/path/to/f2f41456-...-v1.0
    export EBRAINS_STAGE=/path/to/staging
    work/venv/bin/python tools/ingest_ebrains.py
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy import signal

# --------------------------------------------------------------------------
# Prespecified constants. See docs/replication_amendment_01.md.
# --------------------------------------------------------------------------

PRIMARY_TASKS = ("FR", "PAL")           # SP and AP are eye-movement tasks
SENSITIVITY_TASKS = ("SP", "AP")        # reported separately, never pooled
EVENT_GUARD_SECONDS = 2.0               # +/- 2 s clear of any task event
FILTER_PAD_SECONDS = 5.0                # discarded context for filter edges

# The discovery cohort was measured at exactly 1000 Hz with
# butter(1, [70 170]/500). Reproducing that filter response requires the same
# sampling rate, so the wideband path is decimated to 1000 Hz rather than to a
# convenient higher rate. Verified against the stored analysis object:
# obj.Fs = 1000, obj.BandPass = [70 170], obj.ThrX = 5, obj.TimeRng = +/-20000 ms.
TARGET_FS = 1000.0
# Hippocampal contacts are additionally staged at 4 kHz for ripple detection,
# where sampling headroom above the 80-120 Hz band is worth the disk.
RIPPLE_FS = 4000.0
PUPIL_FS = 150.0
ANCHOR_TOLERANCE_US = 1000.0            # records vs SEEG derivation, 1 ms
HIPPOCAMPUS_LABEL = "hippocampus"

# Blink / artifact handling, fixed before any outcome was computed.
BLINK_SHOULDER_SECONDS = 0.100
OUTLIER_MAD = 5.0
RUNNING_MEDIAN_SECONDS = 30.0
MAX_INTERPOLATED_GAP_SECONDS = 0.300

GAZE_CHANNELS = ("LEFT_X_COORD", "LEFT_Y_COORD", "RIGHT_X_COORD",
                 "RIGHT_Y_COORD", "SCREEN_X", "SCREEN_Y")
PUPIL_CHANNELS = ("LEFT_PUPIL_SIZE", "RIGHT_PUPIL_SIZE")

# Typographic repair only. No label is merged across distinct structures.
LABEL_CANON = {
    "middle temporalgyrus": "middle temporal gyrus",
    "suprior frontal gyrus": "superior frontal gyrus",
    "superior frontal gyrus (mesial)": "superior frontal gyrus",
    "post cingulate": "posterior cingulate",
    "mid cingulate": "middle cingulate",
    "parietal operculum (postcentral)": "parietal operculum",
    "planum polare / insula": "planum polare",
}
LESIONAL = re.compile(r"encephalomalacia|encphalomalacia|heterotopia", re.I)
NON_ANATOMICAL = {"wm", "n/a", ""}


def canonical_label(raw: str) -> str:
    """Collapse whitespace and repair the known typographic variants."""
    collapsed = " ".join(str(raw).strip().lower().split())
    return LABEL_CANON.get(collapsed, collapsed)


# --------------------------------------------------------------------------
# Metadata readers
# --------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def discover_runs(dataset: Path, subjects: list[str]) -> list[dict]:
    """Enumerate every run with its sidecar files, in deterministic order."""
    runs = []
    for subject in subjects:
        ieeg_dir = dataset / subject / "ses-001" / "ieeg"
        if not ieeg_dir.is_dir():
            continue
        for mefd in sorted(ieeg_dir.glob("*_ieeg.mefd")):
            if mefd.name.startswith("._"):
                continue
            stem = mefd.name[: -len("_ieeg.mefd")]
            task = stem.split("task-")[1].split("_")[0]
            run = stem.split("run-")[1].split("_")[0]
            runs.append({
                "subject": subject,
                "task": task,
                "run": run,
                "mefd": mefd,
                "events": ieeg_dir / f"{stem}_events.tsv",
                "channels": ieeg_dir / f"{stem}_channels.tsv",
                "electrodes": ieeg_dir / f"{subject}_ses-001_electrodes.tsv",
            })
    return runs


def contact_index(name: str) -> float:
    """Trailing contact number on a shaft, used to find Laplacian neighbours."""
    match = re.search(r"(\d+)\s*$", name)
    return float(match.group(1)) if match else float("nan")


def eligible_contacts(dataset: Path, subject: str,
                      tasks: tuple[str, ...]) -> tuple[dict, dict, dict]:
    """Apply the prespecified contact eligibility rule.

    Returns three things: the analysis frame (gray-matter, good, labelled), the
    wider set of good SEEG channels, and a tally of exclusions.

    The two sets differ on purpose. The discovery cohort rereferences each
    contact against its immediate neighbours on the same shaft
    (`PupilHG.m:323-345`), and a neighbour only has to be a good recording
    channel -- it may well sit in white matter. Restricting the neighbour pool
    to the analysis frame would silently change the reference scheme and make
    the two cohorts incomparable, so every good SEEG channel is staged.
    """
    ieeg_dir = dataset / subject / "ses-001" / "ieeg"
    electrodes = {}
    for row in read_tsv(ieeg_dir / f"{subject}_ses-001_electrodes.tsv"):
        electrodes[row["name"]] = {
            "region": canonical_label(row["anatomy_structure"]),
            "shaft": row["group"],
            "index": contact_index(row["name"]),
            "x": row["x"], "y": row["y"], "z": row["z"],
        }

    good: set[str] | None = None
    for path in sorted(ieeg_dir.glob("*_channels.tsv")):
        if path.name.startswith("._"):
            continue
        task = path.name.split("task-")[1].split("_")[0]
        if task not in tasks:
            continue
        present = {row["name"] for row in read_tsv(path)
                   if row["type"] == "SEEG" and row["status"] == "good"}
        good = present if good is None else (good & present)
    good = good or set()

    keep, reference_pool = {}, {}
    dropped = {"bad_or_absent": 0, "wm_or_na": 0, "lesional": 0}
    for name, info in electrodes.items():
        if name not in good:
            dropped["bad_or_absent"] += 1
            continue
        reference_pool[name] = info
        if info["region"] in NON_ANATOMICAL:
            dropped["wm_or_na"] += 1
        elif LESIONAL.search(info["region"]):
            dropped["lesional"] += 1
        else:
            keep[name] = info
    return keep, reference_pool, dropped


# --------------------------------------------------------------------------
# Time anchor -- see amendment section 2
# --------------------------------------------------------------------------

def session_anchor(session, events: list[dict]) -> dict:
    """Recover the absolute time of session-relative t = 0, and validate it.

    The dataset's own reference notebook uses `earliest_start_time`, which is 0
    in these sessions because the SEEG channels carry `start_time = 0`. Using it
    would silently misplace every event by decades. Two independent derivations
    are computed instead and required to agree.
    """
    onsets = np.array([float(row["onset"]) for row in events], dtype=float)
    records = session.read_records()
    record_times = np.array([r["time"] for r in records], dtype=np.int64)

    result = {"n_events": int(onsets.size), "n_records": int(record_times.size)}

    # Primary derivation: the hardware record stamps, paired in file order with
    # the first rows of events.tsv and reduced by a median so that a minority of
    # bad pairs cannot move the estimate.
    #
    # PAL runs carry exactly 90 more rows in events.tsv than there are records --
    # the sidecar keeps derived events the acquisition hardware never stamped --
    # so the two are paired over their common length rather than required to be
    # equal. Verified across all 17 decodable primary runs: this anchor puts
    # every record within 0.0 us of an event onset.
    anchor_from_records = np.nan
    if record_times.size and onsets.size:
        common = min(record_times.size, onsets.size)
        anchor_from_records = float(
            np.median(record_times[:common] - onsets[:common] * 1e6))

    # Advisory cross-check only. A channel's end stamp less its sample-derived
    # duration assumes the recording is gap-free, and it is not: on
    # sub-006 task-PAL run-01 this derivation lands 81 seconds away from the
    # truth, because the SEEG stream has discontinuities the sample count does
    # not see. It is reported, never used.
    info = session.read_ts_channel_basic_info()
    seeg = [c for c in info
            if c["channel_description"][0] == b"SEEG"]
    anchor_from_signal = np.nan
    if seeg:
        channel = seeg[0]
        duration_us = float(channel["nsamp"][0]) / float(channel["fsamp"][0]) * 1e6
        anchor_from_signal = float(channel["end_time"][0]) - duration_us

    # Validation: every record must land on an event once converted.
    matched_fraction = 0.0
    median_error_us = np.inf
    if np.isfinite(anchor_from_records) and record_times.size and onsets.size:
        relative = (record_times - anchor_from_records) / 1e6
        ordered = np.sort(onsets)
        position = np.clip(np.searchsorted(ordered, relative), 1,
                           ordered.size - 1)
        error = np.minimum(np.abs(relative - ordered[position - 1]),
                           np.abs(relative - ordered[position]))
        median_error_us = float(np.median(error) * 1e6)
        matched_fraction = float(np.mean(error <= ANCHOR_TOLERANCE_US / 1e6))

    result.update({
        "anchor_from_records_uutc": anchor_from_records,
        "anchor_from_signal_uutc": anchor_from_signal,
        "signal_minus_record_us": float(anchor_from_signal - anchor_from_records),
        "record_match_fraction": matched_fraction,
        "record_median_error_us": median_error_us,
    })

    if np.isfinite(anchor_from_records) and matched_fraction >= 0.95 \
            and median_error_us <= ANCHOR_TOLERANCE_US:
        result["anchor_uutc"] = anchor_from_records
        result["anchor_source"] = "records"
        result["anchor_valid"] = True
    else:
        result["anchor_uutc"] = anchor_from_records
        result["anchor_source"] = "records_unvalidated"
        result["anchor_valid"] = False
    return result


def eligible_intervals(events: list[dict],
                       guard: float = EVENT_GUARD_SECONDS) -> list[tuple]:
    """Inter-trial intervals at least `guard` seconds clear of any task event."""
    onsets = sorted(float(row["onset"]) for row in events)
    if not onsets:
        return []
    blocked = []
    for onset in onsets:
        lo, hi = onset - guard, onset + guard
        if blocked and lo <= blocked[-1][1]:
            blocked[-1][1] = max(blocked[-1][1], hi)
        else:
            blocked.append([lo, hi])
    first, last = onsets[0], onsets[-1]
    gaps, cursor = [], first
    for lo, hi in blocked:
        if lo > cursor:
            gaps.append((cursor, min(lo, last)))
        cursor = max(cursor, hi)
    if cursor < last:
        gaps.append((cursor, last))
    return [(a, b) for a, b in gaps if b > a]


def _factorise(factor: int) -> list[int]:
    """Split a decimation factor into stages of at most 13.

    scipy warns and loses accuracy on large single-stage decimations; 32000 Hz
    to 1000 Hz is a factor of 32, which becomes 8 then 4.
    """
    if factor <= 1:
        return [1]
    stages, remaining = [], factor
    for step in (13, 11, 8, 7, 5, 4, 3, 2):
        while remaining % step == 0 and remaining > 1:
            stages.append(step)
            remaining //= step
    if remaining > 1:
        stages.append(remaining)
    return sorted(stages, reverse=True) or [1]


def merge_padded(intervals: list[tuple], pad: float) -> list[tuple]:
    """Coalesce padded read blocks so dense regions become one read."""
    merged = []
    for lo, hi in intervals:
        block = [lo - pad, hi + pad]
        if merged and block[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], block[1])
        else:
            merged.append(block)
    return [(a, b) for a, b in merged]


# --------------------------------------------------------------------------
# Pupil conditioning -- see amendment section 5
# --------------------------------------------------------------------------

def running_median(x: np.ndarray, window: int) -> np.ndarray:
    """Median filter that tolerates NaN, used only to locate outliers."""
    if window % 2 == 0:
        window += 1
    padded = np.pad(x, window // 2, mode="edge")
    strided = np.lib.stride_tricks.sliding_window_view(padded, window)
    return np.nanmedian(strided, axis=-1)


def condition_pupil(left: np.ndarray, right: np.ndarray,
                    fs: float) -> tuple[np.ndarray, dict]:
    """Combine eyes, remove blinks and artifacts, interpolate short gaps."""
    eyes = np.vstack([left, right]).astype(float)
    eyes[~np.isfinite(eyes)] = np.nan
    eyes[eyes <= 0] = np.nan
    with np.errstate(invalid="ignore"):
        pupil = np.nanmean(eyes, axis=0)
    pupil[np.all(np.isnan(eyes), axis=0)] = np.nan

    quality = {"raw_missing_fraction": float(np.isnan(pupil).mean())}

    # Blink shoulders: widen every missing sample by +/- 100 ms.
    shoulder = int(round(BLINK_SHOULDER_SECONDS * fs))
    if shoulder > 0 and np.any(np.isnan(pupil)):
        missing = np.isnan(pupil)
        widened = np.convolve(missing.astype(float),
                              np.ones(2 * shoulder + 1), mode="same") > 0
        pupil[widened] = np.nan
    quality["after_shoulder_missing_fraction"] = float(np.isnan(pupil).mean())

    # Robust outlier rejection against a slow running median.
    window = int(round(RUNNING_MEDIAN_SECONDS * fs))
    if np.isfinite(pupil).sum() > window:
        baseline = running_median(pupil, window)
        residual = pupil - baseline
        mad = np.nanmedian(np.abs(residual - np.nanmedian(residual)))
        if np.isfinite(mad) and mad > 0:
            pupil[np.abs(residual) > OUTLIER_MAD * 1.4826 * mad] = np.nan
    quality["after_outlier_missing_fraction"] = float(np.isnan(pupil).mean())

    # Interpolate gaps up to 300 ms; leave longer gaps missing.
    filled = pupil.copy()
    finite = np.isfinite(pupil)
    if finite.any():
        index = np.arange(pupil.size)
        interpolated = np.interp(index, index[finite], pupil[finite])
        max_gap = int(round(MAX_INTERPOLATED_GAP_SECONDS * fs))
        missing = ~finite
        edges = np.diff(np.concatenate([[0], missing.view(np.int8), [0]]))
        for start, stop in zip(np.flatnonzero(edges == 1),
                               np.flatnonzero(edges == -1)):
            if (stop - start) <= max_gap and start > 0 and stop < pupil.size:
                filled[start:stop] = interpolated[start:stop]
    quality["final_missing_fraction"] = float(np.isnan(filled).mean())
    return filled, quality


def regress_out_gaze(pupil: np.ndarray, gaze: dict[str, np.ndarray]) -> np.ndarray:
    """Project the gaze regressors and their derivatives out of the pupil trace.

    The discovery cohort could not do this; it is the reason this dataset was
    chosen. Fitted on samples where the pupil and every regressor are finite,
    then subtracted everywhere.
    """
    columns = []
    for name in GAZE_CHANNELS:
        if name not in gaze:
            continue
        series = gaze[name].astype(float)
        series[~np.isfinite(series)] = np.nan
        columns.append(series)
        columns.append(np.concatenate([[0.0], np.diff(series)]))
    if not columns:
        return pupil.copy()

    design = np.column_stack(columns)
    usable = np.isfinite(pupil) & np.all(np.isfinite(design), axis=1)
    if usable.sum() < 10 * (design.shape[1] + 1):
        return pupil.copy()

    # Centre the regressors on their own means before fitting, so the fitted
    # values are mean-zero and the residual keeps the pupil's original level.
    # Fitting an explicit intercept and then subtracting only its coefficient
    # leaves the regressor means in the residual as a large constant offset --
    # harmless to the peri-peak measurement, which detrends every trial, but it
    # makes the staged trace impossible to read against the raw one.
    centre = design[usable].mean(axis=0)
    centred = design - centre
    offset = pupil[usable].mean()
    beta, *_ = np.linalg.lstsq(centred[usable], pupil[usable] - offset,
                               rcond=None)
    residual = pupil - centred @ beta
    residual[~np.isfinite(pupil)] = np.nan
    return residual


# --------------------------------------------------------------------------
# Main ingestion
# --------------------------------------------------------------------------

def ingest_run(session_class, run: dict, contacts: dict, reference_pool: dict,
               stage: Path, tasks: tuple[str, ...]) -> dict:
    session = session_class(str(run["mefd"]), None)
    events = read_tsv(run["events"])
    anchor = session_anchor(session, events)
    tag = f"{run['subject']}_task-{run['task']}_run-{run['run']}"
    summary = {"tag": tag, **{k: run[k] for k in ("subject", "task", "run")},
               **anchor}

    if not np.isfinite(anchor["anchor_uutc"]):
        summary["status"] = "dropped_no_time_anchor"
        return summary

    anchor_us = anchor["anchor_uutc"]
    gaps = eligible_intervals(events)
    summary["n_eligible_intervals"] = len(gaps)
    summary["eligible_seconds"] = float(sum(b - a for a, b in gaps))
    if not gaps:
        summary["status"] = "dropped_no_eligible_epochs"
        return summary

    info = {c["name"]: c for c in session.read_ts_channel_basic_info()}
    # Stage every good SEEG channel: the analysis frame plus the neighbours the
    # shaft Laplacian needs. MATLAB decides which rows it models.
    names = [n for n in sorted(reference_pool) if n in info]

    # A channel can be listed in electrodes.tsv and absent from the decoded
    # session if the archive was extracted incompletely -- which happened, on
    # sub-007 task-PAL run-01, where a truncated unzip cost the whole RAB
    # hippocampal shaft and half the eye tracking. Silently intersecting the
    # two lists would have dropped eight hippocampal contacts without a word,
    # so the shortfall is counted and any missing *analysis* contact is fatal.
    absent = sorted(set(reference_pool) - set(names))
    missing_analysis = [n for n in absent if n in contacts]
    summary["n_channels_absent_from_session"] = len(absent)
    summary["n_analysis_contacts_absent"] = len(missing_analysis)
    if missing_analysis:
        summary["status"] = ("dropped_incomplete_session: "
                             f"{len(missing_analysis)} analysis contacts "
                             f"absent, e.g. {missing_analysis[:5]}")
        return summary

    if not any(n in contacts for n in names):
        summary["status"] = "dropped_no_eligible_contacts"
        return summary
    hippocampal = [n for n in names
                   if contacts.get(n, {}).get("region") == HIPPOCAMPUS_LABEL]

    # ---------------- pupil and gaze, whole run at 150 Hz ----------------
    first = min(a for a, _ in gaps)
    last = max(b for _, b in gaps)
    read_lo = int(anchor_us + (first - 40.0) * 1e6)
    read_hi = int(anchor_us + (last + 40.0) * 1e6)

    wanted = [n for n in PUPIL_CHANNELS + GAZE_CHANNELS if n in info]
    traces = session.read_ts_channels_uutc(wanted, [read_lo, read_hi])
    series = {}
    for name, values in zip(wanted, traces):
        if values is None:
            continue
        factor = float(info[name]["ufact"][0])
        series[name] = np.asarray(values, dtype=float) * factor

    if not any(n in series for n in PUPIL_CHANNELS):
        summary["status"] = "dropped_no_pupil_channel"
        return summary

    length = max(v.size for v in series.values())
    left = series.get("LEFT_PUPIL_SIZE", np.full(length, np.nan))
    right = series.get("RIGHT_PUPIL_SIZE", np.full(length, np.nan))
    pupil, quality = condition_pupil(left, right, PUPIL_FS)
    gaze = {k: v for k, v in series.items() if k in GAZE_CHANNELS}
    pupil_gaze_regressed = regress_out_gaze(pupil, gaze)
    pupil_time = (read_lo - anchor_us) / 1e6 + np.arange(pupil.size) / PUPIL_FS
    summary.update({f"pupil_{k}": v for k, v in quality.items()})

    # ---------------- wideband, eligible epochs only ----------------
    blocks = merge_padded(gaps, FILTER_PAD_SECONDS)
    source_fs = float(info[names[0]]["fsamp"][0])
    stage_dir = stage / tag
    stage_dir.mkdir(parents=True, exist_ok=True)

    def decimate_to(trace: np.ndarray, target: float) -> np.ndarray:
        """Anti-aliased rate conversion. Factorised so each stage stays <= 13."""
        factor = int(round(source_fs / target))
        out = np.nan_to_num(trace, nan=0.0)
        for step in _factorise(factor):
            out = signal.decimate(out, step, ftype="fir", zero_phase=True)
        return out.astype(np.float32)

    block_signals, block_times = [], []
    ripple_signals, ripple_times = [], []
    for lo, hi in blocks:
        start = int(anchor_us + lo * 1e6)
        stop = int(anchor_us + hi * 1e6)
        raw = session.read_ts_channels_uutc(names, [start, stop])
        scaled = {}
        for name, values in zip(names, raw):
            if values is None:
                continue
            scaled[name] = np.asarray(values, dtype=float) * float(
                info[name]["ufact"][0])

        columns = [decimate_to(scaled[n], TARGET_FS) if n in scaled else None
                   for n in names]
        widths = [c.size for c in columns if c is not None]
        if not widths:
            continue
        width = min(widths)
        block_signals.append(np.column_stack([
            (c[:width] if c is not None else np.zeros(width, np.float32))
            for c in columns]))
        block_times.append(lo + np.arange(width) / TARGET_FS)

        if hippocampal:
            wide = [decimate_to(scaled[n], RIPPLE_FS) if n in scaled else None
                    for n in hippocampal]
            wide_widths = [c.size for c in wide if c is not None]
            if wide_widths:
                wide_width = min(wide_widths)
                ripple_signals.append(np.column_stack([
                    (c[:wide_width] if c is not None
                     else np.zeros(wide_width, np.float32)) for c in wide]))
                ripple_times.append(lo + np.arange(wide_width) / RIPPLE_FS)

    if not block_signals:
        summary["status"] = "dropped_no_wideband"
        return summary

    ieeg = np.vstack(block_signals)
    ieeg_time = np.concatenate(block_times)

    # Core mask: samples inside an eligible epoch, i.e. not filter padding.
    core = np.zeros(ieeg_time.size, dtype=bool)
    for lo, hi in gaps:
        core |= (ieeg_time >= lo) & (ieeg_time < hi)

    # Sensitivity masks B and C, prespecified in the amendment.
    def window_mask(pre: float, post: float) -> np.ndarray:
        mask = np.zeros(ieeg_time.size, dtype=bool)
        for lo, hi in gaps:
            if (hi - post) > (lo + pre):
                mask |= (ieeg_time >= lo + pre) & (ieeg_time < hi - post)
        return mask

    # HDF5 rather than .npz: MATLAB reads it natively with h5read, and it
    # streams arrays too large for a v5 .mat file. Note that HDF5 is row-major,
    # so a samples-by-channels array here arrives in MATLAB transposed;
    # loadReplicationRun.m transposes it back.
    payload = dict(
        ieeg=ieeg, ieeg_time=ieeg_time.astype(np.float32),
        core_mask=core.astype(np.uint8),
        mask_rule_b=window_mask(5.0, 5.0).astype(np.uint8),
        mask_rule_c=window_mask(20.0, 5.0).astype(np.uint8),
        pupil=pupil.astype(np.float32),
        pupil_gaze_regressed=pupil_gaze_regressed.astype(np.float32),
        pupil_time=pupil_time.astype(np.float32),
        event_onsets=np.array([float(r["onset"]) for r in events],
                              dtype=np.float32),
    )
    if ripple_signals:
        payload["ripple_ieeg"] = np.vstack(ripple_signals)
        payload["ripple_time"] = np.concatenate(ripple_times).astype(np.float32)

    # Stored uncompressed on purpose: these are float32 wideband traces, where
    # gzip returned 7% on a test run for a large share of the wall clock.
    with h5py.File(stage_dir / "run.h5", "w") as handle:
        for key, value in payload.items():
            handle.create_dataset(key, data=value)

    def describe(name: str) -> dict:
        source = contacts.get(name) or reference_pool[name]
        return {"name": name, "region": source["region"],
                "shaft": source["shaft"], "contact_index": source["index"],
                "x": source["x"], "y": source["y"], "z": source["z"],
                "in_analysis_frame": name in contacts}

    meta = {
        "tag": tag, "subject": run["subject"], "task": run["task"],
        "run": run["run"], "anchor_uutc": anchor_us,
        "anchor_source": anchor["anchor_source"],
        "anchor_valid": bool(anchor["anchor_valid"]),
        "ieeg_fs": TARGET_FS, "pupil_fs": PUPIL_FS, "ripple_fs": RIPPLE_FS,
        "n_staged_channels": len(names),
        "n_analysis_contacts": sum(1 for n in names if n in contacts),
        "channels": [describe(n) for n in names],
        "ripple_channels": hippocampal,
        "eligible_seconds": summary["eligible_seconds"],
        "core_seconds": float(core.sum() / TARGET_FS),
        "pupil_quality": quality,
    }
    with open(stage_dir / "meta.json", "w") as handle:
        json.dump(meta, handle, indent=2)

    summary.update({
        "status": "ok", "n_staged_channels": len(names),
        "n_analysis_contacts": meta["n_analysis_contacts"],
        "n_ripple_channels": len(hippocampal),
        "core_seconds": meta["core_seconds"],
        "ieeg_samples": int(ieeg.shape[0]),
    })
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # No machine-local defaults: this file is published in the public mirror.
    # Set EBRAINS_DATASET and EBRAINS_STAGE, or pass the paths explicitly.
    parser.add_argument("--dataset", default=os.environ.get("EBRAINS_DATASET"),
                        help="BIDS root of the EBRAINS dataset "
                             "(default: $EBRAINS_DATASET)")
    parser.add_argument("--stage", default=os.environ.get("EBRAINS_STAGE"),
                        help="directory to write staged runs into "
                             "(default: $EBRAINS_STAGE)")
    parser.add_argument("--subjects", nargs="*",
                        default=["sub-003", "sub-004", "sub-005", "sub-006",
                                 "sub-007"])
    parser.add_argument("--tasks", nargs="*", default=list(PRIMARY_TASKS))
    parser.add_argument("--only", nargs="*", default=None,
                        help="restrict to these run tags, for validation runs")
    args = parser.parse_args()

    if not args.dataset or not args.stage:
        parser.error("--dataset and --stage are required "
                     "(or set EBRAINS_DATASET and EBRAINS_STAGE)")

    try:
        from pymef.mef_session import MefSession
    except ImportError:
        print("pymef is required: work/venv/bin/pip install pymef", file=sys.stderr)
        return 2

    dataset = Path(args.dataset)
    stage = Path(args.stage)
    stage.mkdir(parents=True, exist_ok=True)
    tasks = tuple(args.tasks)

    contact_sets, reference_sets, exclusions = {}, {}, {}
    for subject in args.subjects:
        keep, pool, dropped = eligible_contacts(dataset, subject, tasks)
        contact_sets[subject] = keep
        reference_sets[subject] = pool
        exclusions[subject] = dropped
        hippocampal = sum(1 for v in keep.values()
                          if v["region"] == HIPPOCAMPUS_LABEL)
        print(f"[ingest] {subject}: {len(keep)} analysis contacts "
              f"({hippocampal} hippocampal), {len(pool)} staged for the "
              f"Laplacian; dropped {dropped}", flush=True)

    runs = [r for r in discover_runs(dataset, args.subjects) if r["task"] in tasks]
    if args.only:
        runs = [r for r in runs
                if f"{r['subject']}_task-{r['task']}_run-{r['run']}" in args.only]
    print(f"[ingest] {len(runs)} runs for tasks {tasks}", flush=True)

    summaries = []
    for index, run in enumerate(runs, start=1):
        tag = f"{run['subject']}_task-{run['task']}_run-{run['run']}"
        print(f"[ingest] ({index}/{len(runs)}) {tag}", flush=True)
        try:
            summary = ingest_run(MefSession, run, contact_sets[run["subject"]],
                                 reference_sets[run["subject"]], stage, tasks)
        except Exception as error:  # noqa: BLE001 - recorded, never silently skipped
            summary = {"tag": tag, "status": f"error: {type(error).__name__}: {error}"}
        print(f"[ingest]     -> {summary.get('status')} "
              f"core={summary.get('core_seconds', float('nan')):.1f}s", flush=True)
        summaries.append(summary)

    fields = sorted({k for s in summaries for k in s})
    with open(stage / "ingest_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)
    with open(stage / "contact_exclusions.json", "w") as handle:
        json.dump(exclusions, handle, indent=2)

    ok = sum(1 for s in summaries if s.get("status") == "ok")
    print(f"[ingest] complete: {ok}/{len(summaries)} runs staged -> {stage}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
