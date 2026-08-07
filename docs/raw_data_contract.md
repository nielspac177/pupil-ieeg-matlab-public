# Raw-data contract

The supplied project archive has been inventoried and contains derived MATLAB
tables, but not the synchronized raw iEEG, pupil streams, channel maps, anatomy,
SOZ/IED annotations, or luminance and gaze traces. Accordingly, every session is marked `n/a` in
`config/raw_manifest.tsv`; this is a documented availability conclusion, not a
request for the user to locate additional paths. The raw modules remain in the
repository as prospective, tested code should an independently authorized
source archive ever be recovered.

## Required for core signal analysis

- continuous iEEG (`.ns2`, `.nf3`, or previously verified lossless conversion)
- synchronized pupil samples and validity/interpolation mask
- channel map with native and MNI coordinates
- bad-channel annotations

## Required for confirmatory ripple analysis

- sampling frequency and analog-filter metadata
- anatomical hippocampal contact labels
- adjacent-contact or white-matter reference information
- IED, seizure, and seizure-onset-zone annotations
- enough continuous valid pupil exposure to estimate quantile-state durations

## Required for inter-session pupil prediction

- at least two comparable sessions per participant
- trial/event timing and task identity
- session-wise pupil and iEEG QC
- splits defined by entire blocks/sessions, never overlapping random windows

## Required for BIDS release

- de-identified accepted-format iEEG
- synchronized eye-tracking and events
- electrode and channel metadata
- coordinate-system provenance
- authorized, defaced anatomy when shared
- completed governance and independent PHI review

With the supplied archive, a raw iEEG BIDS release cannot be constructed or
validated. A separately labeled derived-data package containing de-identified
tables and MNI coordinates may be considered only after governance review; it
must not be presented as a raw BIDS iEEG dataset.
