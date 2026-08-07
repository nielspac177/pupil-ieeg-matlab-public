# pupil-ieeg-matlab (public mirror)

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
