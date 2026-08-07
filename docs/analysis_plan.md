# Analysis contract

## Primary research question

Does the *direction* of pupil-coupled high-frequency activity differ by
anatomical region in the human brain? Specifically, does hippocampal high-gamma
activity accompany pupil constriction while neocortical high-gamma activity
accompanies dilation — a pattern incompatible with a single scalar arousal
signal broadcast uniformly to cortex and hippocampus?

## Two facts about the derived data that determine the design

1. **`RespSig` is not a p-value.** It is the fraction of the +/-5 s peri-peak
   window covered by the longest contiguous run in which the pupil response
   confidence band clears the surrogate noise band (`PupilHG.m:517`, `:533`).
   It has no calibrated null distribution. `RespSig > 0.10` selects contacts; it
   does not identify significant ones. Any table conditioned on it is named and
   labelled descriptive, and no confirmatory claim rests on it.

2. **`RespAreaNet` is semicontinuous.** It is exactly zero for 628 of 913
   contacts, precisely those with `RespSig == 0`. A Gaussian model on the pooled
   variable is misspecified.

## Confirmatory hierarchy

The analysis is a two-part hurdle decomposition plus a magnitude model.

1. **Prevalence.** Binomial GLME over all 913 contacts: does a contact show any
   suprathreshold excursion? Region fixed effect, random intercepts for patient
   and shaft-within-patient. Threshold-free.
2. **Direction (primary).** Binomial GLME over the 285 excursion contacts:
   given an excursion, is it dilation- or constriction-linked? The prespecified
   single-degree-of-freedom contrast is hippocampal versus extrahippocampal.
3. **Magnitude (secondary).** Linear MME on an asinh scale over excursion
   contacts, with Satterthwaite denominator degrees of freedom.

## Statistical safeguards

- The patient is the biological replication unit; contacts nest in shafts and
  shafts in patients, and both levels carry random intercepts.
- Linear-model denominator degrees of freedom use Satterthwaite, never the
  residual default, which would put the contact count in the denominator.
- Region contrast families are FDR-controlled (Benjamini-Hochberg) and reported
  as q values.
- Report estimates and confidence intervals, not only p-values.
- The primary result carries two independent robustness checks: a within-patient
  paired comparison using the patient as the unit of analysis (Wilcoxon signed
  rank and an exact sign test), and re-estimation across a range of historical
  selection thresholds.
- Regions with complete separation in the direction model are flagged, not
  reported as effects.
- Confirmatory, secondary, and exploratory outputs live in separately named
  tables.

## Explicitly demoted claims

The following were presented as findings in earlier drafts and are now either
exploratory or withdrawn:

- **Middle temporal gyrus as the principal hub.** MTG leads on raw counts but is
  not more likely than the rest of the implant to show coupling once coverage is
  modelled (q = 0.20). Its prominence was an artefact of electrode placement.
- **Region-specific "spectral fingerprints" and the aperiodic-exponent link.**
  The quantity fitted is a coupling-amplitude spectrum, not a neural power
  spectrum, and the between-region separation is not robust to the estimator.
  Exploratory only.
- **A temporally ordered network cascade.** Stored fit lags describe a
  peripheral latency within each region; they do not order regional engagement.
- **The thalamus result.** Four contacts in one patient. Not reported.

## Stopping rule

Stop adding primary analyses once prevalence, direction, and the two direction
robustness checks are resolved. New post-hoc analyses must be labelled
exploratory and added to the deviation log.
