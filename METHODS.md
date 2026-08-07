# Methods

Exported verbatim from the submitted manuscript by `tools/build_public_repo.py`. Edit the manuscript builder in the private repository, not this file.

![Measurement and model](figures/Fig0_methods_schematic.png)

*Figure 1. Measurement and model. (a) High-gamma peaks are detected in each contact. (b) Pupil diameter is averaged around every peak and compared with a shifted-peak surrogate band; the stored statistic is the fraction of the window spanned by the longest contiguous run above that band, which is a contiguity criterion and not a p-value. (c) A contact with no suprathreshold run contributes an exact zero, so the outcome is semicontinuous and a single regression on it is misspecified; the distribution shown is schematic and the observed proportions are given in the Results. (d) The two processes are therefore modelled separately, with contacts nested in shafts nested in patients. Traces in (a) and (b) are synthetic illustrations, not data.*

## Materials and methods

### Participants

18 adults undergoing stereo-electroencephalography for the localisation of medically refractory epilepsy were recruited between 2019 and 2021. All participants gave informed written consent under a protocol approved by the University of Utah Institutional Review Board (#00069440). Twelve participants performed a visuospatial working-memory task, three a thermal-perception task, and three were recorded during eyes-open quiet wakefulness. Electrode placement was determined solely on clinical grounds.

### Electrodes, recording and localisation

Clinical stereo-electrodes and subdural arrays (4–10 platinum contacts per shaft; Ad-Tech, Racine, WI, USA, or DIXI Medical, Besançon, France) were localised with LeGUI (Davis et al., 2021), which coregisters the postoperative computed-tomography scan to the preoperative T1 MRI, segments contacts, and assigns each to a Neuromorphometrics parcel in native and MNI space. Across participants, 913 contacts on 146 distinct shafts passed artefact criteria and were distributed over 47 anatomical regions. Field potentials were digitised at 1 kHz on a 128-channel NeuroPort system (Blackrock Microsystems, Salt Lake City, UT, USA), bandpass-filtered 0.3–250 Hz, and re-referenced with a Laplacian scheme relative to the nearest two same-shaft contacts. Pupil diameter was recorded at 250 Hz with a Pupil Core infrared eye-tracker (Pupil Labs, Berlin, Germany) and time-locked to the neurophysiological acquisition. Template surfaces for display were generated with LeGUI's own LeG_genSurfaces routine from the bundled tissue probability maps.

### High-gamma peaks and the peri-peak pupil response

Continuous data were notch-filtered at 60 Hz and decomposed into 71 wavelet-based frequency bands between 2 and 256 Hz; the primary band was high-gamma (70–170 Hz). High-gamma power peaks were defined as local maxima exceeding the median plus five interquartile ranges. For each contact we extracted mean pupil diameter in a ±20 s window around every peak, demeaned to the −20 to −10 s baseline, and built a surrogate distribution by randomly shifting peak times (100 iterations).

Within the central ±5 s of the window we identified time points at which the response confidence band lay entirely outside the root-mean-square of the surrogate standard error, separately in the positive and negative directions, and retained the longest contiguous run of such points. Two quantities were stored: the signed area of that run, and the fraction of the ±5 s window it spanned.

This second quantity requires an explicit statement, because it has been described inaccurately in earlier drafts of this work. It is a contiguity statistic: the proportion of the analysis window covered by the longest suprathreshold run. It is not a permutation p-value, it has no calibrated null distribution, and a threshold applied to it does not control any error rate. Throughout this paper, a contact passing the historical criterion (contiguity > 0.10; n = 177) is therefore called a ‘selected’ contact and is used only descriptively. No confirmatory claim in this paper is conditioned on it.

### Statistical analysis

The stored signed response is exactly zero for every contact on which no suprathreshold excursion occurred — 628 of 913 contacts (69%) in this dataset. The variable is therefore semicontinuous, and a Gaussian model fitted to the pooled variable is misspecified. We used a two-part hurdle decomposition.

Part 1 (prevalence) asked whether a contact showed any excursion, using a binomial generalized linear mixed model with region as a fixed effect and random intercepts for patient and for shaft nested within patient, fitted by Laplace approximation over all contacts. Part 2 (direction), the primary confirmatory analysis, asked whether an excursion was dilation- or constriction-linked, using the same random-effects structure over the contacts that showed an excursion. The prespecified fixed effect was hippocampal versus extrahippocampal. Part 3 (magnitude) modelled the signed response among excursion contacts with a linear mixed model on an inverse-hyperbolic-sine scale, with Satterthwaite denominator degrees of freedom so that inference is not carried out at the contact level.

Contacts are not independent observations: neighbouring contacts on one shaft are millimetres apart and share a reference neighbourhood, and the patient is the biological unit of replication. Both levels therefore carry random intercepts. Families of region contrasts were corrected with the Benjamini–Hochberg procedure and are reported as q values. Two additional checks are reported for the primary result: a within-patient paired comparison using the patient as the unit of analysis, and a sensitivity analysis re-estimating the effect across a range of historical selection thresholds. Analyses were performed in MATLAB R2024b (Statistics and Machine Learning Toolbox); code and derived tables are available as described under Data availability.

Figure 2. Analysis pipeline. Each contact yields one signed peri-peak response. The single branch is the one that determines the model: a contact whose longest contiguous run never clears the surrogate band contributes an outcome of exactly zero rather than a small number, so prevalence must be modelled separately from direction and magnitude. Part 1 is fitted to all 913 contacts; parts 2 and 3 are fitted to the 285 with an excursion. All three carry random intercepts for patient and for electrode shaft within patient. This diagram is generated from docs/methods_diagram.md in the project repository, where it is version-controlled as text.
