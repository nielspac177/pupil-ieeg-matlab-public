# Analysis pipeline

The diagram below is the canonical description of how a recording becomes a
model estimate. It is the source of the pipeline figure in the manuscript: the
`.md` is what gets edited and version-controlled, and the image embedded in the
paper is rendered from it by `tools/render_mermaid.py`.

## 🔬 From recording to estimate

The single decision in this pipeline is the one that shapes every downstream
choice. A contact on which the peri-peak pupil trace never leaves the surrogate
band contributes an outcome of *exactly* zero rather than a small number, so the
outcome variable carries a point mass that no single regression can absorb. The
flow therefore splits at that decision and rejoins in a hurdle model whose first
part is fitted to every contact and whose second and third parts are fitted only
to contacts that cleared the band.

```mermaid
flowchart TB
    accTitle: Pupil iEEG Analysis Pipeline
    accDescr: Simultaneous intracranial and pupil recordings are reduced to a per-contact signed response, which is exactly zero when no suprathreshold excursion occurs, and the resulting semicontinuous outcome is then fitted with a three-part hurdle model.

    recordings([📥 Simultaneous sEEG and pupillometry])

    subgraph measurement ["🔬 Per-contact measurement"]
        peaks[⚡ Detect high-gamma peaks] --> window[📊 Average pupil around each peak]
        window --> surrogate[🔄 Compare with shifted-peak surrogate]
        surrogate --> longest[📝 Take longest contiguous run]
    end

    recordings --> measurement
    measurement --> cleared{🔍 Run clears the band?}

    cleared -->|No| zero[❌ Signed area is exactly zero]
    cleared -->|Yes| excursion[✅ Signed area and window fraction]

    prevalence[⚙️ Part 1 prevalence] --> direction[⚙️ Part 2 direction]
    direction --> magnitude[⚙️ Part 3 magnitude]

    zero --> prevalence
    excursion --> prevalence
    excursion --> direction
    magnitude --> estimates([📤 Region estimates with FDR control])

    classDef terminal fill:#ede9fe,stroke:#7c3aed,stroke-width:2px,color:#3b0764
    classDef process fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f
    classDef decision fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12
    classDef outcome fill:#f3f4f6,stroke:#6b7280,stroke-width:2px,color:#1f2937

    class recordings,estimates terminal
    class peaks,window,surrogate,longest,prevalence,direction,magnitude process
    class cleared decision
    class zero,excursion outcome
```

## 📐 What each part asks

Every model in the hurdle carries the same random-effects structure: contacts
nest within electrode shafts, and shafts nest within patients. Neighbouring
contacts on one shaft sit millimetres apart, and the patient is the unit of
biological replication, so treating contacts as independent would badly
overstate the evidence.

| Part | Fitted to | Question | Model |
| --- | --- | --- | --- |
| 1 Prevalence | all contacts | Does this contact couple at all? | binomial GLME |
| 2 Direction | excursion contacts | Dilation or constriction? | binomial GLME (**primary**) |
| 3 Magnitude | excursion contacts | How large is the response? | linear MME, asinh scale, Satterthwaite d.f. |

## ⚠️ The two properties that force this design

| Property | Consequence |
| --- | --- |
| The stored coupling statistic is a **contiguity criterion**, not a p-value — the fraction of the window spanned by the longest suprathreshold run | Thresholding it *selects* contacts; it does not identify significant ones, and no confirmatory claim may be conditioned on it |
| The signed response is **semicontinuous** — exactly zero wherever no run cleared the band | A Gaussian model on the pooled variable is misspecified; prevalence and direction must be modelled separately |

Both are stated in full in [`METHODS.md`](../METHODS.md) and
[`analysis_plan.md`](analysis_plan.md).
