# Questions for Mo: indel counts, poly-T ratios, and hypermutation

Compiled 2026-10-15. All counts below are over the 6,975 tumors in Table S17.
"Spectrum total" means the number of indels in a tumor's catalog in
`unified_indels_data` (`liu_et_al_83_spectra.tsv`). The 83-, 89-, and
476-type catalogs give the same total for every tumor. These catalogs count
repeat lengths only up to 9 (R <= 9).

The main question behind everything below: for each threshold, did we use
indel counts from before or after capping repeat lengths at R <= 9? Either
could be defended, but we need to know which one we used, say so in the
methods, and make the numbers reproducible. If the pre-cap counts were used,
we may have to recompute them from the VCFs.

## 1. Table S17 `mutation_burden`

- The legend calls it "Number of somatic indels in the tumor".
- It equals the spectrum total for only 486 tumors (485 PCAWG, 1 Hartwig).
  For all the others it is larger, by a median of 6% and at most twice the
  spectrum total.
- Guess: these are counts from before capping at R <= 9. This is not
  confirmed.

Questions:

- What data was `mutation_burden` computed from, and with what code?
- Why does it agree with the current catalogs for most PCAWG tumors but
  almost no Hartwig tumors?

## 2. Table S17 poly-T ratio

- The column is "Ratio of 1bp T insertion/deletion at >=5 polyT tracts". The
  legend defines it as the fraction of the tumor's indels that are 1-bp T
  insertions or deletions in poly-T tracts of 5 or more Ts.
- Computed from the current catalogs as (INS:T:1:5+ plus DEL:T:1:5+) divided
  by the spectrum total, it correlates at 0.998 with S17 but matches exactly
  for only 1,249 tumors.
- Dividing by S17's own `mutation_burden` matches only 815 tumors, so the
  S17 ratio was not simply computed from the S17 burden.

Questions:

- What numerator and denominator were used, and from which version of the
  data?
- Is there code for it?

## 3. Routing between SigProfilerAssignment and mSigAct

- The signature-assignment methods in `ms.qmd` say mSigAct was
  used for tumors with 5,000 or more indels, or with a poly-T ratio of 0.5 or
  more. All other tumors went to SigProfilerAssignment on the 172-type
  catalog.
- With S17's values, 43 tumors fall on different sides of the poly-T cutoff
  than with the current catalogs. (`paper_to_do.md` item 6 says 44 for the
  combined routing rule.)
- With S17's `mutation_burden`, 221 tumors have 5,000 or more indels. With
  the current catalogs, 195 do. So 26 tumors differ, all of them MSS.
- The Table S11 column sums show which tumors were apparently fit with
  mSigAct. SigProfilerAssignment fits the 172-type catalog, which stacks the
  83- and 89-type catalogs and counts each indel twice. So its exposures sum
  to twice the spectrum total (6,472 tumors). For the other 503 the sum is
  about the spectrum total (see the note in `ms.qmd` on Table S11 totals).
- These 503 tumors do not match the stated routing rule, whichever counts
  are used. With S17's values, 51 of the tumors with single-counted sums do
  not meet the mSigAct rule, and 1,153 tumors that do meet it have
  double-counted sums. With the current catalogs, the numbers are 46 and
  1,188.

Questions:

- Exactly which tumors went to mSigAct, and what rule selected them?
- Were the routing thresholds applied to pre-cap or post-cap counts?
- Is there a list or code for the routing?

## 4. What defines a hypermutated tumor

Different parts of the paper use different definitions:

| Where | Rule | Count used |
|---|---|---|
| Methods, assignment routing (`ms.qmd`) | 5,000 or more indels | not stated |
| Methods, MSI calling (`ms.qmd`) | more than 4,000 indels, then MSIseq run on these | not stated |
| Figure 4A, MSS-hyper column (`build_fig4/all_cancer_types_enrichment.R`) | not MSI-H by MSIseq, and more than 4,000 | sum of signature exposures, which double-counts most tumors |
| `explore-attributions/reproduce_fig5_panelC_with_MSS_hyper.R` (exploratory, older figure numbering) | MSS and more than 5,000 | Table S11 column sums, which double-count most tumors |

Effect on the Figure 4A MSS-hyper group: by the exposure sums, 723 tumors
have more than 4,000 "indels". By the spectrum totals, only 48 do. For
tumors fit with SigProfilerAssignment, the effective cutoff is about 2,000
indels.

The MSS group for the exploratory script also changes: MSS tumors with more
than 5,000 indels number 475 by Table S11 sums, 130 by S17
`mutation_burden`, and 104 by the spectrum totals.

Questions:

- Which tumors were labeled high-TMB for the MSIseq step, with which counts?
- Was the Figure 4A threshold of 4,000 meant to apply to real indel counts?
  If so, the figure should be recomputed with the spectrum totals (or the
  exposure sums halved for the SigProfilerAssignment tumors).
- Should the paper use a single definition of hypermutated? If so, which
  threshold, and pre-cap or post-cap counts?

## 5. MSI status

`sample_info.tsv` has two MSI calls that disagree. `MSIseq_MSI.H` is TRUE for
202 tumors, and `MSI_status` is "MSI" for 91 (a subset of the 202). Figure 4A
uses MSIseq. Table S9 reports MSIseq and gives `MSI_status` alongside it.

Question:

- Which call should the paper treat as the MSI status, and is that stated
  consistently in the methods and legends?
