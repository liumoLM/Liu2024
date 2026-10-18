# Solid findings on the ID1, ID2, and ID_N splits

This note summarizes the most strongly supported claims about the 89-type
subdivisions of COSMIC ID1, ID2, and ID_N that can be made from per-sample
attributions in Table S11 joined to Table S17 metadata
(`Table_S10_S17_joined_wide.csv`, 6,975 samples).

The supporting analysis script is `analyze_id1_id2_idN.R`. Outputs:
`id1_id2_idN_cooccurrence_focused.pdf`,
`id1_id2_idN_stripplot_by_stratum.pdf`,
`id1_id2_idN_summary.csv`,
`idN_alpha_vs_beta_investigation.csv`,
`idN_alpha_vs_beta_panels.pdf`.

> **Note (2026-10-18).** The numbers below use the paper's definitions
> (MSI-H by MSIseq, hypermutated meaning > 4,000 indels in the R <= 9
> clipped spectra). They are from the outputs of `analyze_id1_id2_idN.R`
> regenerated on 2026-10-17. An earlier version of this note used the
> cohort `MSI_status` and Table S11 column sums >= 5,000, which gave 91
> MSI-H and 476 MSS hypermutated tumors. Several findings changed, notably
> items 2, 6, 11, and 12 below.

## Definitions

- **Signature share** for sample s = (indels attributed to that signature
  in s) / (total indels attributed to all signatures in s).
- A sample is "present" for a signature when share >= 0.05 (5%). This is
  a relative criterion. It is more robust to attribution noise than
  `exposure > 0`, but it causes apparent depletion of clock-like
  signatures in MSI-H samples because the MMR block consumes most of
  those samples' indels. See "Attribution displacement" below.
- **Strata used throughout**: MSS non-hypermutated, not MSI-H by MSIseq
  and <= 4,000 indels (n = 6,725). MSS hypermutated, not MSI-H by
  MSIseq and > 4,000 indels (n = 48). MSI-H, MSI-H by MSIseq (n = 202).
  Indels are counted from the R <= 9 clipped spectra.

## A. ID1 split: InsDel1a (clock-like), InsDel1b (POLE), InsDel1c (liver)

1. **InsDel1a is broadly active and clock-like.** Present at share
   >= 5% in 56.2% of MSS non-hypermutated tumors across essentially
   all cancer types. Mean share 0.21 in MSS non-hyper. Correlates
   with patient age in multiple tissues (Demographic Associations
   section of the manuscript).
2. **InsDel1b is not clearly concentrated in MSS hypermutated tumors.**
   Present at share >= 5% in 1 of 48 MSS hypermutated tumors (2.1%),
   2.1% of MSS non-hypermutated tumors, and 5.4% of MSI-H tumors. Mean
   share is 0.014 in MSS hypermutated, 0.006 in MSS non-hyper, and
   0.027 in MSI-H. With the old strata, InsDel1b looked concentrated in
   MSS hypermutated tumors. Several InsDel1b-high tumors carry mutations
   in the proofreading domain of polymerase epsilon (Table S9),
   consistent with a POLE-proofreading-defect etiology.
3. **InsDel1c is concentrated in MSS non-hypermutated liver
   cancers** and is present in 8.2% of MSS non-hyper tumors overall.
4. **InsDel1a / 1b / 1c are nearly mutually exclusive at the sample
   level.** Co-occurrence phi for InsDel1a vs InsDel1c = -0.33,
   q = 3e-209, and zero samples reach >= 5% share for both. InsDel1a vs
   InsDel1b: phi = -0.16, q = 8e-49, and 2 samples reach 5% for both.
   If the three components were one underlying process partitioned by
   noise, the attributor would distribute them jointly across samples
   instead of partitioning samples cleanly. This is the strongest data
   argument that the ID1 subdivision is biologically real.

## B. ID2 split: InsDel2a (broadly active), InsDel2b / 2c (MSI-restricted)

5. **InsDel2b and InsDel2c are essentially restricted to MSI-H
   tumors.** Present at share >= 5% in 26.7% (InsDel2b) and 22.8%
   (InsDel2c) of MSI-H tumors, versus 0.6% and 0.7% of MSS
   non-hypermutated tumors and none of the MSS hypermutated tumors.
   Mean shares in MSI-H are 0.049 and 0.045.
6. **InsDel2a is enriched in MSI-H but not restricted to it.** Present
   at share >= 5% in 29.2% of MSI-H, 10.4% of MSS hypermutated, and
   2.6% of MSS non-hyper tumors.
7. **The MMR block co-occurs coherently.** C_ID7, InsDel2b, InsDel2c,
   InsDel_D, and InsDel_J co-occur at high phi (top pair C_ID7 vs
   InsDel_J phi = 0.55, q = 2e-97), consistent with a shared
   MMR-deficiency etiology. The 83-type system collapses InsDel2a
   with the MMR-driven InsDel2b/2c and so cannot resolve the
   MMR-specific contribution.

## C. ID_N split: InsDel_Nβ (broadly active) and InsDel_Nα (rare)

8. **InsDel_Nβ is broadly distributed and one of the most commonly
   detected indel signatures in the cohort.** Present at share
   >= 5% in 39.0% of MSS non-hypermutated and 20.8% of MSS
   hypermutated tumors. Its behavior is consistent with a
   broadly-active, near-clock-like process.
9. **InsDel_Nα is rare**, present at share >= 5% in only 60 of
   6,975 tumors (0.9%). Only the 89- and 476-type schemes
   separate it from InsDel_Nβ.
10. **InsDel_Nα and InsDel_Nβ are mutually exclusive at the sample
    level.** Zero samples reach the 5% threshold for both. Evidence
    that the N split is real.
11. **InsDel_Nβ is rare in MSI-H tumors at the 5% threshold** (10 of
    202, 5.0%). With the old strata it was 0 of 91. This is probably an
    attribution-displacement effect: in MSI-H tumors the MMR block
    consumes most of the indel budget, so InsDel_Nβ falls below 5%
    share even though it still contributes some mutations in absolute
    terms.

## D. Cohort stratification structure

12. **Different 89-type components dominate each of the three strata**:
    - **MSI-H** is dominated by C_ID7 (73.3% present), InsDel_D
      (37.1%), InsDel_J (35.6%), InsDel2a (29.2%), InsDel2b (26.7%),
      and InsDel2c (22.8%).
    - **MSS hypermutated** (48 tumors) is too small for a clear
      pattern. InsDel_Nβ (20.8%) and InsDel2a (10.4%) are the most
      common of the focused signatures. With the old strata (476
      tumors), this group seemed to surface InsDel1b and an elevated
      InsDel_Nβ, which no longer holds.
    - **MSS non-hypermutated** is where clock-like (InsDel1a,
      InsDel5a, InsDel_B, InsDel_Nβ), liver-associated (InsDel1c,
      InsDel_B), and tissue-restricted exogenous (InsDel_H in
      skin, InsDel3a in bladder, C_ID18 in colon / esophagus)
      signatures are visible without competition from high-burden
      samples.
    - This is the narrative behind the 4-panel Figure 4
      reproduction in `Figure4_dotplot_panels.pdf`.

## E. Attribution displacement caveat

13. **The 5% relative-share threshold causes apparent depletion of
    clock-like signatures in MSI-H tumors.** InsDel1a's median
    share drops from 0.164 in MSS non-hypermutated to 0 in MSI-H
    (present in 15.8% of MSI-H), even though InsDel1a still contributes
    indels in MSI-H samples in absolute terms. InsDel_Nβ shows the same
    pattern. Clock-like analyses should therefore be restricted to MSS.

## What the data does NOT support

- **No "POLE- and MMR-restricted InsDel_Nα" claim.** Only 4 of 202
  MSI-H tumors and 1 of 48 MSS hypermutated tumors reach 5% share for
  InsDel_Nα. The Spearman correlation between InsDel_Nα share and
  C_ID7 share across the cohort is 0.016.
- **No "ID2 mutations reassigned to InsDel_Nβ" claim at the
  per-sample level.** Spearman correlation between InsDel_Nβ share
  and the sum of InsDel2a + InsDel2b + InsDel2c share is -0.029
  across all samples. Any reassignment that may have happened
  between 83-type ID2 and 89-type InsDel_Nβ would need to be
  demonstrated at the mutation-class level, not the sample level.
- **No InsDel_Nα biological etiology claim** beyond "rare and
  resolved only by 89/476-type schemes". The biology of InsDel_Nα
  is not yet established.

## Recommended visuals for the manuscript

- **`Figure4_dotplot_panels.pdf`**: 4-panel dot plot (All / MSI-H /
  MSS hypermutated / MSS non-hypermutated). Single best figure
  for the stratification story.
- **`id1_id2_idN_cooccurrence_focused.pdf`**: 11x11 phi heatmap
  restricted to ID1, ID2, and ID_N children plus the MMR block.
  Makes the mutual exclusion of InsDel1a/1b/1c and InsDel_Nα/Nβ,
  and the coherence of the MMR block, visible in one panel.
- **`id1_id2_idN_stripplot_by_stratum.pdf`**: per-sample share by
  stratum for each focused signature. Makes the
  MSS-non-hyper / MSS-hyper / MSI-H stratification visible at the
  single-sample level rather than as cancer-type averages.
