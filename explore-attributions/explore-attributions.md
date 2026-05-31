# Exploring the indel-signature attributions in Table S10

This folder contains exploratory analyses of the per-sample 83-type / 89-type
indel signature attributions (`Sup Tables/Table S10 ...tsv`) joined to the
sample metadata (`Sup Tables/Table S17 ...xlsx`). The goal is to assess how
robust the attributions are and to surface hints about the biology behind
each signature.

All scripts in this folder read inputs from `../Sup Tables/` and write their
outputs alongside themselves.

## Data inputs

- `Sup Tables/Table S10 83-type and 89-type signature assignment.tsv`. Per
  sample (column) mutation counts attributed to each joint 83/89-type
  signature (row). Row labels have the form `C_ID1/InsDel1a`. The 89-type
  exposures equal the corresponding rows of
  `liu_2025_draft_release/Manuscript_data/finalized_cap9/liu_et_al_89_assignment.tsv`
  exactly, and the 83-type exposures equal the sums of joint rows that share
  the same 83-type prefix (verified across all 6975 samples, zero
  discrepancies).
- `Sup Tables/Table S17 metadata of 6975 samples.xlsx`. Sample metadata with
  `Patient`, `Major Cancer Type`, `Cancer Type`, `cohort`, `MSI_status`
  ("MSS", "MSI-H", or a few "MSI"), `mutation_burden`, `MSIseq`, and `Ratio
  of 1bp T insertion/deletion at >=5 polyT tracts` (per-sample fraction of
  indels that are 1-bp T at long poly-T tracts; ranges 0 to 0.85, median
  0.30, used as the routing criterion to mSigAct when >= 0.5).

## Scripts and outputs

| Script | Output |
|---|---|
| `build_assignment_csv.R` | `Table_S10_S17_joined_wide.csv`: one row per sample, columns `cancer.type`, `sample`, then the 42 signature exposures. |
| `plot_signature_assignments.R` | `Table_S10_hamburger_plots_mf0.pdf`: per-signature "hamburger" plot of exposures across cancer types (one page per signature, each dot is a sample, red dashed line is per-cancer-type median, MSI-H samples are red triangles). `--min-fraction` filters samples whose signature share is below the cutoff. |
| `analyze_signatures.R` | `analysis1_msi_vs_mss.csv`, `analysis2_cooccurrence_global.csv`, `analysis2_cooccurrence_global_phi.pdf`, `analysis2_cooccurrence_phi_by_tissue.pdf` (see below). |
| `plot_exposure_dotplot.R` | `Figure4_dotplot_panels.pdf`: 4-page reproduction of the Figure-4 dot plot (size = active proportion, color = median mut/Mb), shown for All / MSI / MSS-hypermutated / MSS-non-hypermutated subsets. |

A sample is "exposed" to a signature, throughout the analyses below, when
that signature accounts for >= 5% of the sample's total assigned indels.
This relative threshold is more robust to attribution noise than `exposure
> 0` (which counts isolated mutations).

## Analysis 1: MSI vs MSS within each Major Cancer Type

Fisher's exact test on a 2x2 table of (MSI, MSS) x (signature present,
absent), computed separately for each signature within each Major Cancer
Type. BH-adjusted across the full set of tests. Output:
`analysis1_msi_vs_mss.csv`.

**Top positive associations (signature enriched in MSI vs MSS).** These
confirm the expected MMR / POLE block:

- `C_ID7/InsDel7` is the strongest MSI marker: present in 100% of MSI-H
  samples in Colon, Prostate, Uterus, Breast, Stomach, CNS, and Other, vs 0
  to 3% of MSS samples in the same tissues. q-values down to 1e-96.
- `C_ID2/InsDel2b`, `C_ID2/InsDel2c`, `ID_D/InsDel_D`, `ID_J/InsDel_J`,
  `ID_M/InsDel_M`, and `ID_N/InsDel_N_alpha` follow the same pattern.
- The Methods note that `C_ID2` deletions of single T from poly-T tracts of
  length >= 6 dominate MSI spectra is supported by the data, and the
  89-type split into 2a (broadly active) and 2b/2c (MSI-restricted) is
  recovered cleanly.

**Top negative associations (signature *depleted* in MSI vs MSS).**

- `C_ID1/InsDel1a` is dramatically depleted in MSI samples of Colon (4%
  MSI-positive vs 59% MSS-positive) and Uterus (0% vs 60%). This is almost
  certainly an attribution-displacement artifact, not biology: when MSI
  drives the spectrum, the residual signal available for slow-clock
  attribution falls below the 5% threshold.
- `ID_N/InsDel_N_beta` is depleted in MSI Colon for the same reason.

The displacement signal is a warning that the binary "presence" of
clock-like signatures should not be interpreted in MSI-H samples; restrict
clock-like analyses to MSS.

## Analysis 2: pairwise signature co-occurrence

For every pair of the 42 joint signatures, compute the 2x2 co-presence
table, phi coefficient, Haldane-corrected log odds ratio, Fisher's p, and
BH-adjusted p. Output: `analysis2_cooccurrence_global.csv`. Heatmaps:
`analysis2_cooccurrence_global_phi.pdf` (global) and
`analysis2_cooccurrence_phi_by_tissue.pdf` (per Major Cancer Type).

**MMR block (strong positive co-occurrence).** `C_ID7/InsDel7`,
`ID_J/InsDel_J`, `ID_D/InsDel_D`, `C_ID2/InsDel2b`, and `C_ID2/InsDel2c`
all co-occur at high phi (top pair: C_ID7-ID_J, phi = 0.55, q = 4e-96).
This corroborates analysis 1, MMR drives a coherent block of indel
signatures, and supports a shared MMR-deficiency etiology for the
89-type-novel members of the block (ID_D, ID_J).

**UV pair.** `C_ID13/InsDel13` and `ID_H/InsDel_H` co-occur (phi = 0.19,
q = 5e-21). This supports a shared UV etiology for InsDel_H and is
consistent with InsDel_H being skin-restricted (see hamburger plot for
ID_H/InsDel_H).

**NHEJ / DSB pair.** `C_ID6/InsDel6` and `C_ID8/InsDel8` co-occur (phi =
0.18, q = 1e-32), consistent with the known mechanistic relationship
between HR-defect deletions (ID6) and NHEJ-mediated DSB repair (ID8).

**GI-ROS adjacency.** `C_ID14/InsDel14` and `ID_N/InsDel_N_beta` co-occur
strongly (phi = 0.24, q = 2e-85, n_present_both = 288). This pair has
enough samples to support a real biological link, and is consistent with
the GI-ROS cluster described in the manuscript Figure 5 narrative.

**The 89-type split is real.** The largest negative association in the
whole table is `C_ID1/InsDel1a` vs `C_ID1/InsDel1c` (phi = -0.33, q = 4e-208,
zero samples co-present). `C_ID1/InsDel1a` vs `C_ID1/InsDel1b` is the same
pattern. If the 89-type split of COSMIC ID1 were just a noisy
re-partitioning of one underlying process, these subcomponents would
co-occur, instead, the attributor consistently assigns each sample to a
single component. The same logic supports the split of `C_ID5` (InsDel5a
strongly anti-correlated with `ID_B/InsDel_B`, phi = -0.18, q = 4e-72).
**This is the strongest argument in the data for the biological reality of
the 89-type subdivisions.**

**Attribution-displacement signal.** Many negative phi values involve
`C_ID1/InsDel1a` paired with MSI signatures (InsDel2a, C_ID7, InsDel_N_beta,
InsDel3b). The displacement effect from analysis 1 is visible here as
well: clock-like InsDel1a "leaves" the attribution as MSI mutations
dominate. Anti-correlation is therefore an artifact of attribution, not
biological mutual exclusion.

## Figure-4 reproduction with subset panels

`Figure4_dotplot_panels.pdf` reproduces the original Figure 4a/4b dot plot
(rows = signatures interleaving 83-type aggregates with their 89-type
children, columns = Major Cancer Type, size = active proportion, color =
median mutations / Mb on a log scale), in four panels:

1. **All samples (n = 6975).** Matches the published Figure 4 layout.
2. **MSI samples only (n = 91).** Isolates the MMR-driven exposures so
   they're not diluted by the MSS background.
3. **MSS hypermutated samples, total indels >= 5000 (n = 476).** Surfaces
   POLE-proofreading and other non-MSI high-burden processes that would
   otherwise blend with MSI.
4. **MSS non-hypermutated samples (n = 6408).** The "background" panel,
   where clock-like, tissue-specific exogenous, and rare-but-real
   signatures are visible without competition from high-burden samples.

The hypermutation threshold (5000 total indels) matches the high-TMB
cutoff used in the Methods to route samples to mSigAct.

## What this exploration suggests for the manuscript

1. **Confirmed-biology block.** The MSI/MMR signatures (C_ID7, C_ID2b/c,
   ID_D, ID_J, plus their 89-type analogues) are robust on both
   analyses, they are MSI-enriched within multiple tissues *and* co-occur
   coherently across the cohort. These are safe to interpret causally.
2. **Real 89-type subdivisions.** The mutual-exclusion pattern of
   InsDel1a/1b/1c and of InsDel5a/InsDel_B is direct evidence that the
   89-type split captures genuinely distinct processes, not noise.
3. **UV etiology for InsDel_H** is supported by both its skin-restricted
   tissue distribution (hamburger plot) and its co-occurrence with the
   COSMIC UV indel ID13 (analysis 2). This is a defensible biological
   claim.
4. **Clock-like signatures should be analyzed within MSS only.** The
   displacement artifact (InsDel1a depleted in MSI tumors) cannot be
   distinguished from real biology in mixed-MSI analyses.
5. **Rare signatures (InsDel_C, InsDel_I, InsDel_K, InsDel_L, InsDel_M)
   need verified counts before tissue claims.** The to-do list values
   ("9 tumors only pancreas/prostate/colon" etc.) didn't match the actual
   joined data. The corrected counts are visible in the hamburger plots
   and the dot-plot panels.
6. **Some signatures may not be robust enough to interpret.** Anything
   appearing predominantly in one cohort, or whose only enrichment is
   borderline q-value in a single tissue, should be down-weighted.
   Analysis 1 (`analysis1_msi_vs_mss.csv`) and a planned cohort-effect
   analysis (not yet run) can be used to triage these.
