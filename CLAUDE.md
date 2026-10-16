# Local notes for the 0Liu2024 repo

## Key data sources for sample-level analyses

- **89-type / 83-type indel signature exposures per sample**:
  `Sup Tables/Table S11 83-type and 89-type signature assignment.tsv`
  - Rows are signatures named `C_ID*/InsDel*` (83-type/89-type pair).
  - Columns are 6975 samples (Patient IDs).
  - Cell values are exposures (indel counts). Column sums are NOT total
    indels per sample. For the 6,472 samples fit with SigProfilerAssignment
    to the 172-type catalog (83-type and 89-type stacked), each indel is
    counted twice, so the sum is twice the number of indels. For the other
    503 (presumably fit with mSigAct), the sum is about the number of
    indels. For indel counts, use the column sums of
    `unified_indels_data/liu_et_al_83_spectra.tsv`.

- **Sample metadata (cancer type, MSI status, etc.)**:
  `Sup Tables/Table S17 metadata of 6975 samples.xlsx`
  - Key columns: `Patient`, `Major Cancer Type`, `MSI_status` (values `MSI` / `MSS`),
    `mutation_burden`, `cohort`, `Cancer Type`.
  - All 6975 Patient IDs match the column names in Table S11.

## Sample grouping convention for the MSS-hyper enrichment analysis

The MSS-hyper enrichment analysis is now Figure 4A. The name "Figure 5
panel C" below and in `explore-attributions/` is an older figure number.
The current Figure 5C is genic vs. intergenic enrichment. The published
Figure 4A is built by `build_fig4/all_cancer_types_enrichment.R` in the code
repo with a different rule: MSI-H by MSIseq, then MSS-hyper if the sum of
signature exposures is > 4000. Because those sums double-count most samples,
that is in effect > 2000 indels for most of them. Which threshold and which
counts are intended is an open question (see `QUESTIONS_FOR_MO.md`).

The exploratory rule, used in the reference script below: when stratifying
samples into mutually exclusive groups for the Fisher enrichment analysis,
assign in this priority order so each sample lands in exactly one column:

1. **MSI-H**: `MSI_status == "MSI"`.
2. **MSS-hyper**: `MSS_status == "MSS"` AND total indels (column sum in
   Table S11) > 5000. Table S11 column sums double-count most samples, so
   this is in effect > 2500 indels for most of them.
3. **Major cancer type** (`Bladder`, `Colon`, `Esophagus`, `Kidney`,
   `Liver`, `Lung`, `Ovary`, `Prostate`, `Skin`): MSS, total indels <= 5000,
   and `Major Cancer Type` matches.

Samples in other cancer types (Breast, Pancreas, CNS, etc.) are not
displayed as their own column but still contribute to the "remaining"
pool of each Fisher's exact test.

Reference script: `explore-attributions/reproduce_fig5_panelC_with_MSS_hyper.R`.

## Source of truth for numbers in ms.qmd

All numbers quoted in `ms.qmd` (signature counts, proportions, sample counts,
etc.) must come from the code and data in `~/github/Liu2026_code_and_data/`,
not from memory or from earlier drafts. Locate the script or input file for
the relevant figure or table there (for example `build_fig5/` for Figure 5,
`ID83_topography_analysis_code/` for the topography pipeline behind Table S12)
and derive the number from it. When a number in `ms.qmd` disagrees with that
repo, the repo wins and the manuscript is corrected.

## Supplementary tables vs. the code repo

`sup_tables_manifest.tsv` maps each file in `Sup Tables/` to its source in
`~/github/Liu2026_code_and_data/` (blank if none). `Rscript check_sup_tables.R`
compares them by cell contents, and `ms.qmd` runs the same check at render
time. When a table is added, renamed, or renumbered, update the manifest.

## Supplementary figures vs. the code repo

The supplementary figures are in `sup_figs/`. The folder was renamed from
`Sup Figures/` on 2026-10-07, so older notes, commits, and scripts may still
use the old name.

`sup_figures_manifest.tsv` maps each file in `sup_figs/` to the plots in
`~/github/Liu2026_code_and_data/` it is built from, one row per source plot.
The paper's figures are hand-assembled or edited from those plots, so
`Rscript check_sup_figures.R` compares change times rather than contents. It
flags a figure when a source plot's contents changed (git commit time,
ignoring pure renames) after the paper's copy did. After checking a flagged
figure and finding it needs no update, set `reviewed_through` to that date in
the manifest. `ms.qmd` runs the same check at render time.

## Main figures vs. the code repo

`main_figures_manifest.tsv` maps each file at the top level of
`main_figures/` to the plots in `~/github/Liu2026_code_and_data/` it is built
from, in the same format as `sup_figures_manifest.tsv`.
`Rscript check_main_figures.R` runs the same change-time comparison as
`check_sup_figures.R`, and `ms.qmd` runs it at render time. When a main figure
is added or renamed, update the manifest and `crop_figs.sh`.

## The older/ folder

Ignore everything under `older/` (and `main_figures/older/`). These are
superseded files kept only for history. Do not use them as data sources or
as references when checking figures, tables, or numbers.
