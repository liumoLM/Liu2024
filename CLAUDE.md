# Local notes for the 0Liu2024 repo

## Key data sources for sample-level analyses

- **89-type / 83-type indel signature exposures per sample**:
  `sup_tables/sup_table_s11_83_89_signature_assignment.tsv`
  - Rows are signatures named `C_ID*/InsDel*` (83-type/89-type pair).
  - Columns are 6975 samples (Patient IDs).
  - Cell values are exposures (indel counts). Column sums are NOT total
    indels per sample. Since the assignments of 2026-10-19 were adopted on
    2026-10-22, each indel is counted twice (83-type and 89-type), so every
    column sum is about twice the number of indels. For indel counts, use
    the column sums of `unified_indels_data/liu_et_al_83_spectra.tsv`.

- **Sample metadata (cancer type, MSI status, etc.)**:
  `sup_tables/sup_table_s16_sample_metadata.xlsx`
  - Key columns: `Patient`, `Major Cancer Type`, `MSIseq` (values `MSI-H` /
    `Non-MSI-H`), `mutation_burden`, `cohort`, `Cancer Type`. `MSI_status`
    (the cohort calls) was dropped on 2026-10-30.
  - All 6975 Patient IDs match the column names in Table S11.
  - Built by `build_sup_table_s16/build_sup_table_s16.R` in the code repo
    from `sample_info.tsv` and `liu_et_al_83_spectra.tsv`.

## MSI-H and hypermutation definitions

Decided on 2026-10-16 and 2026-10-17 (with Mo), and used everywhere in the
paper and in `explore-attributions/`:

- **MSI-H** is the MSIseq call: `MSIseq` in Table S16 (`MSI-H` /
  `Non-MSI-H`), or `MSIseq_MSI.H` in `unified_indels_data/sample_info.tsv`
  (204 MSI-H tumors since 2026-10-30, from `msiseq/call_msi_h.R` in the
  code repo; 202 before). Never use `MSI_status`, which holds the cohort calls
  from PCAWG and Hartwig clinical data (91 MSI). Those are not reliable.
- **Hypermutated (MSS-hyper)**: not MSI-H by MSIseq and >= 4,000 indels
  (51 tumors; 48 before 2026-10-30).
- **Indel counts** always come from the spectra in `unified_indels_data`,
  which are clipped to repeat length R <= 9 (column sums of
  `liu_et_al_83_spectra.tsv`). Never use Table S11 column sums (which count
  most indels twice) as indel counts. Since 2026-10-17, Table S16
  `mutation_burden` and the poly-T ratio are computed from these spectra by
  `build_sup_table_s16/build_sup_table_s16.R` in the code repo. Earlier
  versions of S16 held larger counts from an older version of the data.

Figure 4A is built with these rules by `build_fig4/all_cancer_types_enrichment.R`
in the code repo. Tumors were refit with mSigAct if they had >= 4,000 indels,
or a poly-T ratio > 0.5 and a SigProfilerAssignment reconstruction cosine
similarity < 0.9 (731 tumors, listed in `mSigAct.reanalyzed.samples.txt`).

The MSS-hyper enrichment analysis is now Figure 4A. The name "Figure 5
panel C" in `explore-attributions/` is an older figure number. The current
Figure 5C is genic vs. intergenic enrichment. In the exploratory version,
samples are assigned to mutually exclusive groups for the Fisher
enrichment analysis in this priority order:

1. **MSI-H**: MSI-H by MSIseq.
2. **MSS-hyper**: not MSI-H and >= 4,000 indels.
3. **Major cancer type** (`Bladder`, `Colon`, `Esophagus`, `Kidney`,
   `Liver`, `Lung`, `Ovary`, `Prostate`, `Skin`): not MSI-H, < 4,000
   indels, and `Major Cancer Type` matches.

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

## File names for figures and tables

Since 2026-10-17, each figure and table in this repo has the same file name
as its source in the code repo: number first, then a short description, all
lower case with underscores. Examples: `main_figures/fig5.pdf` from
`build_fig5/fig5.pdf`, `sup_figs/sup_fig_s4_replication_timing.pdf` from
`build_fig_s4/`, and `sup_tables/sup_table_s16_sample_metadata.xlsx` from
`build_sup_table_s16/`. A figure edited by hand comes from a source with a
`_draft` suffix. Keep this rule for new or renumbered figures and tables, and
rename in both repos. `sup_tables/` was `Sup Tables/` before that date.

## Supplementary tables vs. the code repo

`sup_tables_manifest.tsv` maps each file in `sup_tables/` to its source in
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
