# Local notes for the 0Liu2024 repo

## Key data sources for sample-level analyses

- **89-type / 83-type indel signature exposures per sample**:
  `Sup Tables/Table S11 83-type and 89-type signature assignment.tsv`
  - Rows are signatures named `C_ID*/InsDel*` (83-type/89-type pair).
  - Columns are 6975 samples (Patient IDs).
  - Cell values are exposures (indel counts). Column sums = total indels per sample.

- **Sample metadata (cancer type, MSI status, etc.)**:
  `Sup Tables/Table S18 metadata of 6975 samples.xlsx`
  - Key columns: `Patient`, `Major Cancer Type`, `MSI_status` (values `MSI` / `MSS`),
    `mutation_burden`, `cohort`, `Cancer Type`.
  - All 6975 Patient IDs match the column names in Table S11.

## Sample grouping convention for Figure 5 panel C (with MSS-hyper)

When stratifying samples into mutually exclusive groups for the panel C
Fisher enrichment analysis, assign in this priority order so each sample
lands in exactly one column:

1. **MSI-H**: `MSI_status == "MSI"`.
2. **MSS-hyper**: `MSS_status == "MSS"` AND total indels (column sum in
   Table S11) > 5000.
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
`ID83_topography_analysis_code/` for the topography pipeline behind Table S13)
and derive the number from it. When a number in `ms.qmd` disagrees with that
repo, the repo wins and the manuscript is corrected.
