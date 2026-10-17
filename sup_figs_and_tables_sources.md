# Sources of the supplementary tables and figures

Survey of 2026-10-15 across this repo, `~/github/Liu2026_code_and_data/`
(the code repo), and `~/github/unified_indels_data/` (the data repo).

Every supplementary table and figure in the paper exists, and
`check_sup_tables.R` and `check_sup_figures.R` pass. Several items are built
from intermediate files checked into the code repo that nothing in the repo
regenerates.

## Built only from checked-in intermediates (upstream data missing)

1. **Table S14** (89-type replication timing).
   - Mini sent the r.04 matrix
     `04_createRTgroupInf/rtGroup/tissue.rtGroup.Mutation_50cutoff.xlsx` on
     2026-10-16. It is now in `ID89_topography_analysis_code/`, and
     `r.09_barplot.R` there reads it (changed from the `.Rds` copy that r.04
     also writes).
   - Re-run on 2026-10-16: the `trend.xlsx` written by `r.09_barplot.R`
     reproduces the paper's Table S14 exactly (945 rows, same slopes, P
     values, and trends). The script then stops in the per-signature plots,
     which the paper does not use.
   - `build_sup_table_s14/build_sup_table_s14.R` reformats a copy of that
     file, `trend_89.xlsx`, into the table. The paper's copy was replaced by
     the built table, which only corrects the title and names the first
     column "Signature", as in Table S13.
   - The data for the steps before r.04 are still missing. The main input is
     `all.indel.ID89.partial.credit.txt` (per-indel 89-type signature
     probabilities, made by Mo Liu), which
     `ID89_topography_analysis_code/README.md` says is too large to include.

2. **Table S13, Figure S4, and Figure 5D.**
   - `ID83_topography_analysis_code/r.09_barplot.R` (from Mini Huang) makes
     `plot_data.xlsx` and `all_lm_trend.xlsx`. It reads
     `04_create4InfMatrix/rtGroup/tissue.rtGroup.{Mutation,Simulated}_50cutoff.xlsx`.
     Mini sent these two files, and they were added to the code repo on
     2026-10-16.
   - Re-run on 2026-10-16: `plot_data.xlsx` and `all_lm_trend.xlsx` are
     identical to the copies saved as `09_barplot_data.xlsx` and
     `09_all_lm_trend.xlsx` in `build_fig5/` and `build_fig_s4/`. Rebuilt from
     them, `FigureS4.png`, `Figure5D.png`, `Figure5ABC.png`, and `Figure5.png`
     are pixel-identical to the checked-in copies. Table S13 rebuilt from them
     matches the paper's copy (labels identical, numbers agree to a relative
     difference of at most 1e-12).
   - So steps r.04 onward for replication timing now run from the repo. The
     data for the steps before r.04 are missing, as for item 3.
3. **Table S12, Figure S3, and Figure 5 panels A to C.**
   - The 50%-cutoff matrices from step r.04 are checked in at
     `ID83_topography_analysis_code/04_create4InfMatrix/` for `DNA.region`,
     `repli.strand`, and `trans.strand`. So steps r.05 (Fisher test) and r.06
     (odds ratio) can be rerun from the repo.
   - The data for the steps before r.04 are missing. See
     "ID83 topography steps before r.04" below.
   - Re-run of r.05 and r.06 on 2026-10-15, from the checked-in r.04
     matrices:
     - r.06: all 9 output files are identical to the copies in
       `ID83_topography_analysis_code/06_computeOddsRatio/` and
       `build_fig_s3/06_computeOddsRatio/`.
     - r.05: all 3 outputs are identical to the copies in
       `ID83_topography_analysis_code/05_FisherTest/`. Against
       `build_fig_s3/05_FisherTest/` (the input to Table S12), `repli.strand`
       and `trans.strand` are identical. `DNA.region` differs only in the
       `odds_ratio` column, with a mean relative difference of 1e-7.
     - Table S12 rebuilt from the re-run matches the paper's copy. The
       cancer-type and signature columns are identical. P values agree to a
       relative difference of at most 1e-13. Odds ratios agree to at most
       1e-12, except Colon, ID_J, intergenic/genic, which is 2.893254 in the
       re-run against 2.893205 in the paper (relative difference 1.7e-5). The
       cause of that one larger gap was not checked. It does not matter at
       any precision the paper reports.
4. **Tables S15 and S16, Figure S6.** Mo Liu added
   `build_fig_s6/build_fig_s6_stats.R` on 2026-10-16. It writes
   `gender.bias.txt` and `correlation.between.age.and.signature.activity.txt`
   from the data repo's `sample_info.tsv` and assignment tables. Re-run on
   2026-10-16:
   - The age table (Table S16) reproduces exactly. Only the row order of the
     "Overall" rows differs.
   - The sex-bias table (Table S15) does not reproduce. The script tests 186
     cancer-type and signature groups, but the checked-in file and Table S15
     have 160. The 26 extra rows are all ID_N and InsDel_N_beta. The 160
     shared rows have identical counts, P values, and odds ratios. The BH
     FDRs differ (by up to 16%), because more tests are corrected. The same 10
     groups pass FDR < 0.2, so Figure S6A shows the same signatures, with
     slightly different color intensities.
   - The `bonev` column needs the Bioconductor-dependent package BonEV,
     which is not installed here. The script then copies the old values.
5. **Figure S1.** Mo Liu added `build_fig_s1/build_all_meta.R` and
   `build_fig_s1/build_assignment_high_similarity_samples.R` on 2026-10-16.
   Re-run on 2026-10-16, both reproduce the checked-in `all.meta.txt` and
   `assignment.high.similarity.samples.txt`. `build_all_meta.R` had a
   hard-coded `D:/` path and wrote CRLF line endings. After fixing both, its
   output is byte-identical to `all.meta.txt`.
   Duplicated data sources:
   - `all.meta.txt` repeated the data repo's `sample_info.tsv` (age, sex,
     cancer type, and cohort agree for all 7,013 of its tumors). On
     2026-10-16, `Signature Assignment.qmd` was changed to read the cancer
     types from `sample_info.tsv`, which gives identical plot data, and
     `all.meta.txt` was removed. `build_all_meta.R` was moved to the data
     repo as a record of how `sample_info.tsv` was assembled, and its inputs
     were deleted from `build_fig_s1/`.
   - `build_fig_s1/liu_et_al_{83,89}_assignment.tsv` were identical copies
     of the data repo files. On 2026-10-16, `Signature Assignment.qmd` was
     changed to read them from the data repo (identical plot data), and the
     copies were deleted.
   - `build_fig_s1/HMF_metadata.tsv` was the Hartwig clinical metadata
     (patient IDs, biopsy dates and sites), probably covered by the Hartwig
     data access agreement. It was removed from the code repo and its git
     history on 2026-10-16.

## No build script, but the data probably exist

6. **Table S17.** All 6975 patients are in
   `unified_indels_data/sample_info.tsv`, which has 7016 rows. The
   `mutation_burden` column and the ratio of 1-bp T insertions to deletions at
   poly-T tracts of 5 or more should be computable from the spectra. This has
   not been checked.
7. **Supplementary Figure Legend.docx.** This is hand-written, so no source is
   expected.

## Present locally but not publishable

8. **Table S1.** `build_sup_table_s1/build_sup_table_s1_step1.R` reads about
   16.7k annotated VCFs from `~/MEGA/important_mut_sig_data/`, which is
   outside the data repo.
9. **Table S9.** This table reads controlled-access variant tables in
   `build_sup_table_s9/DELETE_BEFORE_PUBLISHING/`.

## Complete from the data repo

Tables S2 to S8, S10, and S11, and Figures S2 and S5.

## ID83 topography steps before r.04

There are six steps before r.04. Their scripts are all in
`ID83_topography_analysis_code/` in the code repo. What is missing is the data
they read, plus the code that made two of those inputs. All the input files
sit under `/public/data/Topography_analysis*/` on the Huang-lab server.

1. **`normalize.sh`** rescales each indel's signature probabilities so they
   sum to 1. It reads `all.indel.ID83.partial.credit.txt`, the per-indel ID83
   signature probabilities, and writes `all.indel.ID83.normalized.txt`. No
   code in either repo makes the input file. The 89-type counterpart came from
   Mo Liu.
2. **`r.01_createBed.R`** builds the BED files for the real indels. It reads:
   - the normalized table from step 1
   - `per_base_territories_20kb.mat` (replication-timing territories)
   - `single.indel.aggregate.strand.added.matrix.txt` (strand reference)
   - `all.annotate.ID.vcfs.RData`
3. **`generate_partialcredit_match_file.R`** makes a lookup table of signature
   probabilities for each sample and indel class, built from step 2's output.
4. **`add_partial_credit.R`** gives each simulated indel the probabilities of
   the matching real sample and indel class. It reads
   `all_simulated_with_annotation.txt`. No code in either repo makes the
   simulated indels.
5. **`r.02.intersect_no_loop.sh`** runs bedtools to overlap the real and
   simulated indels with the replication-timing and replication-strand BEDs.
   It also needs `ENCODE/repli_strand_info_new.bed`.
6. **`r.03_mergeInf.R`** (real) and **`r.03_mergeInf_simulation.R`**
   (simulated) merge the overlaps with tissue and strand information. They read
   `all.meta.0528.rds`, which is checked in under
   `ID89_topography_analysis_code/` but not the ID83 folder.

## What is needed from collaborators

- From Mini Huang (likely): the raw topography inputs (items 2 and 3). The
  `rtGroup` matrices arrived 2026-10-16.
- The code that made `all.indel.ID83.partial.credit.txt` and the code that
  simulated the background indels (see "ID83 topography steps before r.04").
- From Mo Liu: `all.indel.ID89.partial.credit.txt` (item 1). This is no
  longer needed to build Table S14, only to rerun the steps before r.04.
- From Mo Liu: why the checked-in sex-bias table lacks ID_N and
  InsDel_N_beta (item 4), and where `HMF_metadata.tsv` came from (item 5).
