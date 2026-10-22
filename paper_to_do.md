# Journal requirements (cell resources)

- 7 display items allowed (We will have 1 table and 6 main figures)
- Word Count: The main text should be under 7,000 words. This count includes the figure legends but excludes the STAR Methods, supplemental legends, and references.
- Supplemental Items: You may include additional display items as supplemental information at the editor's discretion. For Cell Research Articles, this is often up to 7 supplemental figures.
- Mandatory Items: You must also include a Graphical Abstract (exactly 1,200 x 1,200 pixels) and 3–4 Highlights (each 85 characters or fewer). [1, 2, 3, 4]
http://8.138.15.71/


# To do

## Current (updated 2026-10-22)

- Done (2026-10-17): Item 6, `build_fig4/all_cancer_types_enrichment.R` now
  calls a tumor hypermutated if it has > 4,000 indels in the R <= 9 clipped
  spectra (column sums of `liu_et_al_83_spectra.tsv`). Figure 4A was
  recomputed and copied to `main_figures/`. MSS-hyper now has 48 tumors
  instead of 723, and its only dot is InsDel3b instead of InsDel14. The
  legend already says "> 4,000 indels", and the results text does not
  mention the MSS-hyper column.
- Done (2026-10-18): Item 6, the same definitions are used everywhere
  else. The scripts in `explore-attributions/` were switched on 2026-10-17
  (ce9a4f3). On 2026-10-18 the numbers in `id1_id2_idN_findings.md` and
  `more-on-attributions.md` were updated from the regenerated outputs
  (several findings about MSS hypermutated tumors no longer hold), and the
  two scripts that reproduce Mo's original figure with `MSI_status` say so
  in a comment. `ms.qmd` and the code repo already used them.
- Done (2026-10-17): Item 5, Table S17 is now built by
  `build_sup_table_s17/build_sup_table_s17.R` in the code repo.
  `mutation_burden` and the poly-T ratio are computed from the clipped
  spectra, and the other columns reproduce exactly. The legend now gives the
  channels in the poly-T ratio (INS:T:1:5+ and DEL:T:1:5+, so deletions
  are from runs of 6 or more Ts) and no longer says that tumors with a ratio
  of 0.5 or more went to mSigAct. The methods wording about poly-T tracts
  was fixed on 2026-10-18.
- New signature assignments (Mo, email "Question about S11",
  2026-10-19). Mo could not reproduce some of the current assignments
  (probably because no random seed was set), so the assignments must be
  redone, some figures regenerated, and the topography analysis redone.
  Mo reconstructed the original workflow: SigProfilerAssignment on the
  172-type catalogs, then mSigAct reruns for (1) tumors with > 4,000
  indels and (2) tumors with poly-T ratio > 0.5 and a SigProfilerAssignment
  reconstruction cosine similarity < 0.9 in either the 83-type or the
  89-type fit. That is 731 tumors, and 784 tumors have different
  attributions than now. Mo thinks the overall effect is limited.
  - Done (2026-10-22): the four attachments of Mo's 2026-10-19 email are
    in `unified_indels_data/assignments_2026-10-19/` (not yet committed):
    `ID83.SPA.assignment.1020.txt`, `ID89.SPA.assignment.1020.txt`,
    `assignment.high.similarity.samples.for.figureS1.txt` (6,516 tumors),
    and the new `mSigAct.reanalyzed.samples.txt` (731 tumors). The new
    assignment tables have the same 6,975 tumors and signatures as
    `liu_et_al_{83,89}_assignment.tsv`, and 819 tumors have exposures that
    differ by more than 0.5 (Mo said 784). The list has 248 of the 249
    tumors with > 4,000 clipped indels. The missing one is the exome tumor
    `WES::TCGA-EW-A2FV-01A-11D-A17D-09` (4,198 indels). Ask Mo whether
    leaving it out was deliberate. The copies of
    `mSigAct.reanalyzed.samples.txt` in this repo and in the code repo
    (untracked) are the OLD 2026-10-18 list (1,535 tumors, S17
    `mutation_burden` >= 4,000 or poly-T ratio > 0.5) and should be
    replaced.
  - Steve agreed (2026-10-19) to rerun the figure and table pipeline with
    the new assignments, as an independent check of how much the results
    change. Mo has not yet pushed them to the code repo (checked
    2026-10-22).
  - Then update the methods to the new routing rule (they now give the old
    1,535-tumor rule).
  - Mo will redo the topography partial credit (needed for Table S14,
    `all.indel.ID89.partial.credit.txt`) once Steve says the new
    attributions are OK.
  - The Table S11 question (503 tumors with single-counted sums,
    `s11_single_counted_503.tsv`) should be rechecked against the new
    assignments.
- Meeting with Mo: Thursday 2026-10-22, 8 PM US time
  (https://meet.google.com/gnh-mwye-jgd).
- MSIseq. Mo said (2026-10-18) that MSIseq was run on all tumors, but Mini
  wrote (2026-10-19) that it was run on only 220 tumors (probably those
  above some mutation-load cutoff), with `uniform.seq.len = 3200` and the
  repeats annotation `Hg19repeats.rda` from the MSIseq package. Mini's code
  is in the code repo at `msiseq/run-MSIseq.R` (not yet committed). Mo is
  preparing SBS and indel VCFs so Mini can run MSIseq on all tumors. No
  news since 2026-10-19. Ask Mo whether the rerun is needed. If it is, the
  methods must say what MSIseq was run on, and the MSI-H calls
  (`MSIseq_MSI.H` in `sample_info.tsv`, Table S17, Figure 4A, Table S9)
  may change.
- Answered (Mo, 2026-10-18): the high-TMB set in the de novo extraction
  used > 4,000 indels after clipping to R <= 9, and the extraction was
  redone after clipping. Still to do: give this cutoff in the methods
  (`ms.qmd:969`, `:975`, `:981`, `:991`).
- Earlier items from Mo's replies:
  - Answered and done (2026-10-17): ID_N and InsDel_N_beta in Table S15
    (Table S15 and Figure S6 regenerated), and `donor.all_projects.tsv.gz`
    (the ICGC portal download of 9 May 2024 named in the data repo README).
- Steve: find the Hartwig Medical Foundation Data Transfer Agreement (or
  similar), probably the older one signed in Singapore, and get the
  release version or data request ID for the Hartwig metadata and VCFs.
  Mo (2026-10-17): both `HMF_metadata.tsv` and the VCFs came from the
  Duke-NUS server, where Arnoud Boot had downloaded them, and Mo does not
  know the release or request ID. Then add it to the header of
  `build_all_meta.R` and to the `sample_info.tsv` entry of the data repo
  README, which has a "TO DO" in capitals for it.

Steve review Mini's methods -

Steve Review subsection "Demographic Associations with Gender and Age"


check Shiny app persistent doi: 10.5281/zenodo.20092288
XmWU123 Github name


Steve: check "self" referencing from web page to shiny app doi

steve check: Sup Fig 1 legend

steve: review all sup fig legends

- Done (2026-10-22): Table S9 no longer has `MSI_status` ("MSI status" in
  the main sheet, "MSI status = MSI" in the summary sheet), following the
  decision to use only MSIseq. `build_sup_table_s9/build_sup_table_s9_pol_proofreading.R`
  in the code repo no longer writes them, and its comments and README say
  not to add them back. It also now spells "tumor". The table was rebuilt
  and copied to `sup_tables/`. Apart from the dropped columns and the
  spelling, the contents are unchanged (248 tumors, 200 MSI-H), and
  `check_sup_tables.R` passes. "Total somatic indels" is the column sum of
  the clipped 89-type spectra.

# To do: figures, sup figures, and sup tables vs. the code repo (2026-10-07)

Run `Rscript check_main_figures.R`, `Rscript check_sup_figures.R`, and
`Rscript check_sup_tables.R` to see the current state.

1. Closed (2026-10-17): read Mo's email replies. Follow-ups are in the
   "Tomorrow" section above.
2. Topography annotations and simulations. Answered (Mo, 2026-10-07 and
   2026-10-08, and Steve, 2026-10-17): `dna.region` and `trans.strand` came
   from ICAMS with `trans.ranges.GRCh37` (GENCODE 30 CCDS genes), genome
   build hg19 (GRCh37). SigProfilerSimulator 1.2.2 was run as
   `sigSim.SigProfilerSimulator("project", "/path/to/project/", "GRCh37",
   contexts=["ID"], simulations=10)`. The strand column (+, -, Q) is
   computed by SigProfilerSimulator (`mutational_simulator.py`, lines
   1134-1151). Decision: 10 simulations per tumor is enough. Mo asked that
   we cite Otlu et al. 2023 and Bergstrom et al. 2020.
   Still to do:
   - Waiting on Mo: the topography code keeps indels in genes on both
     strands (114,297 indels, 0.6%), although the methods say they were
     ignored. Mo is probably redoing the simulations and analysis with code
     that drops them, with Mini. Rerun the downstream topography tables and
     figures when that arrives.
   - Done (2026-10-17): filled the placeholders in the green "Simulating
     synthetic cancer datasets" and "Annotating somatic indels as genic or
     intergenic" paragraphs in `ms.qmd`, and updated the blue strand note
     with Mo's answer.
3. Mini replied 2026-10-08. `ID83_topography_analysis_code/r.09_barplot.R`
   (now in the code repo) creates `09_barplot_data.xlsx` and
   `09_all_lm_trend.xlsx`, and reproduces all 714 trends in Table S13.
   Table S13 is confirmed to come from `09_all_lm_trend.xlsx`. The 89-type
   scripts are in `ID89_topography_analysis_code/`. Mini later sent the r.04
   inputs `tissue.rtGroup.{Mutation,Simulated}_50cutoff.xlsx`, now in the
   code repo. With them, r.09 reproduces Table S13, Figure S4, and Figure 5D
   exactly (2026-10-16). Mini also sent the 89-type r.04 matrix, and the
   ID89 r.09 reproduces Table S14 exactly. `build_sup_table_s14.R` now
   builds the table (2026-10-16). Still missing, but only needed to rerun
   the steps before r.04: the 89-type per-indel input
   `all.indel.ID89.partial.credit.txt` (too big to email). In `ms.qmd`, the
   89-type replication-timing result now cites only Table S14, since Figure
   S4 shows only 83-type signatures (2026-10-16). The
   replication-timing methods and the Fig 5D
   legend now follow r.09_barplot.R (500-indel minimum, trend rule), in
   green for review (2026-10-08).
4. Closed (2026-10-07): Fig S2. Mo has no assembly code, so we wrote
   `build_fig_s2/make_figure_s2_layout.R` in the code repo. It builds
   `sup_fig_s2_msighdp_signatures_draft.pdf` in the paper's layout from the Table S2, S4, and
   S6 signatures. The paper's `sup_figs/sup_fig_s2_msighdp_signatures.pdf` is edited by hand
   from the draft, as recorded in `sup_figures_manifest.tsv`.
5. Table S17 (sample metadata) is only partly reproducible. Patient, age,
   sex, MSI status, major cancer type, and cohort come from
   `unified_indels_data/sample_info.tsv`, which also has 41 samples not in
   S17. `mutation_burden` and the polyT ratio do not match the current
   spectra (mutation_burden equals the spectrum total for only 486 of 6,975
   tumors and is otherwise larger). They probably come from an earlier,
   pre-cap-9 version of the data. `MSIseq` matches `MSIseq_MSI.H` in
   `sample_info.tsv` for all 6,975 tumors. `Cancer Type` matches
   `Original_Cancer_Type` there, except that Colon, Head, and
   Bone-SoftTissue are renamed Colon/Rectum, Head and neck, and Bone/Soft
   tissue (2026-10-07). S17 has a "legend" sheet.
   Decision (Mo, 2026-10-16): mutation counts come from the spectra in
   `unified_indels_data`, which are clipped to repeat length R <= 9 (column
   sums of `liu_et_al_83_spectra.tsv`). Done (2026-10-17):
   `build_sup_table_s17/build_sup_table_s17.R` in the code repo builds S17
   with `mutation_burden` and the poly-T ratio from those spectra.
   Routing: Mo's list of mSigAct tumors (2026-10-18) follows the 4,000
   threshold on the clipped counts, and the methods now say so (see the
   routing entry above).
6. Figure 4A MSS-hyper threshold. Decision (Mo, 2026-10-16): a tumor is
   hypermutated if it has > 4,000 indels, counted from the R <= 9 clipped
   spectra in `unified_indels_data`. Done (2026-10-17): Figure 4A was
   recomputed with this rule (48 MSS-hyper tumors instead of 723), and the
   project CLAUDE.md records it. Done (2026-10-18): the same definition
   is used everywhere else (see the "Tomorrow" section).
7. Closed (2026-10-17): Figures S1 and S6, and Tables S15 and S16.
   - Figure S1 reproduces. It reads `sample_info.tsv` and the assignment
     tables from the data repo, and `build_all_meta.R` moved to the data
     repo as a record of how `sample_info.tsv` was assembled.
   - `HMF_metadata.tsv` (controlled-access Hartwig clinical data) was
     removed from the code repo and its git history, and the repo was
     force-pushed. Mo, Mini, and Xueming were emailed to reset their clones.
   - Table S16 reproduces. Table S15 and Figure S6A were regenerated with
     ID_N and InsDel_N_beta, which Mo had left out by mistake (186 groups,
     same 10 at FDR < 0.2).
   - Figure S6 is now assembled in R by
     `build_fig_s6/make_figure_s6_layout.R`, replacing the Illustrator
     version. The sex-bias text in `ms.qmd` was updated to match.
   Details in `sup_figs_and_tables_sources.md`, items 4 and 5.
8. Review the STAR Methods for the topography analysis.
9. Closed (2026-10-17): SigProfilerTopography
   (https://github.com/AlexandrovLab/SigProfilerTopography) uses the same
   indel strand rule as we do. It reads the strand from the `seqInfo`
   output of SigProfilerMatrixGenerator (`MutationMatrixGenerator.py`,
   lines 1427-1443 and 1789-1799 at commit e67c5c4). There, an indel whose
   bases are only C and T is on the + strand, only G and A on the - strand,
   and mixed (for example GAGAGA->GAGA) gets strand 0 (Q), which
   SigProfilerTopography excludes (`ReplicationStrandBiasAnalysis.py`,
   lines 389-396 at commit bff7a86). Recorded in the blue strand note in
   `ms.qmd`.
10. Closed (2026-10-17): file names aligned across the two repos. Each
   paper figure or table has the same name as its source in the code repo,
   number first, then a short description: `main_figures/fig5.pdf` from
   `build_fig5/fig5.pdf`, `sup_figs/sup_fig_s3a_transcription_strand.pdf`
   from `build_fig_s3/`, and `sup_tables/sup_table_s2_msighdp_83_type_signatures.xlsx`
   from `build_sup_table_s2/`. Hand-edited figures come from a `_draft`
   source (for example `sup_fig_s5_rnase_h2_spectra_draft.pdf`).
   `Sup Tables/` was renamed `sup_tables/`. Renumbering a figure or table
   now means renaming it in both repos and in the manifests.

# Data sources for Degaspari SBS data; make sure we cite this
  https://zenodo.org/records/5571551/files/SBS_v2.03.zip?download=1                                                                                                                        
