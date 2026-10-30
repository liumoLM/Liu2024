# Journal requirements (cell resources)

- 7 display items allowed (We will have 1 table and 6 main figures)
- Word Count: The main text should be under 7,000 words. This count includes the figure legends but excludes the STAR Methods, supplemental legends, and references.
- Supplemental Items: You may include additional display items as supplemental information at the editor's discretion. For Cell Research Articles, this is often up to 7 supplemental figures.
- Mandatory Items: You must also include a Graphical Abstract (exactly 1,200 x 1,200 pixels) and 3–4 Highlights (each 85 characters or fewer). [1, 2, 3, 4]
http://8.138.15.71/


# To do (updated 2026-10-23)

## Before submission and before posting on bioRxiv

- Release the code as a new GitHub repo, not `Liu2026_code_and_data` itself
  (decided 2026-10-24). The working repo's git history still holds
  controlled-access per-variant data
  (`build_fig6/ID4_IDF_linking_tumor_indel_annotated_vcf.txt`, untracked
  since 2026-10-24), so the working repo stays private and keeps its history
  for internal use.
  - Make the new repo from the tracked files only, with no history:
    ```
    git -C ~/github/Liu2026_code_and_data archive --prefix=Liu2026_release/ HEAD | tar -x -C ~/github
    git -C ~/github/Liu2026_release init && git -C ~/github/Liu2026_release add -A && git -C ~/github/Liu2026_release commit -m "Initial release"
    ```
  - Use `git archive`, not a copy of the folder without `.git`: it copies only
    tracked files (278 MB of the folder's 5.5 GB on 2026-10-24), and leaves
    out untracked and ignored files, such as the old Figure 6 file, scratch
    output, and renders. The repo has no Git LFS files, so no LFS pointer
    files are copied.
  - Not a shallow clone: it keeps the old commit IDs, can be deepened with
    `git fetch --unshallow`, and GitHub often refuses a push from a shallow
    clone to a new remote ("shallow update not allowed"). It gives one
    commit anyway.
  - Before making it, check `git ls-files` in the working repo for anything
    else controlled access, because every tracked file is copied.
  - Before running `git archive`, replace Steve's local paths
    (`~/MEGA/important_mut_sig_data/...`) in `local_dir` in
    `find_annotated_vcfs()` and `find_partial_credit()` in `find_data.R`
    with placeholders, such as
    `"path/to/hartwig/annotated/vcfs"`. Also check
    the other scripts that hard-code `~/MEGA` paths (`rg -l MEGA` in the
    code repo, for example in `build_sup_table_s9/` and
    `build_sup_table_s1/`).
  - Leave out `ID89_topography_analysis_code/likely_dead_ID89_topography_code/`
    (plots from the ID89 r.09 that the paper does not use, moved there on
    2026-10-24). The code repo's `.gitattributes` marks it `export-ignore`,
    so `git archive` already skips it. The same holds for
    `msiseq/likely_dead_msiseq_code/` (the old MSIseq script, moved there on
    2026-10-24) and `build_sup_table_s1/likely_dead_sup_table_s1_files/` (two
    unread Del3:U1:R1 example files, moved there on 2026-10-25). Check that they are absent from the release repo, and treat
    any other `likely_dead_*` folder the same way.
  - Fixes made later in the working repo must be copied to the release repo
    again (rerun `git archive` and commit the differences).
  - Decide whether the data repo `unified_indels_data` needs the same
    treatment.

## Open: Steve

- `sample_info.tsv` in the data repo has CPCT02220078T (Hartwig, Bladder),
  which has no spectrum in `liu_et_al_83_spectra.tsv` or
  `liu_et_al_89_spectra.tsv` (found 2026-10-26). Decide whether to drop it
  from `sample_info.tsv` or find out why it has no spectrum.
- Dig into the age and sex models in detail (the doc "Age and sex models
  for Tables S15 and S16" (old numbers), https://claude.ai/artifact/X8iSQz1TiCAw89HgKkwHdM,
  exported as `age-and-sex-model.md`), and then review and rewrite the
  subsection "Demographic Associations with Sex and Age" in `ms.qmd` to
  match.
- Presence thresholds differ: Figures 4A and S1 count a signature as
  present above 5% of the tumor's indels, but Tables S14 and S15
  (`presence_threshold <- 0.01` in `build_fig_s6/build_fig_s6_stats.R`) use
  1%, and the S15 Spearman correlations use exposure > 0. Checked
  2026-10-29 with exposures under 5% set to 0: all 10 S14 hits at q < 0.2
  remain (plus Other / C_ID19), 39 of the 41 S15 by-cancer-type hits remain
  (Colon C_ID4 / InsDel4 lost), and the age effects in the text are
  unchanged. In the S15 cancer-type-adjusted models some FDR < 0.05 calls
  change (exposure 31 to 26, presence 17 to 15, burden-adjusted presence 16
  to 22), none cited in `ms.qmd`. Decide with the age and sex items whether
  to use 5% throughout.
- Review Mini's methods.
- Get the Hartwig data download date or version (release) number, for the
  data availability statement and for the header of `build_all_meta.R` and
  the `sample_info.tsv` entry of the data repo README.
- MSI status stated consistently (from `QUESTIONS_FOR_MO.md`, 2026-10-15).
  Following the decision to use only MSIseq: drop `MSI_status` from Table
  S16 (as was done for Table S9; `build_sup_table_s16/build_sup_table_s16.R`
  in the code repo), and cut the methods section "MSI/MSS status and
  high/low TMB status" in `ms.qmd`, which still describes the PCAWG and
  Hartwig cohort MSI calls (it is already marked "NEEDS DRASTIC CUTTING").
- Check the Figure S1 legend, then review all supplementary figure legends.
- Review the STAR Methods for the topography analysis.
- Add the STAR Methods "Resource availability" section (lead contact,
  materials availability, data and code availability), listing the code
  repo, the data repo, and the browser DOI.

- Problems found while building `data_flow/` (2026-10-23):
  - Fixed (2026-10-24): the example SigProfiler scripts in
    `Extraction_and_assignment_code/` could not run. Named
    `SigProfilerExtractor.py` and `SigProfilerAssignment.py`, they imported
    themselves instead of the packages. Renamed `run_SigProfilerExtractor.py`
    and `run_SigProfilerAssignment.py`. In the extractor, only the import was
    inside the cancer-type loop, so only `All` would have run. Fixed and
    checked with stub packages, then with the real packages (2026-10-24).
    On 2026-10-25 Mo confirmed v1.2.1 and 20 to 40 signatures for all
    genomes, all PCAWG, and all Hartwig genomes (2 to 20 per cancer type and
    high-TMB), and the methods were corrected (they said v1.1.24 and 10 to
    40). The extractor example now uses only the 198 bladder tumors (the
    6,975-tumor `All` catalog was removed), has a `--quick` option, and
    `Extraction_and_assignment_code/pixi.toml` pins Mo's versions. Mo's
    numpy 1.26.4 and pandas 2.0.3 cannot be installed with these
    SigProfiler versions, which require numpy 2, so they are not pinned.
    Email to Mo about this sent 2026-10-25. Decided 2026-10-26: numpy 2.2.6
    and pandas 2.3.3 reproduce Mo's SPA exposures exactly, so the methods
    give them for assignment. For extraction they are unknown (numpy 2 or
    later), and the methods give no numpy or pandas versions for it. Still open (blue notes in the
    methods): the extraction script and code repo READMEs define high-TMB
    as at least 5,000 indels, but the methods say more than 4,000. Resolved
    2026-10-26: Mo confirmed 4,000, and the code repo comments were
    corrected. Fixed
    (2026-10-25): the example assignment scripts and READMEs gave the old
    mSigAct criteria (at least 5,000 indels, or poly-T at least 90%, with SPA
    skipping those genomes). They now give the adopted ones (SPA on all
    genomes, then mSigAct for more than 4,000 indels, or poly-T more than 50%
    and cosine less than 0.9).
    Fixed (2026-10-25): `Example_data/assignment.example.data/` held copies
    of the full 172-type catalog (all 6,975 tumors) and signatures. The
    catalog copy had no row names, so SigProfilerAssignment silently dropped
    CPCT02010349T. Removed. The new
    `Extraction_and_assignment_code/make_172_type_inputs.R` builds both from
    the data repo, and both assignment scripts read its output. Rerun:
    6,243 of the 6,244 tumors not refit with mSigAct match Mo's
    assignments (`liu_et_al_89_assignment.tsv`) exactly. CPCT02010349T (144 indels,
    cosine 0.47) is InsDel5a in the rerun and InsDel19a in Mo's assignments.
    On Mo's advice the rerun's values were adopted on 2026-10-26. Open:
    `mSigAct.R` refits all tumors and does not merge results. The code that
    merged the refits is not in the code repo (the selection rule
    reproduces the 731 tumors exactly).
    See the two items of 2026-10-25 under "Open: Steve with Mo".
  - Fixed (2026-10-24): the data repo had a copy of
    `make_connection_table.R` that sourced a `code/` folder that no longer
    exists. Deleted. The data repo README now points to the code repo's
    `vignette/make_connection_table.R`, which reproduces
    `connection_table.tsv` exactly.
  - Fixed (2026-10-24): `msiseq/run-MSIseq.R` was Mini's earlier
    220-tumor script, which cannot run (undefined objects, a `source()` path
    from `/`, inputs not in the repo) and which nothing reads. Moved to
    `msiseq/likely_dead_msiseq_code/prev-run-MSIseq.R` with a header, and
    `msiseq/README.md` added. The script behind the adopted calls goes in
    `msiseq/` after the MSIseq item under "Open: Steve with Mo" is settled.
  - Fixed (2026-10-24): three topography scripts stopped with errors or
    pointed at the wrong data. `r.03_mergeInf_simulation.R` (ID83) ended in a
    shell command, now commented out. The ID89 `r.09_barplot.R` now ends
    after writing `trend.xlsx` and `plot_data.xlsx` (identical contents) and
    exits cleanly. Its unused plots are in
    `ID89_topography_analysis_code/likely_dead_ID89_topography_code/`. The
    ID89 `r.02.intersect_no_loop.sh` was a copy of the ID83 one with ID83
    paths. Fixed (2026-10-25): replaced with the corrected version from
    Mini, which uses the `results0322/ID89/` folders that
    `r.03_mergeInf.R` reads. Only `01.rtGroup.bed` comes from the ID83
    results.
  - Open (2026-10-24), emailed Mini and Mo, waiting on them: the ID83 topography steps
    `normalize.sh` to r.04 ran on the Huang-lab server and disagree on input
    paths, `results/` vs `results0322/`, and columns (`ID.class` vs
    `COSMIC_83`, `ID_O`/`ID_P`). Documented in "Known inconsistencies" in
    `ID83_topography_analysis_code/README.md`, script headers, and
    `data_flow/`. Fix the paths and columns once they confirm which files
    and script versions made the r.04 matrices.
  - Fixed (2026-10-24): Figures 5, S3, and S4 and Tables S12 to S14 read
    hand copies of topography outputs (`build_fig5/09_*.xlsx`,
    `build_fig_s4/09_*.xlsx`, `build_sup_table_s14/trend_89.xlsx`,
    `build_fig_s3/05_FisherTest/`, `build_fig_s3/06_computeOddsRatio/`).
    They now read the pipeline outputs in `ID83_topography_analysis_code/`
    and `ID89_topography_analysis_code/` directly, and the copies are
    deleted. Both `r.09_barplot.R` scripts were rerun in the repo, and their
    `09_barplot/*.xlsx` tables are now tracked. The figures are
    pixel-identical. Tables S12 and S13 differ only in the 13th to 15th
    significant digit (the copies came from an earlier run), and the paper's
    copies were updated. After the topography rerun, rerun r.05, r.06, and
    r.09, then the figure and table scripts.
  - Fixed (2026-10-23): `build_fig6/Plot signatures and motifs.qmd` read
    `liu_et_al_{83,89,476}_signatures.tsv` from
    `external_data/reijns_et_al_2022/` in the data repo, where they are not,
    so it failed. It now reads them from the top level and renders. The
    sequence logos are identical, and the signature plots differ only by
    font, so the committed PDFs were kept.
  - Fixed (2026-10-24): `build_fig6/ID4_IDF_linking_tumor_indel_annotated_vcf.txt`
    held per-variant indels (with coordinates and sequence context) of
    controlled-access tumors, labeled by signature, so it may be
    re-identifiable. Figure 6 now draws its logos from
    `fig6_logo_counts.tsv` (per-position base counts), written from the
    annotated VCFs by `make_fig6_logo_counts.R`. The file is no longer
    tracked, but it is still in the code repo's git history. So the code
    repo itself will not be released. The release will be a new repo made
    from the tracked files with no history (`git archive HEAD`, then
    `git init`). Data repo `linking_tumors/README.md` now says the file is
    not released. Matching its
    rows to the VCFs showed it came from 3 tumors per signature (C_ID4:
    CPCT02100089T, SP48008, SP124364; ID_F: DRUP01030028T, SP103894,
    SP5559), not the two linking tumors, and not the 6 C_ID4 tumors listed
    in the notebook. 6 of its 4,066 logo deletions are not in the current
    VCFs, so 3 of the 4 logos differ slightly from the paper's panel C.

## Open: Steve with Mo

- mSigAct refit (added 2026-10-25, spot checked 2026-10-26). Mo's code
  (`mSigAct_for_rerun.R`, from his email of 2026-10-26) is cleaned up
  as `Extraction_and_assignment_code/run_mSigAct_refit.R` in the code repo,
  replacing `mSigAct.R`. It refits only the 89-type part of the spectrum and
  doubles the exposures. Run on 16 of the 731 refit tumors, it reproduced
  all their exposures in `liu_et_al_8?_assignment.tsv`, except InsDel1b of
  CPCT02340022T (13,340 vs 13,342). Open: confirm with Mo that he used
  mSigAct2 1.0.0 (github.com/Rozen-Lab/mSigAct2, branch v3.0.3-branch,
  commit b3e11fc) and mSigTools 1.0.8, then fix the methods (blue note).
- Age and sex analyses (Tables S14 and S15, Figure S6), discussed in the
  doc (see "Open: Steve"). Decide:
  - whether to use the joint models (age, sex, cohort, MSI-H, cancer type,
    `build_fig_s6/age_sex_glm_prototype.R`) for Table S15, or keep the
    age + cancer type models now in `build_fig_s6_stats.R`;
  - whether to replace the Table S14 Fisher tests with per-type Firth
    logistic models adjusted for age, cohort and MSI-H
    (`build_fig_s6/sex_glm_prototype.R`). The S14 associations at q < 0.2
    all hold after adjustment;
  - what to do with sex and age associations in CNS, lymphoid and
    bone/soft tissue, which mix subtypes that differ in age, sex and
    signature frequency (drop them, adjust for subtype where known, or
    report with a caveat);
  - whether to follow up ID_N and InsDel_N_beta (strong age x cancer type
    interaction, exposure and presence in opposite directions);
  - who writes the methods paragraph and the Table S15 legend.
- MSIseq. Mini reran MSIseq on all 6,975 tumors (Mo's email "The new
  MSIseq result", 2026-10-23). The results, Mini's script
  (`run_MSIseq_overlap_6975.R`), the earlier 220-tumor script
  (`prev-run-MSIseq.R`), and the emails are in
  `unified_indels_data/new-mini-msi-calls/` (committed 2026-10-23 as
  provisional: if Mini and Mo say these are not the right files, move them
  and delete the folder). The new
  calls use a rule set by hand, S.ind > 0.4 and T.ind > 1.2: 204 MSI-H
  (202 before, 195 in both, 7 dropped, 9 added). Checked 2026-10-23: the
  6,975 tumors are exactly those of Tables S11 and S16, and the indel count
  per tumor behind T.ind (T.ind x 3,200) equals the cap-9 spectrum total
  for all of them, so the indel input matches the spectra. Steve emailed
  Mini and Mo (2026-10-23) asking whether `Hg19repeats.rda` is the MSIseq
  package file (answered 2026-10-29: the package help points to
  https://steverozen.net/data/Hg19repeats.rda, and with that file
  `Compute.input.variables()` reproduces Mini's T.ind and S.ind exactly
  for 199 tumors, 104 of them near the cutoffs), whether all Hartwig SBS and indel rows
  are FILTER = PASS (answered by Steve 2026-10-29: yes, only FILTER = PASS
  rows are used for SBS and indels), and (most important) why the hand-set rule was used rather than
  `MSIseq.classify()`, how the two compare, and for
  `MSI_classified_overlap_6975.tsv`. Waiting for the answer. Then: adopt
  the new calls in `sample_info.tsv`, rebuild Figure 4A, Tables S9, S10,
  and S16, update `ms_numbers.R` and `ms.qmd` (202 MSI-H, 200 of 202
  refit), the project CLAUDE.md, and the MSIseq methods (run on all
  tumors, the inputs, the classification rule, and Mo's request to drop
  the sentence on how the 4,000 cutoff was chosen). Also copy the script
  that made the adopted calls into the code repo's `msiseq/`, with its paths
  set through `find_data.R`, and add it to `data_flow/`.
  Found 2026-10-29: no code applies the S.ind > 0.4 and T.ind > 1.2 rule.
  Mini's script calls `MSIseq.classify()`, and the rule column was added by
  hand. `MSIseq.classify()` has no threshold argument. Its built-in
  classifier (MSIseq 1.0.0) is the single split S.ind > 0.394737, which gives
  274 MSI-H on our tumors. The rule's 204 are a subset of these, so the added
  T.ind > 1.2 drops 70 borderline tumors. To ask Mo at the meeting of
  2026-10-29 (`questions_for_mo_2026-10-29.md`). Once decided, the methods
  and the code repo's `msiseq/README.md` must give the rule and its
  justification, and the rule must be applied in code, not by hand. TO DO
  stubs are in both (blue note in the `ms.qmd` section "MSI/MSS status and
  high/low TMB status").

## Open: Mo (and Mini)

- Mo: update the assignment code for the assignments of 2026-10-19. The
  code is in the code repo (`Extraction_and_assignment_code/
  SigProfilerAssignment.py` and `mSigAct.R`, both on the 172-type catalog,
  run on example data), but:
  - the comments in both scripts and their README still give the old
    routing rule (TMB >= 5,000, or Ins:1:T:5+ and Del:1:T:5+ >= 90%), not
    the new one (> 4,000 indels, or poly-T ratio > 0.5 and reconstruction
    cosine similarity < 0.9 in the 83-type or 89-type fit);
  - there is no code for selecting the 731 tumors (`mSigAct.R` says "here
    only shows an example");
  - `SigProfilerAssignment.py` sets no random seed, which Mo thinks is why
    the earlier assignments could not be reproduced.
  (From `QUESTIONS_FOR_MO.md`, 2026-10-15: "Is there a list or code for the
  routing?")

- Redo the topography partial credit with the new assignments (needed for
  Tables S12 and S13, `all.indel.ID83.partial.credit.txt`). Mo sent version
  1026 on 2026-10-26. Check
  that it uses the new assignments. Mo knows the new
  assignments are approved. Steve asked Mo (2026-10-23) whether Mini has the
  partial-credit information, how it is computed, and whether we have the
  intermediate results. Waiting for the answer. Mo (2026-10-23): Mini is
  working on the topography analysis, with results in a few days. The
  partial-credit code may be `ID83_topography_analysis_code/
  add_partial_credit.R` and `generate_partialcredit_match_file.R` in the
  code repo. Then Mini reruns r.01 to
  r.04 (the raw inputs are on the Huang-lab server), and we rerun r.05
  onward (Tables S12 and S13, Figure 5, Figures S3 and S4).
- The topography code keeps indels in genes on both strands (114,297
  indels, 0.6%), although the methods say they were ignored. Mo is probably
  redoing the simulations and analysis with code that drops them, with
  Mini. Rerun the downstream topography tables and figures when that
  arrives.

## Open: Xueming

- In about a week (around 2026-10-29), once the browser is final: ask
  Xueming to tag a new GitHub release of XmWU123/indelsig (for example
  V1.1.0), so that Zenodo archives the current app under the same concept
  DOI (10.5281/zenodo.20092288), and to add a "How to cite" line (this
  paper and the DOI) to https://indelsig.net and the GitHub README.

## Done

- 2026-10-30: closed (Table S14 dropped). Steve (2026-10-25): get the inputs
  for running the ID89 topography steps
  r.01 to r.04 locally, so we can test them, including Mini's corrected
  `r.02.intersect_no_loop.sh`.
  - Done (2026-10-26): Mo sent the ID83 and ID89 partial-credit files
    (version 1026, `all.indel.ID{83,89}.partial.credit.1026.txt.gz`, kept
    outside the repos because they are controlled access).
    `find_partial_credit()` in `find_data.R` in the code repo locates them.
  - Done (2026-10-29): got `01.rtGroup.bed` and
    `repli_strand_info_new.bed`, the `-a` inputs of r.02. They are in the
    code repo under `ID83_topography_analysis_code/01_createBed/` and
    `ID83_topography_analysis_code/reference/ENCODE/`, and both r.02 scripts
    read them from there.
- 2026-10-30: Table S14 (89-type replication timing) dropped, with the
  89-type replication-timing sentences in `ms.qmd` and
  `ID89_topography_analysis_code/` and `build_sup_table_s14/` in the code
  repo. The analysis had no simulated control, so it could not separate a
  signature's trend from the general late-replication enrichment of
  indels. Tables S15 to S17 became S14 to S16 in both repos, the manifest,
  `ms_numbers.R` (`ms_num$s15` is now `ms_num$s14`), `data_flow/`, the
  project CLAUDE.md, and `explore-attributions/`. Table numbers in the Done
  items below and in older notes are the old ones. The question whether the
  ID89 partial credit was normalized (found 2026-10-29: the old r.04 matrix
  had about 3 counts per indel, so it could not come from per-indel
  fractions) no longer matters for the paper. The local ID89 run (28 GB,
  controlled access) is in
  `~/MEGA/important_mut_sig_data/id89_topography_run_2026-10-29/`.

- 2026-10-29: Figure S1 dots. The mSigAct refits added 118 dots and removed
  12 (923 to 1,029), 84 of the added from 1 to 3 tumors. 102 of the 119
  tumor-signature pairs behind them were mSigAct refits (69 MSI-H), with a
  median exposure of 196 indels, 3.6% of the tumor's indels, which passed
  the absolute 0.03 mutations/Mb cutoff because the tumors are large.
  Steve and Mo decided to count a signature in Figure S1 only if it is at
  least 5% of the tumor's indels, as in Figure 4A (code repo 9a8a2bc, paper
  896b920, legend updated). Figure 4B is left as is. With the 5% rule it
  keeps the quoted numbers and the MMR cluster (InsDel7, InsDel_J, InsDel2c,
  InsDel2b, InsDel_D, C_ID2 with SBS44, SBS6, SBS26), but InsDel11 and
  InsDel_M drop out, and its hand-placed boxes would need adjusting.
  Precedent for a share cutoff: deconstructSigs drops signatures under 6%
  (Rosenthal et al. 2016, https://doi.org/10.1186/s13059-016-0893-4).

- 2026-10-29: Figure S1, stomach InsDel1a rose more than tenfold (0.13 to
  1.60 mean mutations/Mb) with the new assignments. Almost all of the rise
  is one tumor, CPCT02450014T (267,486 indels, refit with mSigAct). On the
  same double-counted scale, InsDel1a went from 0 to 116,852 and InsDel_D
  fell from 126,822 to 15,520 (InsDel_L 368,017 to 335,958, InsDel7 40,134
  to 58,926). The 14 tumors are not quite the same: CPCT02130122T (1,230
  indels) lost InsDel1a and CPCT02450014T gained it. The other 13 have about
  the same InsDel1a as before. InsDel_D and InsDel1a are similar (cosine
  0.85, both mostly T insertions into runs of 5 or more), but the tumor's T
  insertion peaks are A-flanked, as in InsDel1a, and InsDel_L covers its
  Ins C(4,6) and L(2+)R(5+) peaks. The new fit reconstructs the spectrum
  better (cosine 0.9992 vs. 0.9919). Steve decided to keep the new fit.
  Stomach C_ID1 (0.09 to 0.70) probably has the same cause (83-type ID_D
  fell as C_ID1 rose), not checked per tumor.

- 2026-10-29: reran the outputs that read the assignment tables after
  CPCT02010349T was changed on 2026-10-26, as Mo advised, to the exposures
  from the SPA rerun of 2026-10-25 (InsDel5a / C_ID5 = 288, was InsDel19a /
  C_ID19 = 288; data repo 706124e). It has 144 indels and a poly-T ratio of
  0.097, so it is not an mSigAct tumor, and it is not in the Figure S1 list.
  Mo says (2026-10-26) the partial-credit files of 2026-10-26 already have
  the new assignment. Rerun: Tables S9, S10, S11, S15, S16, and S17,
  `build_fig4/`, `build_fig_s1/`, `build_fig_s6/` (with the two
  prototypes), and `ms_numbers.R` (code repo e6de415). `vignette/` reads
  only the file's modification time, so it was not rerun. Tables S9, S10,
  and S17 and Figures 4, S1, and S6 are unchanged in content. Table S11
  changed only in CPCT02010349T. In Tables S15 and S16 the Bone/Soft tissue
  C_ID19, C_ID5, InsDel5a, and InsDel19a rows, the pooled C_ID5, InsDel5a,
  C_ID19, and InsDel19a rows, and some FDRs changed slightly. No call at
  q < 0.2 changed. Tables S11, S15, and S16 were copied to `sup_tables/`.

- 2026-10-24: Steve decided on the Figure S1 dot cutoff: a dot is hidden
  when <= 10% of tumors have the signature and the mean is < 0.03
  mutations/Mb, as in `build_fig_s1/make_figure_s1_layout.R`. The deleted
  notebook used 1%. 10% reproduces the earlier Illustrator figure.

- 2026-10-23: `QUESTIONS_FOR_MO.md` (compiled 2026-10-15: Table S17
  `mutation_burden` and poly-T ratio, mSigAct routing, the hypermutation
  definition, and MSI status) was merged into this list and deleted (it is
  in the git history). Everything in it was answered and is recorded
  below, except updating the routing code ("Open: Mo") and stating MSI status
  consistently ("Open: Steve").

- 2026-10-22: meeting with Mo (8 PM US time).

- 2026-10-22: Figures S1a and S1b are built by
  `build_fig_s1/make_figure_s1_layout.R` in the code repo, replacing the
  Illustrator assembly, with etiologies from
  `build_fig_s1/figure_s1_etiologies.tsv`. The color is the mean (not the
  median) mutations per Mb among tumors with the signature, as the code
  always computed. The legend and the caption in `sup_figs/sup_figs.qmd` now
  say so.

- 2026-10-22: new signature assignments adopted. Mo (email "Question about
  S11", 2026-10-19) could not reproduce some of the earlier assignments
  (probably because no random seed was set) and reconstructed the original
  workflow: SigProfilerAssignment on the 172-type catalogs, then mSigAct
  reruns for (1) tumors with > 4,000 indels and (2) tumors with poly-T
  ratio > 0.5 and a SigProfilerAssignment reconstruction cosine similarity
  < 0.9 in either the 83-type or the 89-type fit (731 tumors).
  - The four attachments of that email are in
    `unified_indels_data/assignments_2026-10-19/`:
    `ID83.SPA.assignment.1020.txt`, `ID89.SPA.assignment.1020.txt`,
    `assignment.high.similarity.samples.for.figureS1.txt` (6,516 tumors),
    and `mSigAct.reanalyzed.samples.txt` (731 tumors). The first two were
    removed on 2026-10-25 because they are identical to
    `liu_et_al_83_assignment.tsv` and `liu_et_al_89_assignment.tsv`. 819 tumors have
    exposures that differ from the earlier ones by more than 0.5 (Mo said
    784).
  - The list has all 248 assigned tumors with > 4,000 clipped indels. The
    one other such tumor, the exome tumor
    `WES::TCGA-EW-A2FV-01A-11D-A17D-09` (4,198 indels), is one of the 3
    tumors without assignments, because it is the only tumor with its
    signature, which is its own spectrum (Steve).
  - `mSigAct.reanalyzed.samples.txt` in this repo holds the 731-tumor list.
    The old 2026-10-18 list (1,535 tumors) is only in the git history.
  - The new tables replace `liu_et_al_{83,89}_assignment.tsv` in the data
    repo, whose README describes them. Figure 4, Figure S6, and Tables S9,
    S10, S11, S15, and S16 were rebuilt and copied here. Table S17 did not
    change.
  - The methods give the new routing rule (731 tumors refit with mSigAct,
    including 200 of the 202 MSI-H tumors), with the counts from
    `ms_numbers.R`. The blue Table S11 note was replaced: every column sum
    is now twice the indel count, which resolves the question of the 503
    tumors with single-counted sums.
  - The sex-bias paragraph follows the new Table S15 (lymphoid and "Other"
    hits gone, Pancreas InsDel12/C_ID12 female-enriched, Liver C_ID5
    added).
- 2026-10-22: age analysis adjusted for cancer type. The pooled "Overall"
  Spearman correlation in Table S16 and Figure S6B confounded age with
  cancer type (InsDel5b was significant overall but had no age effect
  within cancer types). It was replaced by regression models with cancer
  type as a covariate (`age_glm_table()` in
  `build_fig_s6/build_fig_s6_stats.R`, Table S16 sheet "Adjusted for cancer
  type"): negative binomial exposure ~ age + cancer type in tumors with the
  signature, and logistic presence ~ age + cancer type, with and without
  log10 total indels. The InsDel5b sentences were deleted, and a green
  sentence in the age paragraph gives the adjusted effects for C_ID1,
  InsDel1a, C_ID5, and InsDel5a. Prototypes of the joint age/sex models
  and of the adjusted Table S15 are in `build_fig_s6/` (see "Open: Steve
  with Mo").
- 2026-10-22: the methods give the high-TMB cutoff (Mo, 2026-10-18: > 4,000
  indels after clipping to R <= 9, and the extraction was redone after
  clipping) in the de novo extraction section (249 of 6,978 genomes, from
  `ms_num$tmb`) and in the MSI/TMB section.
- 2026-10-22: Table 1's first column is "89/476-Type Classification
  Signature ID" in the code repo's `vignette/table_1_autogenerated.xlsx`,
  its copy here, and the hand-edited `table_1_autogenerated_edited.xlsx`
  that `ms.qmd` reads.
- 2026-10-22: the vignette, the single-signature pages and Table 1 were
  rerendered with the new input data. Nothing substantive changed: the
  vignette does not read the assignments.
- 2026-10-22: `ms.qmd` cites the Degasperi et al. Zenodo data record
  (10.5281/zenodo.5571551, `SBS_v2.03.zip`) for the SBS activities in the
  Figure 4B methods, and "(Table S11)" now refers to our indel activities.
- 2026-10-22: the browser DOI was checked. 10.5281/zenodo.20092288 is the
  concept DOI (always the latest version) of the Indel Signature Browser
  by Xueming Wu, Steven G. Rozen, and Mo Liu. Its only version is V1.0.0
  (2026-05-09), archiving GitHub release V1.0.0 of
  XmWU123/indel-signature-browser, since renamed XmWU123/indelsig (GitHub
  redirects; XmWU123 is the GitHub name). The repo was pushed after that
  (2026-07-13) with no new tag, so the archive is older than the running
  app. `ms.qmd` now cites the DOI after https://indelsig.net. indelsig.net
  and the GitHub README do not give the DOI or the paper (see "Open:
  Xueming").
- 2026-10-22: Table S9 no longer has `MSI_status` ("MSI status" in the main
  sheet, "MSI status = MSI" in the summary sheet), following the decision
  to use only MSIseq. `build_sup_table_s9/build_sup_table_s9_pol_proofreading.R`
  in the code repo no longer writes them, and its comments and README say
  not to add them back. It also now spells "tumor". "Total somatic indels"
  is the column sum of the clipped 89-type spectra.
- 2026-10-18: the same MSI-H and hypermutation definitions are used
  everywhere. The scripts in `explore-attributions/` were switched on
  2026-10-17 (ce9a4f3). On 2026-10-18 the numbers in
  `id1_id2_idN_findings.md` and `more-on-attributions.md` were updated from
  the regenerated outputs (several findings about MSS hypermutated tumors
  no longer hold), and the two scripts that reproduce Mo's original figure
  with `MSI_status` say so in a comment.
- 2026-10-17: `build_fig4/all_cancer_types_enrichment.R` calls a tumor
  hypermutated if it has > 4,000 indels in the R <= 9 clipped spectra
  (column sums of `liu_et_al_83_spectra.tsv`). Figure 4A was recomputed.
  MSS-hyper now has 48 tumors instead of 723, and its only dot is InsDel3b
  instead of InsDel14. The legend already says "> 4,000 indels", and the
  results text does not mention the MSS-hyper column.
- 2026-10-17: Table S17 is built by
  `build_sup_table_s17/build_sup_table_s17.R` in the code repo.
  `mutation_burden` and the poly-T ratio are computed from the clipped
  spectra, and the other columns reproduce exactly. The legend gives the
  channels in the poly-T ratio (INS:T:1:5+ and DEL:T:1:5+, so deletions
  are from runs of 6 or more Ts) and no longer says that tumors with a
  ratio of 0.5 or more went to mSigAct. The methods wording about poly-T
  tracts was fixed on 2026-10-18.
- 2026-10-17: answered by Mo and done: ID_N and InsDel_N_beta in Table S15
  (Table S15 and Figure S6 regenerated), and `donor.all_projects.tsv.gz`
  (the ICGC portal download of 9 May 2024 named in the data repo README).

# Figures, sup figures, and sup tables vs. the code repo (2026-10-07)

Run `Rscript check_main_figures.R`, `Rscript check_sup_figures.R`, and
`Rscript check_sup_tables.R` to see the current state. Open items from this
list are in the sections above.

1. Closed (2026-10-17): read Mo's email replies. Follow-ups are in the
   sections above.
2. Topography annotations and simulations. Answered (Mo, 2026-10-07 and
   2026-10-08, and Steve, 2026-10-17): `dna.region` and `trans.strand` came
   from ICAMS with `trans.ranges.GRCh37` (GENCODE 30 CCDS genes), genome
   build hg19 (GRCh37). SigProfilerSimulator 1.2.2 was run as
   `sigSim.SigProfilerSimulator("project", "/path/to/project/", "GRCh37",
   contexts=["ID"], simulations=10)`. The strand column (+, -, Q) is
   computed by SigProfilerSimulator (`mutational_simulator.py`, lines
   1134-1151). Decision: 10 simulations per tumor is enough. Mo asked that
   we cite Otlu et al. 2023 and Bergstrom et al. 2020. Done (2026-10-17):
   filled the placeholders in the green "Simulating synthetic cancer
   datasets" and "Annotating somatic indels as genic or intergenic"
   paragraphs in `ms.qmd`, and updated the blue strand note with Mo's
   answer. Still open: indels in genes on both strands (see "Open: Mo").
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
5. Closed (2026-10-17): Table S17 (sample metadata). Patient, age, sex,
   MSI status, major cancer type, and cohort come from
   `unified_indels_data/sample_info.tsv`, which also has 41 samples not in
   S17. `MSIseq` matches `MSIseq_MSI.H` in `sample_info.tsv` for all 6,975
   tumors. `Cancer Type` matches `Original_Cancer_Type` there, except that
   Colon, Head, and Bone-SoftTissue are renamed Colon/Rectum, Head and
   neck, and Bone/Soft tissue (2026-10-07). S17 has a "legend" sheet.
   Decision (Mo, 2026-10-16): mutation counts come from the spectra in
   `unified_indels_data`, which are clipped to repeat length R <= 9 (column
   sums of `liu_et_al_83_spectra.tsv`). The earlier `mutation_burden` and
   poly-T ratio came from an older, pre-cap-9 version of the data.
6. Closed (2026-10-18): Figure 4A MSS-hyper threshold (see "Done").
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
8. Review the STAR Methods for the topography analysis (see "Open: Steve").
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
