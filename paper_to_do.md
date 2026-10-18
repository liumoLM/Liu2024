# Journal requirements (cell resources)

- 7 display items allowed (We will have 1 table and 6 main figures)
- Word Count: The main text should be under 7,000 words. This count includes the figure legends but excludes the STAR Methods, supplemental legends, and references.
- Supplemental Items: You may include additional display items as supplemental information at the editor's discretion. For Cell Research Articles, this is often up to 7 supplemental figures.
- Mandatory Items: You must also include a Graphical Abstract (exactly 1,200 x 1,200 pixels) and 3–4 Highlights (each 85 characters or fewer). [1, 2, 3, 4]
http://8.138.15.71/


# To do

## Tomorrow (2026-10-17)

- Done (2026-10-17): Item 6, `build_fig4/all_cancer_types_enrichment.R` now
  calls a tumor hypermutated if it has > 4,000 indels in the R <= 9 clipped
  spectra (column sums of `liu_et_al_83_spectra.tsv`). Figure 4A was
  recomputed and copied to `main_figures/`. MSS-hyper now has 48 tumors
  instead of 723, and its only dot is InsDel3b instead of InsDel14. The
  legend already says "> 4,000 indels", and the results text does not
  mention the MSS-hyper column.
- Item 6: use the same definition everywhere else in the paper, and in
  `explore-attributions/`.
- Done (2026-10-17): Item 5, Table S17 is now built by
  `build_sup_table_s17/build_sup_table_s17.R` in the code repo.
  `mutation_burden` and the poly-T ratio are computed from the clipped
  spectra, and the other columns reproduce exactly. The legend now gives the
  channels in the poly-T ratio (INS:T:1:5+ and DEL:T:1:5+, so deletions
  are from runs of 6 or more Ts) and no longer says that tumors with a ratio
  of 0.5 or more went to mSigAct. The methods wording about poly-T tracts
  was fixed on 2026-10-18.
- Routing between SigProfilerAssignment and mSigAct (item 5). Mo replied
  2026-10-18 and sent `mSigAct.reanalyzed.samples.txt`: 1,535 tumors with
  `mutation_burden` >= 4,000 or poly-T ratio > 0.5 in Table S17 (clipped
  counts). The list matches that rule applied to the current S17 exactly
  (all 202 MSI-H tumors are on it). Done (2026-10-18): the methods now give
  this rule. Still open:
  - The list matches the current S17 (built 2026-10-17), not the earlier
    S17, where the same rule gives 1,517 tumors (39 on the list do not meet
    it, 21 that meet it are not on the list). Ask Mo whether the list was
    regenerated from the current S17 or is the set the fits were actually
    run on.
  - Table S11 sums. Mo asked for the names of the 503 tumors with
    single-counted sums. They are in `s11_single_counted_503.tsv` (447 on
    the mSigAct list, 56 not, and 1,088 tumors on the list have
    double-counted sums). Steve emailed Mo the list (2026-10-18), asking
    why these 503 are not doubled. Waiting for the reply.
- Review Mo's replies and act on them. Still open with Mo:
  - MSIseq: Mo replied 2026-10-18. MSIseq was run on all tumors, the
    `MSIseq` column is MSIseq's own output (the cohort MSI calls are a
    subset of it), and the 4,000 threshold is after clipping. Mini ran
    MSIseq. Waiting for Mini's details (inputs, VCFs and filters, repeats
    annotation, sequence length) and code. Done: the methods no longer say
    that MSIseq was run only on high-TMB tumors or how the 4,000 cutoff
    was chosen.
  - Mo's partial-credit table for Table S14 (`all.indel.ID89.partial.credit.txt`),
    promised for the weekend.
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
   project CLAUDE.md records it. Still to do: use the same definition
   everywhere else in the paper and in `explore-attributions/`, which
   still uses a > 5000 rule on Table S11 sums.
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

# From here down old


# To do steve

- Look at fig 5 (? indel vs SBS associations) and write up

- Look at Fig 6 and corresponding section of the MS.

- Consider ID14 trans. strand. bias and interaction w/ SBS93 --> potential etiology of id14

- consider info on assignments below and review that section of the ms

- H -- how common are the 2 deletions that make up H across all tumors? how constant is the ratio in skin tumors? Scatter plot?


# Decisions

# Status: 


# Under consideration / nice to have

Electronically (via link) Examples of mutations in the different (somewhat obscure) classification

The Sig F / Reijn annotated VCFs need to be redone with "cap 9"

Review SBS / indel associations

# Data sources for Degaspari SBS data
  https://zenodo.org/records/5571551/files/SBS_v2.03.zip?download=1                                                                                                                        
                                                                                                                                                                                            
  (From Zenodo record https://zenodo.org/records/5571551)
  https://doi.org/10.5281/zenodo.5571551 

The following is unreliable:
# Assignments

- 1a everywhere
- 1b a bit enriched in breast
- 1c everywhere 
- ID* much rarer than in COSMIC -- to partial credit in 476 for the Del T 6+ peak
- 2a some enrichment in breask lunk eso, colon
- 2b MSI per Mo, high in uterus, breast,colon
- 2c MSI per Mo,  high prosate, colon, and no msi lymphoid
- 3a prevalent in bladder
- 3b liver and lung
- 4
- 5
- 6
- 7 MSI per Mo
- 8
- 9
- 10
- 11 super rare
- 12 prostate
- 13 skin (expected)
- 14 colon and esophagus (and stomach)
- 17 ** breast, low counts but top2A **** worth revisting in light of our paper / stomach
- 18 colon (colibactin)
- 19a rare -- assoc w MSI in colon and esophagus
- 19b everywhere
- 23 rare; but should not be in breast or prostate
- A alpha beast / colon, but not MSI
- A beta strong breast
- B everywhere except stomach
- C super rare -- some very high in Colon, Breast
- D MSI per Mo, esp prostate and colon
- E generally rare..
- F rare; we think we have an explanation
- G rare
- H ** skin only, consists of basically of Del(6,9):U(4,):R(2,9) and Del(7,9):M(6,)
- I super rare - 10 tumors
- J MSI per Mo in prostate, breask, colon
- K super rare; 9 tumors nly pancrease, prostate, colon
- L super rare; about 20 tumors only CNS, SKin, pancrease
- M super rare, high MIS assoc, only lung breast colon protate
- N alpha, somwhat rare, slightly more in stomach?
- N beta, everywhere except thymus kindey liver *** look for complimentary signature 
