# Journal requirements (cell resources)

- 7 display items allowed (We will have 1 table and 6 main figures)
- Word Count: The main text should be under 7,000 words. This count includes the figure legends but excludes the STAR Methods, supplemental legends, and references.
- Supplemental Items: You may include additional display items as supplemental information at the editor's discretion. For Cell Research Articles, this is often up to 7 supplemental figures.
- Mandatory Items: You must also include a Graphical Abstract (exactly 1,200 x 1,200 pixels) and 3–4 Highlights (each 85 characters or fewer). [1, 2, 3, 4]
http://8.138.15.71/


# To do

## Tomorrow (2026-10-17)

- Done (2026-10-17): Item 7, `build_fig4/all_cancer_types_enrichment.R` now
  calls a tumor hypermutated if it has > 4,000 indels in the R <= 9 clipped
  spectra (column sums of `liu_et_al_83_spectra.tsv`). Figure 4A was
  recomputed and copied to `main_figures/`. MSS-hyper now has 48 tumors
  instead of 723, and its only dot is InsDel3b instead of InsDel14. The
  legend already says "> 4,000 indels", and the results text does not
  mention the MSS-hyper column.
- Item 7: use the same definition everywhere else in the paper, and in
  `explore-attributions/`.
- Item 6: recompute Table S17 `mutation_burden` and the poly-T ratio from
  the clipped spectra, ideally with a build script in the code repo.
- Item 6: decide what to do about the 43 tumors whose SigProfilerAssignment
  vs. mSigAct routing changes with the clipped counts, and whether the
  routing threshold should be >= 5,000 (as in the methods) or 4,000.

Steve review Mini's methods -

Steve Review subsection "Demographic Associations with Gender and Age"


check Shiny app persistent doi: 10.5281/zenodo.20092288
XmWU123 Github name


Steve: check "self" referencing from web page to shiny app doi

steve check: Sup Fig 1 legend

# To do: figures, sup figures, and sup tables vs. the code repo (2026-10-07)

Run `Rscript check_main_figures.R`, `Rscript check_sup_figures.R`, and
`Rscript check_sup_tables.R` to see the current state.

1. Read Mo's email replies.
2. Ask Mo to confirm that the `dna.region` (G or I) and `trans.strand`
   columns in the annotated indel file read by r.01_createBed.R came from
   ICAMS, with which version and which transcript ranges
   (`trans.ranges.GRCh37`, GENCODE 30 CCDS genes, or GRCh38). ICAMS sets
   `dna.region` to G when `trans.strand` is + or -, else I
   (`CreateOneColIDMatrix()`, `R/ID_functions.R`). The green "Annotating
   somatic indels as genic or intergenic" methods paragraph in `ms.qmd`
   assumes this and marks the genome build as unconfirmed. Also ask for
   the genome build and SigProfilerSimulator settings, to fill the placeholders in the green "Simulating synthetic
   cancer datasets" paragraph (Mo already gave 10 simulations per tumor and
   version 1.2.2). Both questions are in the draft email "Topography:
   indels in genes on both strands". Also decide whether 10 simulations
   per tumor, rather than the 100 Otlu et al. recommend, is acceptable (Mo
   asked).
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
4. Waiting on Mo (email sent 2026-10-07) about Fig S5. Mo has no
   assembly code, so we wrote `build_fig_s5/make_figure_s5_layout.R` in
   the code repo, which builds `figure_s5_draft.pdf` in the paper's layout
   from the current catalogs. The paper's `sup_figs/figure_s5.pdf` has a
   stale 83-type top row (e.g. 271 single-T deletions and 253 single-T
   insertions in the RNase H2 null cells, now 99 and 129). Sent Mo both
   PDFs and asked Mo either to update the paper version with the correct
   top row from the draft, or to edit the draft as needed. Replace
   `sup_figs/figure_s5.pdf` with the result.
5. Closed (2026-10-07): Fig S2. Mo has no assembly code, so we wrote
   `build_fig_s2/make_figure_s2_layout.R` in the code repo. It builds
   `figure_s2_draft.pdf` in the paper's layout from the Table S2, S4, and
   S6 signatures. The paper's `sup_figs/figure_s2.pdf` is edited by hand
   from the draft, as recorded in `sup_figures_manifest.tsv`.
6. Table S17 (sample metadata) is only partly reproducible. Patient, age,
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
   sums of `liu_et_al_83_spectra.tsv`). So recompute `mutation_burden` and
   the poly-T ratio from those spectra, ideally with a build script in the
   code repo. With the current spectra, 43 tumors would be routed
   differently between SigProfilerAssignment and mSigAct (recount of
   2026-10-15). Decide whether the routing should be redone with the
   clipped counts, or whether the methods should say which counts the
   routing used. The routing threshold in the methods is >= 5,000 indels,
   not 4,000. Check whether that is intended. Background is in
   `QUESTIONS_FOR_MO.md`.
7. Figure 4A MSS-hyper threshold. Decision (Mo, 2026-10-16): a tumor is
   hypermutated if it has > 4,000 indels, counted from the R <= 9 clipped
   spectra in `unified_indels_data`. `build_fig4/all_cancer_types_enrichment.R`
   instead uses signature exposure sums > 4000. Those sums count each indel
   twice for tumors fit with SigProfilerAssignment, so its effective cutoff
   is about 2,000 indels for most tumors, giving 723 MSS-hyper tumors
   against 48 with the decided rule (2026-10-15 count). To do: change the
   script to use the spectrum totals, recompute Figure 4A, update its
   legend and text, and use this definition of hypermutated throughout
   the paper (also `explore-attributions/` and the project CLAUDE.md, which
   still describe a > 5000 rule on Table S11 sums).
8. Figures S1 and S6 (2026-10-16). Mo added `build_fig_s6/build_fig_s6_stats.R`
   and the Figure S1 metadata scripts.
   - Figure S1 reproduces. It now reads `sample_info.tsv` and the assignment
     tables from the data repo, and the copies in `build_fig_s1/` were
     deleted (identical plot data). `build_all_meta.R` moved to the data
     repo as a record of how `sample_info.tsv` was assembled.
   - `HMF_metadata.tsv` (controlled-access Hartwig clinical data) was
     removed from the code repo and its git history, and the repo was
     force-pushed. Mo, Mini, and Xueming were emailed to reset their clones.
   - Table S16 reproduces from the script. Table S15 does not: the script
     adds 26 ID_N and InsDel_N_beta rows, which changes the BH FDRs by up
     to 16%, but the same 10 groups pass FDR < 0.2.
   - Emailed Mo (2026-10-16) why ID_N and InsDel_N_beta were left out of
     Table S15, and where `HMF_metadata.tsv` and the PCAWG donor table came
     from. Details in `sup_figs_and_tables_sources.md`, items 4 and 5.
9. Review the STAR Methods for the topography analysis.
10. Check the Otlu et al. code (SigProfilerTopography,
   https://github.com/AlexandrovLab/SigProfilerTopography) to see whether
   they normalize indels to the pyrimidine strand, for example GAGAGA->GAGA
   to TCTCTC->TCTC, and compare with what we do.
11. Later: Fig 6 panels A to C are built from `build_fig6/ID4_IDF_signature.pdf`
   and `build_fig6/sequence_logos.pdf` (checked 2026-10-07) and carry the
   same information, but they are not visually identical to the source
   plots. For example, the panel C logos are rescaled to 0-2 bits and the
   DNA diagrams are drawn by hand. Come back to this.
12. Review the naming of the figure-generating code and folders. On
   2026-10-17 the supplementary figures were renumbered to S1 to S10 (old
   S2A to S2C are now S2 to S4, old S3A to S3C are now S5 to S7, and old S4
   to S6 are now S8 to S10). The code repo folders (`build_fig_s2/` to
   `build_fig_s6/`) and file names (`FigureS3_*.pdf`, `FigureS4.pdf`,
   `figure_s5_draft.pdf`) still use the old numbers, as do the mixed-style
   file names in `sup_figs/` (`Figure S5.pdf`, `figure_s9.pdf`,
   `figure_s2_to_s4.pdf`). `sup_figures_manifest.tsv` maps them correctly.
   Consider names that describe content rather than figure numbers, so
   that future renumbering needs no renames.

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
