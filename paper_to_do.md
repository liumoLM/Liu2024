# Journal requirements (cell resources)

- 7 display items allowed (We will have 1 table and 6 main figures)
- Word Count: The main text should be under 7,000 words. This count includes the figure legends but excludes the STAR Methods, supplemental legends, and references.
- Supplemental Items: You may include additional display items as supplemental information at the editor's discretion. For Cell Research Articles, this is often up to 7 supplemental figures.
- Mandatory Items: You must also include a Graphical Abstract (exactly 1,200 x 1,200 pixels) and 3–4 Highlights (each 85 characters or fewer). [1, 2, 3, 4]
http://8.138.15.71/


# To do

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
   scripts are in `ID89_topography_analysis_code/`. Still missing: the r.04
   inputs `tissue.rtGroup.{Mutation,Simulated}_50cutoff.xlsx` (Mini said
   r.09 uses them but did not attach them), and the 89-type per-indel input
   `all.indel.ID89.partial.credit.txt` (too big to email; Mini asks whether
   we already have it). Methods to fix: the replication-timing minimum is
   500 real indels, not 1,000, and "None" also covers cancer types where
   neither the real nor the simulated slope is significant. "Flat" means
   both slopes are significant in the same direction and the decile term in
   the third model is not.
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
   pre-cap-9 version of the data. With the current spectra, 44 tumors would
   be routed differently between SigProfilerAssignment and mSigAct (methods,
   "Table S17" paragraph). `MSIseq` matches `MSIseq_MSI.H` in
   `sample_info.tsv` for all 6,975 tumors. `Cancer Type` matches
   `Original_Cancer_Type` there, except that Colon, Head, and
   Bone-SoftTissue are renamed Colon/Rectum, Head and neck, and Bone/Soft
   tissue (2026-10-07). S17 now has a "legend" sheet. Check its description
   of `mutation_burden` ("Number of somatic indels in the tumor") once the
   source of that column is known. Overall, S17 needs a lot of work.
   Decide whether we need it (probably yes, as an overall description of
   the data set), and trace where each column's numbers come from,
   ideally so that S17 can be built by a script in the code repo.
7. Review the STAR Methods for the topography analysis.
8. Check the Otlu et al. code (SigProfilerTopography,
   https://github.com/AlexandrovLab/SigProfilerTopography) to see whether
   they normalize indels to the pyrimidine strand, for example GAGAGA->GAGA
   to TCTCTC->TCTC, and compare with what we do.
9. Later: Fig 6 panels A to C are built from `build_fig6/ID4_IDF_signature.pdf`
   and `build_fig6/sequence_logos.pdf` (checked 2026-10-07) and carry the
   same information, but they are not visually identical to the source
   plots. For example, the panel C logos are rescaled to 0-2 bits and the
   DNA diagrams are drawn by hand. Come back to this.

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
