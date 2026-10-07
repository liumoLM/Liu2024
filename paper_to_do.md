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

1. Waiting on Mini (email sent 2026-10-07) for the missing replication-timing code
   and data. No script in the code repo creates `09_barplot_data.xlsx` and
   `09_all_lm_trend.xlsx` (inputs to Fig 5D and Fig S4, and probably the
   source of Table S13), and there is no code for the 89-type
   replication-timing analysis (Table S14). Also asked for the
   `tissue.rtGroup.{Mutation,Simulated}_50cutoff.xlsx` matrices from r.04.
   When they arrive, add them to the code repo and to the manifests.
2. Waiting on Mo (email sent 2026-10-07) about Fig S5. The 83-type top row
   is stale. Its category counts (for example
   mouse single-T deletions 370) match `build_fig_s5/F83.pdf` from before
   e6f1c34 (2026-08-14), which switched to
   `mSigSpectra::annot_vcf_to_83_catalog(clip_le_9 = TRUE)`. The current
   F83.pdf (297) agrees with the 89-type row. The 89-type and 476-type rows
   match the current plots, which changed only in layout on 2026-08-15.
   Asked Mo for any code that assembles Fig S5
   from its constituents, which have different aspect ratios, and which code
   made the top-row panels.
3. Waiting on Mo (email sent 2026-10-07) for the code that makes the
   constituent plots of Fig S2 (83-type, 89-type, and 476-type zoom panels
   for C_ID13/ID_C, C_ID11/C_ID16, and C_ID9/ID_A), and any notes on how the
   panels were assembled. Put it in a new `build_fig_s2/` in the code repo
   and add it to `sup_figures_manifest.tsv`.
4. Table S17 (sample metadata) is only partly reproducible. Patient, age,
   sex, MSI status, major cancer type, and cohort come from
   `unified_indels_data/sample_info.tsv`, which also has 41 samples not in
   S17. `mutation_burden` and the polyT ratio do not match the current
   spectra (mutation_burden equals the spectrum total for only 486 of 6,975
   tumors and is otherwise larger). They probably come from an earlier,
   pre-cap-9 version of the data. With the current spectra, 44 tumors would
   be routed differently between SigProfilerAssignment and mSigAct (methods,
   "Table S17" paragraph). `MSIseq` and `Cancer Type` have no known source.
5. Replication strand text in `ms.qmd`: "12 out of the 26 signatures with
   sufficient data". The 12 matches the Fig 5 inputs, but the 26 still needs
   to be derived from the r.06 outputs.
6. Review the STAR Methods for the topography analysis.
7. Cite Otlu et al. 2023, "Topography of mutational signatures in human
   cancer", Cell Reports,
   https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00941-5,
   DOI 10.1016/j.celrep.2023.112930 (https://doi.org/10.1016/j.celrep.2023.112930).
8. Check the Otlu et al. code (SigProfilerTopography,
   https://github.com/AlexandrovLab/SigProfilerTopography) to see whether
   they normalize indels to the pyrimidine strand, for example GAGAGA->GAGA
   to TCTCTC->TCTC, and compare with what we do.
9. In the topography methods, cite SigProfilerSimulator for the simulated
   genomes: Bergstrom et al. 2020, "Generating realistic null hypothesis of
   cancer mutational landscapes using SigProfilerSimulator", BMC
   Bioinformatics, https://pubmed.ncbi.nlm.nih.gov/33028213/,
   DOI 10.1186/s12859-020-03772-3 (https://doi.org/10.1186/s12859-020-03772-3).
10. In the topography methods, state that the simulated genomes are
    available on request, or can be deposited in Zenodo on request.
11. Later: Fig 6 panels A to C are built from `build_fig6/ID4_IDF_signature.pdf`
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
