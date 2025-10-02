# Koh 2025 summaries

## Claude Sonnet 4.5

# Summary of "A redefined InDel taxonomy provides insights into mutational signatures"

## Abstract & Introduction

- **Study focus**: Small insertions and deletions (InDels) have received less attention than substitutions despite their harmful effects
- **Approach**: Created CRISPR-edited human cell models with postreplicative repair dysfunction (PRRd), including DNA mismatch repair (MMR) and replicative polymerase defects
- **Key finding**: Current InDel classification (COSMIC-83) cannot discriminate different InDel signatures from each other or background
- **Solution**: Developed new 89-channel classification system incorporating flanking sequences and longer homopolymers
- **Results**: Discovered 37 InDel signatures (27 new) in seven tumor types from 100,000 Genomes Project
- **Application**: Created PRRDetect classifier for identifying PRRd status in tumors with immunotherapy implications

## Background

- **InDel importance**: Second most common genetic variation after substitutions, reflecting underlying mutational processes
- **Previous limitations**: Most mutational signature research focused on substitutions; only 18 small InDel signatures previously identified
- **Clinical relevance**: Accurate InDel characterization crucial for detecting homologous recombination deficiency and microsatellite instability (MSI)
- **Study goals**: Establish "ground truth" experimental InDel signatures and improve classification to uncover new etiologies

## Results - Diversity of InDel patterns in PRRd

- **Cell models**: Generated 10 CRISPR-edited RPE1 cell lines with MMR knockouts and polymerase mutations
- **Mutation burden**: All gene edits except ΔSETD2 showed elevated InDel burdens (2-300 fold increases)
- **Signature diversity**: Each genotype showed unique InDel profiles with variations between MMRd and polymerase dysfunction
- **Strand bias**: POLE mutants showed leading strand enrichment, POLD1 mutants showed lagging strand bias for specific InDels
- **Key observation**: Despite signature diversity, COSMIC-83 taxonomy showed high similarity (>0.9 cosine similarity) between different genotypes

## Results - Limitations of current InDel taxonomy

- **Misclassification**: MMRd signatures (ΔMSH2, ΔMLH1) resembled normal replication error signatures (ID1, ID2) rather than MMRd-associated ID7
- **Root cause**: COSMIC-83 aggregates InDels at homopolymers >5bp into single channels, losing discriminatory power
- **Missing signal**: ID7 lacks signal in the most informative homopolymer channels (>5bp)
- **Poor discrimination**: Polymerase mutant signatures were indistinguishable from each other and from ID1
- **Contrast with substitutions**: PRRd substitution signatures show clear, distinct patterns unlike corresponding InDel signatures

## Results - New InDel classification framework

- **Design principles**: Incorporated flanking sequence context (5' and 3'), expanded homopolymer channels, and genome-wide motif prevalence
- **Channel development**: Started with 476 channels, consolidated to 89 based on signal distribution across 18,522 tumors
- **Key improvements**: Expanded 1bp A/T InDel channels, condensed rare/uninformative channels
- **Improved discrimination**: Mean cosine similarity to background dropped from 0.89 to 0.68 (p=1.917×10⁻⁷)
- **Better separation**: Gene-edit signatures more distinguishable from each other (mean similarity 0.57 vs 0.64)

## Results - Mechanistic insights from new taxonomy

- **MMRd patterns**: Deletions amplified at longer homopolymers (8-9bp > 5-7bp > 0-4bp)
- **Polymerase patterns**: Insertions elevated at shorter homopolymers (5-7bp > 8-9bp > 0-4bp)
- **Biological explanation**: Pattern reflects 5-7bp "footprint" where polymerases contact duplex DNA
- **De novo extraction**: 89-channel format yielded 4 signatures vs only 2 with COSMIC-83
- **Algorithm comparison**: Three different extraction algorithms consistently extracted more signatures with 89-channel format

## Results - New InDel signatures in seven cancer types

- **Dataset**: Analyzed 4,775 tumors from seven high-TMB cancer types in 100,000 Genomes Project
- **Discovery**: Identified 37 consensus InDel signatures (InDS); 10 matched known signatures, 27 were new
- **Exogenous sources** (5 signatures): Tobacco (InD3a/b), UV (InD13), colibactin (InD18), platinum (InD32)
- **Endogenous sources** (20 signatures): Normal replication errors (InD1, InD2a/b), TOP1 (InD4a), HR deficiency (InD6), NHEJ (InD8), APOBEC (InD9a/b/c)
- **PRRd-specific** (8 signatures): MMRd (InD7, InD19), POLE defects (InD14, InD15), POLD1 defects (InD14, InD21), combined defects (InD16a/b, InD20)

## Results - APOBEC mechanism

- **Pattern**: InD9a shows 1bp C deletions at TCT and TCA motifs at short poly-T tracts
- **Proposed mechanism**: APOBEC deaminates C to U, UNG removes uracil creating abasic site, template strand slippage occurs, resulting in C deletion
- **Supporting evidence**: Corroborated by APOBEC overexpression experiments in DT40 cells
- **Related signatures**: InD9b/c show similar C deletions without T preference, suggesting alternative mechanism

## Results - Uncertain etiology signatures

- **Artifacts** (5 signatures): InD27, InD28, InD28m (related to SBS57), InD5, InD10
- **Correlated with other mutations** (3 signatures): InD31 (with SBS105), InD24 (with DBS8), InD12 (with DBS25)
- **C-insertion signatures**: InD26 and InD30 show different distributions at poly-C tracts
- **Tissue-specific**: InD4b, InD29 may be variants of known signatures requiring further investigation
- **Unique pattern**: InD23 shows tandem duplications ≥5bp, predominantly in bladder and colorectal cancers

## Results - PRRDetect classifier development

- **Rationale**: Current MSI detection methods insufficient, especially for polymerase mutants and non-epithelial tissues
- **Training set**: 571 GEL cancers (214 MMRd, 36 Pol-dys, 41 mixed, 280 PRR-proficient)
- **Model selection**: Tested multiple feature combinations; selected multinomial elastic net regression
- **Key features**: SBS and InD signature exposures for MMRd, Pol-dys, and mixed phenotypes; total InDel/SNV ratio
- **Validation performance**: AUROC=1.0, AUPRC=0.99 on independent validation cohort (1,351 samples)

## Results - PRRDetect performance comparison

- **Superior specificity**: Outperformed MSIseq, MMRDetect, and TMB-based detection
- **ICGC/Hartwig cohorts**: Identified 50/1,335 (3.7%) PRRd cases; correctly identified all Pol-dys and most MMRd
- **Missed by other methods**: MSIseq missed 6 MMRd and all 7 Pol-dys cases; MMRDetect missed 7 MMRd cases
- **Driver-independent detection**: 66% of PRRDetect-positive cases lacked identifiable driver mutations
- **TMB limitations**: In TMB-high samples (>10 mut/Mb), only 10.9% had actual PRRd; 89% had high TMB from other causes

## Results - PRRDetect application to GEL cohort

- **Large-scale analysis**: Applied to 4,775 GEL tumors; 1,371 were TMB-high
- **PRRd prevalence**: 49.4% of TMB-high cases had predicted MMRd/Pol-dys; ~50% of these lacked identified drivers
- **Tissue distribution**: Found PRRd in expected types (colorectal 19%, uterine 37%) and unexpected types (stomach 6%, bladder 1%, CNS 1%, lung 1%)
- **Clinical implications**: TMB-high status non-specific for immunotherapy response; PRRDetect provides more accurate biological stratification
- **Tumor-agnostic utility**: WGS with PRRDetect can identify PRRd across all cancer types

## Discussion

- **Classification impact**: Mutation classification directly affects signature analysis accuracy more than extraction algorithms
- **Key innovation**: Incorporating flanking sequence context and distributing signals across more channels improves discrimination
- **Biological validation**: New taxonomy captures true diversity of PRRd signatures matching corresponding substitution signature diversity
- **New discoveries**: 27 new InDel signatures identified; offered mechanistic insights for several
- **Future improvements**: Long-read sequencing could enable detection of InDels at longer simple repeats currently limited by short-read technology

## Discussion - Clinical applications

- **PRRDetect advantages**: Only tool distinguishing PRRd subtypes; highly sensitive and specific
- **Biomarker limitations**: TMB lacks specificity; >50% of TMB-high cancers lack immunotherapy-relevant biology
- **Hidden cases**: PRRDetect identifies ~50% of PRRd cases without detectable driver mutations
- **Genotype discrimination**: Currently doesn't distinguish MLH1 vs MSH2 vs MSH6 vs PMS2, but could be adapted if clinically needed
- **WGS value**: Comprehensive sequencing approach identifies actionable biology that targeted panels miss

## GPT-5

### Abstract & Introduction

- Small insertions/deletions (InDels) are the 2nd most common genomic variants after substitutions but less studied.

- CRISPR-edited human cellular models with mismatch repair (MMR) and polymerase ε/δ defects showed distinct InDel footprints.

- The standard COSMIC-83 InDel classification failed to distinguish these signatures.

- Authors developed a new classification system (89 subtypes) considering flanking sequences and motifs.

- Analyzed 100k Genomes Project data across 7 cancers, identifying 37 InDel signatures (27 novel).

- Built *PRRDetect* to classify tumors with postreplicative repair deficiency (PRRd), with clinical implications for immunotherapy.

---

### Background

- InDels (<100bp) are substantial, nonrandom, and reflect diverse mutational processes.

- Historically overshadowed by substitution signature studies.

- COSMIC-83 defined 18 InDel signatures but has limitations.

- Accurate InDel classification is clinically important for HR deficiency and microsatellite instability (MSI).

---

### Results – PRRd CRISPR Models

- Created isogenic CRISPR cell models: MMR knockouts (MLH1, MSH2, MSH3, SETD2), polymerase ε/δ proofreading mutants, and combinations.

- Most edits (except SETD2) greatly increased InDel burden (up to 300-fold for combined defects).

- COSMIC-83 profiles could not separate PRRd signatures from background.

- Different edits showed distinct InDel patterns (e.g., deletions at long poly-T tracts in MMRd; insertions in polymerase mutants).

---

### Limitations of COSMIC-83

- Aggregates long homopolymers (>5bp) into single bins, losing discriminatory signal.

- Leads to misattribution of MMRd patterns to “normal replication error” signatures ID1/ID2.

- Polymerase and MMRd edits appeared indistinguishable under COSMIC-83.

---

### New 89-Channel Classification

- Developed an expanded taxonomy using:
  
  - Insertion/deletion type and size.
  
  - Flanking bases.
  
  - Homopolymer length.
  
  - Microhomology at junctions.

- Started with 476 categories, consolidated to 89 informative channels.

- Better separation of CRISPR gene-edit signatures from controls.

- Revealed biologically meaningful patterns (e.g., MMR deletions enriched in long homopolymers; polymerase insertions enriched in shorter 5–7bp repeats).

---

### Improved Signature Extraction

- 89-channel scheme extracted more and cleaner de novo signatures compared with COSMIC-83.

- Showed distinct MMRd, polymerase, and combined mutational processes.

- Tested across algorithms and datasets (colorectal, ICGC) – always higher resolution.

---

### Application to Human Cancers

- Applied to 4,775 tumors from 7 cancer types (bladder, CNS, colorectal, endometrial, lung, stomach, skin).

- Identified 37 consensus InDel signatures (InDS):
  
  - 10 mapped to known COSMIC IDS.
  
  - 27 novel.

- Exogenous signatures: tobacco (InD3a/3b), UV (InD13), colibactin (InD18), platinum therapy (InD32).

- Endogenous signatures: replication slippage (InD1/2), HR deficiency (InD6), NHEJ (InD8), APOBEC (InD9a).

- Novel gene-specific signatures tied to MMRd and polymerase defects (InD14–21).

- 12 of uncertain/possible artifact origin.

---

### PRRDetect Classifier

- Built multinomial regression model using SBS + InDel signatures.

- Trained on 571 GEL tumors with known PRR status.

- Validated on 1,351 cancers: AUROC = 1, outperforming MSIseq, MMRDetect, and TMB.

- Detected PRR deficiency in tumors lacking identifiable driver mutations.

- Showed TMB is a poor surrogate biomarker: most TMB-high tumors not PRR deficient.

- PRRDetect identified PRRd across multiple tissues, not just colorectal/uterine cancers.

---

### Discussion

- Mutation classification choice strongly influences ability to discern true mutational processes.

- COSMIC-83 underestimates diversity; 89-channel taxonomy provides better separation and biological insight.

- Reveals diverse PRRd signatures, mirroring substitution SBS diversity.

- Future work: expand taxonomy further as sequencing improves.

- PRRDetect is highly specific and sensitive, and could refine immunotherapy patient stratification.

- Key clinical implication: TMB is overly broad; signature-based classifiers like PRRDetect give biologically grounded predictions.

## GPT-5 "table"

- Title: *A redefined InDel taxonomy provides insights into mutational signatures*

- Authors: Koh, Nanda, Rinaldi, Boushaki, Degasperi, Badja, Pregnall, Zhao, Chmelova, Black, Heskin, Dias, Young, Memari, Shooter, Czarnecki, Brown, Davies, Zou, Nik-Zainal.

- Abstract:
  
  - Small insertions and deletions (InDels) are common but understudied relative to substitutions.
  
  - Generated CRISPR-edited human cell models with postreplicative repair dysfunction (PRRd) including mismatch repair (MMR) and polymerase (Pol ε/δ) edits.
  
  - Found diverse InDel mutational footprints, but existing COSMIC-83 classification could not distinguish them.
  
  - Developed new classification using flanking sequence and motif context → 89 InDel subtypes.
  
  - Analyzed seven tumor types (100k Genomes Project) → uncovered 37 InDel signatures, 27 novel.
  
  - Built PRRDetect classifier to identify PRRd tumors, with potential immunotherapy relevance.

---

### Introduction

- InDels (<100 bp) are the second most common variation after substitutions.

- InDel mutagenesis is substantial, shaped by both physiological and pathological processes.

- Past decade focused mostly on substitutions; but advances in detection/annotation revealed 18 InDel signatures using COSMIC-83.

- COSMIC-83 system: 83 InDel subtypes defined by size, repeat context, and homology at junctions.

- A later reanalysis reported nine additional signatures.

- Accurate InDel classification is critical for biology and clinical prediction (e.g., HR deficiency, MSI detection).

---

### Results: InDel Diversity in PRRd Models

- Generated “ground truth” set by CRISPR editing in hTERT-RPE1 TP53-null cells.

- Edits included:
  
  - MMR knockouts: ΔMLH1, ΔMSH2, ΔMSH3, ΔSETD2.
  
  - Polymerase mutants: POLE P286R, POLE L424V, POLD1 S478N (exonuclease), POLD1 R689W (polymerase).
  
  - Double mutants combining polymerase proofreading with MMRd.

- Cultured for 45–50 days, then sequenced 2–5 daughter clones/genotype.

- InDel burden increased:
  
  - ~2× for ΔMSH3 and POLD1 R689W
  
  - 10× for POLD1 S478N, POLE mutants
  
  - 55× for ΔMSH2, ΔMLH1
  
  - 200–300× for combined edits.

- Different edits produced distinctive InDel profiles:
  
  - ΔMLH1/ΔMSH2/ΔMSH3 → 1 bp T deletions at long poly-T tracts.
  
  - POLE/POLD1 exonuclease mutants → 1 bp T insertions at poly-T tracts.
  
  - Double mutants mixed patterns.

- Same gene, different domain mutants gave distinct signatures (e.g., POLD1 exonuclease vs polymerase domain).

---

### Substitution Profiles

- MMRd → lower substitution/InDel ratio.

- Polymerase-dysfunction → higher ratio.

- Supports idea that genome instability differs by repair pathway: MMRd mainly InDels, Pol-dys mainly substitutions.

---

### Limitations of COSMIC-83

- COSMIC-83 lumps long homopolymer InDels into broad channels (e.g., “T6+ deletion”), losing discriminatory resolution.

- Leads to systematic misattribution:
  
  - ΔMSH2 and ΔMLH1 signatures misclassified as ID1/ID2 (replication slippage) instead of true MMRd.
  
  - Polymerase mutants indistinguishable from ID1.
  
  - POLD1 R689W had no match in COSMIC signatures.

- Net effect: COSMIC-83 cannot distinguish biological distinct InDel processes.

---

### Development of New Taxonomy

- Inspired by substitution classification, authors included:
  
  - Flanking sequence context (5′ and 3′ bases).
  
  - Motif size (1 bp vs ≥2 bp).
  
  - Homopolymer tract length.
  
  - Microhomology usage at junctions.
  
  - Repetitive unit features.

- Initial taxonomy: 476 categories (“channels”).

- Tested on ~18,500 tumors (ICGC/TCGA/Hartwig/GEL).

- Found many channels uninformative → consolidated to final **89-channel system (“InDel-89”)**.

- Compared to COSMIC-83:
  
  - Expanded informative 1 bp A/T InDels into more subcategories.
  
  - Condensed rare/uninformative categories.

- Result: higher resolution without overwhelming complexity.
