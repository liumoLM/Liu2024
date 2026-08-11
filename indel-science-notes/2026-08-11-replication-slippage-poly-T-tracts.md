# Replication slippage at homopolymer tracts: insertions vs deletions of T, and A/T vs G/C

Date: 2026-08-11

## Q1. Can slippage produce both insertions and deletions of T in poly-T tracts?

Yes, and this is one of the better-supported facts in the field. The two directions
are separate signatures with separate proposed mechanisms.

- **COSMIC ID1** = 1 bp **insertion** of T at poly-T. Proposed aetiology: "slippage
  during DNA replication of the replicated (nascent) DNA strand."
  https://cancer.sanger.ac.uk/signatures/id/id1/
- **COSMIC ID2** = 1 bp **deletion** of T at poly-T. Proposed aetiology: "slippage
  during DNA replication of the template DNA strand."
  https://cancer.sanger.ac.uk/signatures/id/id2/

  Both were extracted in Alexandrov LB et al. (2020) "The repertoire of mutational
  signatures in human cancer," *Nature* 578:94-101.
  https://www.nature.com/articles/s41586-020-1943-3 , doi:10.1038/s41586-020-1943-3

  The mechanistic logic: a loop-out on the **nascent** strand leaves an extra base in
  the product (insertion), a loop-out on the **template** strand skips a base
  (deletion). Same slippage event, opposite strand, opposite sign.

- Koh G et al. (2025) "A redefined InDel taxonomy provides insights into mutational
  signatures," *Nature Genetics* 57:1132-1141.
  https://www.nature.com/articles/s41588-025-02152-y , doi:10.1038/s41588-025-02152-y
  Experimentally separates the two directions: MMR-deficient lines are dominated by
  1 bp T **deletions** at poly-T6+, while polymerase-dysfunction lines give
  **exclusively** 1 bp T **insertions** at poly-T5+. Note the tract-length split:
  insertion excess peaks at 5-7 bp, deletion excess grows with longer tracts. They
  attribute this to proofreading acting only within ~5-7 bp of the active site.

- Germline confirmation in a real pedigree: Happ HC, Sasani TA, Warner D, Neklason DW,
  Quinlan AR (2026) "AVITI sequencing of a four-generation CEPH/Utah pedigree confirms
  low mutation rates at homopolymer loci despite their low sequence complexity,"
  *Genome Biology* 27:215.
  https://link.springer.com/article/10.1186/s13059-026-04099-7 , doi:10.1186/s13059-026-04099-7
  "Most mutations are 1 bp expansions or contractions, with a bias toward the former"
  (P = 1.4e-8), in both maternal and paternal de novo mutations. So both directions
  occur de novo in humans, with insertions modestly favoured.

- Mechanism at the polymerase level, two distinct routes to a 1 bp indel:
  Garcia-Diaz M, Kunkel TA (2006) "Mechanism of a genetic glissando: structural
  biology of indel mutations," *Trends Biochem Sci* 31:206-214.
  https://www.cell.com/trends/biochemical-sciences/fulltext/S0968-0004(06)00058-2 ,
  doi:10.1016/j.tibs.2006.02.004
  Classic Streisinger strand slippage plus **dNTP-stabilized misalignment**, where an
  incoming dNTP pairs correctly with a downstream template base across a looped-out
  template base. Evidence that the dNTP-stabilized route dominates at short tracts and
  gives way to slippage as tract length grows:
  https://www.sciencedirect.com/science/article/abs/pii/S156878641500052X

## Q2. Why does it hit poly-A/poly-T harder than poly-G/poly-C?

The honest answer: **mostly because of target abundance and tract length, not because
A:T base pairs are intrinsically more slippage-prone.** The intuitive 2-hydrogen-bond
explanation is contradicted by the per-locus data.

- Lujan SA, Clark AB, Kunkel TA (2015) "Differences in genome-wide repeat sequence
  instability conferred by proofreading and mismatch repair defects," *Nucleic Acids
  Research* 43:4067-4074.
  https://academic.oup.com/nar/article/43/8/4067/2414528 , doi:10.1093/nar/gkv271
  Far more indels are seen in A/T sequence genome-wide, but that reflects the greater
  abundance of A/T runs. Normalized for target size, "indel rates per base pair per
  generation are similar and sometimes slightly higher in runs of G/C base pairs than
  in runs of A/T base pairs of corresponding lengths." They note this is
  counterintuitive under a pure thermodynamic model and conclude slippage is limited
  "not only thermodynamically, but also kinetically and/or structurally."
  Also from this paper: indel rates rise **>100,000-fold** with homopolymer run length.
  Length is by far the dominant variable.

- Happ et al. 2026 (above) makes the same point in the human germline, more strongly:
  G/C homopolymers are under 1% of homopolymer loci but show an **18-fold higher
  per-locus mutation rate** than A/T homopolymers (rate ratio 0.056, 95% CI
  0.042-0.077, P < 2.2e-16). Caveat: their locus set is long homopolymers (>=10 bp),
  and a >=10 bp G/C run is an unusual sequence, so this is not a clean matched
  comparison across all lengths.

So the observed poly-A/poly-T dominance in human indel spectra is largely a **supply**
effect. The human genome is full of long A/T homopolymers, mainly the poly-A tails of
Alu and L1 insertions, while G/C runs are short and rare. Combine a
super-exponential length dependence with a length distribution shifted far to the
right for A/T, and A/T homopolymer indels dominate the genome-wide spectrum even at
equal or lower per-locus rates.

### Structural/dynamic arguments that are still made

- Bacolla A, Zhu X, Chen H, Howells K, Cooper DN, Vasquez KM (2015) "Local DNA
  dynamics shape mutational patterns of mononucleotide repeats in human genomes,"
  *Nucleic Acids Research* 43:5065-5080.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4446427 , doi:10.1093/nar/gkv364
  A-tracts and G-tracts have qualitatively different mutation profiles, not just
  different rates. A-tract mutations cluster at the first and last base pair of the
  tract and track local flexibility (peak fluctuation at the 5'-TA-3' junction).
  G-tract mutations cluster internally at G2-G3 and G5, and are explained by
  long-range hole transfer after ionization, that is oxidative chemistry, not slippage.
  Indels are more abundant at flanking base pairs than within tracts, ~10-fold more so
  for G-tracts. This argues the two tract types have different dominant mutagenic
  mechanisms, which is a cleaner framing than "A/T slips more."

- Structural: poly(dA).poly(dT) has a narrow, deep minor groove, high propeller twist
  and cross-strand hydrogen bonding, and is unusually straight and rigid. G-runs
  instead form G-quadruplex-type secondary structure, which is actively resolved by
  FANCJ/DOG-1 helicase, and loss of that helicase causes deletions initiating in G/C
  tracts. See the C. elegans work: https://pmc.ncbi.nlm.nih.gov/articles/PMC2211496/
  So G/C runs may be protected by a dedicated pathway rather than being intrinsically
  stable.

- Base stacking, not hydrogen bond count, dominates sequence-dependent duplex
  stability, which further undercuts the naive "2 vs 3 H-bonds" story.

## Which human polymerases are involved

Humans encode roughly 15 DNA-dependent DNA polymerases, catalogued in
Lange SS, Takata K, Wood RD (2011) "DNA polymerases and cancer," *Nature Reviews
Cancer* 11:96-110. https://www.nature.com/articles/nrc2998 , doi:10.1038/nrc2998
Only a handful carry out bulk synthesis.

### Bulk nuclear replication (B-family)

| Pol | Gene | Role | 3'->5' proofreading |
|---|---|---|---|
| Pol alpha | `POLA1` + primase | RNA primer then ~20 nt of DNA, at every origin and every Okazaki fragment | **No** |
| Pol delta | `POLD1` (+ POLD2/3/4) | Bulk lagging strand, maturation of Pol alpha tracts, most repair resynthesis | Yes |
| Pol epsilon | `POLE` (+ POLE2/3/4) | Bulk leading strand | Yes |

Division of labour reviewed in Burgers PMJ, Kunkel TA (2017) "Eukaryotic DNA
Replication Fork," *Annual Review of Biochemistry* 86:417-438.
https://www.annualreviews.org/doi/abs/10.1146/annurev-biochem-061516-044709 ,
doi:10.1146/annurev-biochem-061516-044709

**Pol gamma** (`POLG`) is the only polymerase that replicates mitochondrial DNA, so it
is the relevant enzyme for any mtDNA homopolymer indel analysis. PrimPol and Pol theta
also act on mtDNA to a minor extent.

### The rest, by job

- Gap filling and base excision repair: Pol beta (`POLB`), Pol lambda (`POLL`)
- NHEJ and V(D)J: Pol mu (`POLM`), Pol lambda, TdT (`DNTT`)
- Translesion synthesis: Pol eta (`POLH`), Pol iota (`POLI`), Pol kappa (`POLK`),
  Pol zeta (`REV3L`/`REV7`), Rev1 (`REV1`), Pol nu (`POLN`)
- Microhomology-mediated end joining: Pol theta (`POLQ`)
- Repriming past lesions: PrimPol (`PRIMPOL`)

### Which ones matter for homopolymer slippage

Pol delta and Pol epsilon, and to a lesser extent Pol alpha. Pol alpha is notably
error-prone and has no proofreading exonuclease, and it lays down a short DNA stretch
at the start of every Okazaki fragment, which is a plausible source of some
lagging-strand indels. The "Pol-dys" lines in Koh et al. 2025 that gave exclusively
1 bp T insertions at poly-T5+ are `POLE`/`POLD1` exonuclease-domain defects, and this
is the same system that motivates the proofreading-reach argument, that proofreading
only rescues errors within about 5-7 bp of the active site. The translesion
polymerases have far lower per-nucleotide fidelity but synthesize a tiny fraction of
the genome, so they contribute little to the genome-wide indel spectrum.

## Measured 1 bp insertion and deletion rates for Pol delta and Pol epsilon

### In vitro, purified enzyme

| Enzyme | Base substitution | 1 bp deletion | 1 bp insertion | Source |
|---|---|---|---|---|
| Yeast Pol delta, WT, 3-subunit | < 1.3e-5 | **3.0e-4** in homopolymeric runs | not reported as a notable class | Fortune 2005 |
| Yeast Pol epsilon, WT, 4-subunit | <= 2e-5 | **<= 5e-7** | not reported separately | Shcherbakova 2003 |
| Human Pol epsilon, exo-deficient | 4.4e-5 | **0.71e-5** | **0.16e-5** | Korona 2011 |

- Fortune JM, Pavlov YI, Welch CM, Johansson E, Burgers PMJ, Kunkel TA (2005)
  "Saccharomyces cerevisiae DNA polymerase delta: high fidelity for base substitutions
  but lower fidelity for single- and multi-base deletions," *J Biol Chem*
  280:29980-29987. doi:10.1074/jbc.M505236200
  https://www.jbc.org/article/S0021-9258(20)56429-6/fulltext
  Wild-type Pol delta "inefficiently proofreads single nucleotide deletion mismatches
  in homopolymeric runs, such that the error rate is 30 single nucleotide
  deletions/100,000 nucleotides polymerized." The paper's own conclusion: "strand
  slippage during replication by wild type Pol delta may be a primary source of
  insertion and deletion mutagenesis in eukaryotic genomes."

- Shcherbakova PV, Pavlov YI, Chilkova O, Rogozin IB, Johansson E, Kunkel TA (2003)
  "Unique error signature of the four-subunit yeast DNA polymerase epsilon,"
  *J Biol Chem* 278:43770-43780. doi:10.1074/jbc.M306893200
  https://pubmed.ncbi.nlm.nih.gov/12882968/
  Wild-type Pol epsilon proofreads at least 92% of base substitution errors and **at
  least 99% of frameshift errors**.

- Korona DA, LeCompte KG, Pursell ZF (2011) "The high fidelity and unique error
  signature of human DNA polymerase epsilon," *Nucleic Acids Research* 39:1763-1773.
  doi:10.1093/nar/gkq1034 https://academic.oup.com/nar/article/39/5/1763/2408989
  Exonuclease-deficient human Pol epsilon: deletion:insertion ratio about **4.4:1**.
  Also, human Pol epsilon was "nearly as accurate for 4-5 homonucleotide runs as it was
  for shorter runs or non-iterated nucleotides," unlike other B-family polymerases.

**The headline contrast:** wild-type Pol delta makes single-base deletions in
homopolymers at ~3e-4, roughly **600-fold** more often than wild-type Pol epsilon
(~5e-7). Both have comparable base-substitution fidelity. The difference is entirely
in frameshift proofreading, which Pol epsilon does well and Pol delta does badly. If
that in vitro difference holds in vivo, homopolymer indels should be a
predominantly **lagging-strand** phenomenon.

### In vivo, yeast whole genomes

Lujan SA, Clark AB, Kunkel TA (2015) *NAR* 43:4067-4074, doi:10.1093/nar/gkv271
(Table 1, MMR-deficient strains):

| Run type | 1 bp deletions | 1 bp insertions | del:ins |
|---|---|---|---|
| A/T | 4,686 | 424 | ~11:1 |
| G/C | 311 | 82 | ~3.8:1 |

"The number of deletions far exceeds the number of insertions, and the number of
deletions and insertions is much higher for A/T as compared to G/C base pairs." The
ratio is strongly length-dependent, reaching about **100-fold** deletion bias at 12 bp
A/T runs. Their proposed reason is thermodynamic, that forming the insertion-type
misaligned intermediate requires disrupting more hydrogen bonds than the deletion-type
one. Note these counts are from MMR-deficient backgrounds, so they measure the
polymerase error spectrum minus whatever MMR would have caught.

### Does MMR repair 1 bp deletions more efficiently than 1 bp insertions?

Short answer: no, and where an asymmetry exists it runs the other way for the complex
that handles most single-nucleotide loops.

Terminology, since it is easy to invert. In the MMR literature an **insertion loop** is
an extra nucleotide extrahelical on the **primer/nascent** strand, which if unrepaired
becomes an insertion mutation. A **deletion loop** is an extra nucleotide on the
**template** strand, which if unrepaired becomes a deletion.

- Romanova NV, Crouse GF (2013) "Different roles of eukaryotic MutS and MutL complexes
  in repair of small insertion and deletion loops in yeast," *PLoS Genetics* 9:e1003920.
  doi:10.1371/journal.pgen.1003920
  https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003920
  There **is** an asymmetry, and it is complex-dependent. Median repair ratios for 1-nt
  loops (their Table 4):

  | MutS complex | 1-nt insertion loops | 1-nt deletion loops | P |
  |---|---|---|---|
  | MutS-alpha (Msh2-Msh6, msh3 strains) | 190 | 73 | 0.017 |
  | MutS-beta (Msh2-Msh3, msh6 strains) | 3 | 22 | 0.029 |

  "MutSalpha has a consistently greater activity toward 1-nt insertion mismatches,
  whereas the MutSbeta activity is the reverse." They argue the pattern likely holds in
  mammalian cells. Since MutS-alpha is the main complex acting on single-nucleotide
  loops, the net effect of MMR should be to remove **insertion**-type errors somewhat
  preferentially, which would *deepen* rather than offset the polymerase deletion bias.

- Lujan SA, Clark AB, Kunkel TA (2015) *NAR* 43:4067-4074, doi:10.1093/nar/gkv271
  At the genome scale they see no evidence of differential MMR efficiency. The
  deletion-to-insertion ratio is similar in MMR-proficient and MMR-deficient strains,
  single-base deletion rates exceed insertion rates across genotypes, and "a high
  deletion to insertion ratio is also observed during synthesis *in vitro* by DNA
  polymerases." They attribute the bias to the polymerase step, specifically that
  "base-base hydrogen bonding must be disrupted for one more base pair in order to
  create an insertion mismatch...as compared to a deletion mismatch."

**Correction to an earlier note in this file.** I had suggested that MMR might
preferentially remove deletion-type intermediates, leaving an insertion-biased spectrum
in MMR-proficient cells. Both papers above contradict that. Lujan et al. find no
differential MMR efficiency at all, and Romanova and Crouse find MutS-alpha favouring
insertion loops, the opposite direction.

### The discordance that remains

So the insertion-biased human observations are still unexplained by differential MMR:

- Happ et al. 2026 germline de novo homopolymer mutations favour expansions
  (P = 1.4e-8).
- Koh et al. 2025 polymerase-dysfunction cell lines gave *exclusively* 1 bp T
  insertions at poly-T5+.
- COSMIC ID1 (insertion) is near-universal and clock-like, while ID2 (deletion) is the
  one that explodes in MMR deficiency.

Against in vitro and yeast in vivo data that are deletion-biased at every genotype
tested. Candidate explanations not yet checked here: the yeast reporter and genome
assays sample much shorter homopolymers than the human loci in question, short-read
alignment is known to call insertions and deletions at homopolymers with different
sensitivity, and the human germline number comes from a specific long-allele locus set
(>=10 bp). A direct test in our data is the del:ins ratio stratified by homopolymer
length in MMR-proficient versus MMR-deficient tumours, which shows whether the human
insertion bias is real, length-specific, or a calling artefact. That test has now been
run, see the next section.

## Our own data: the discordance is length-specific

The test proposed above has now been run on the Liu et al. ID476 cohort. Script,
augmented table and figures are in
`~/github/zz-liu_2025_draft_release/del-vs-ins-analysis/` (commit 11ff2ae).

Method. Transpose the 476-channel spectra (6,978 tumours), aggregate the 342 one-base
homopolymer channels across their 9 flanking contexts into per-tract-length deletion
and insertion counts, and split by `MSIseq_MSI.H`, which gives 202 MMR-deficient
(MSI-H), 6,773 MMR-proficient and 3 unknown. Panel Rn compares Del(base):Rn with
Ins(base):Rn, both meaning a tract of n bases existed **before** the event, so the two
channels describe the same substrate, which is confirmed by the ICAMS classifier: `R`
is the tract length before the mutation for insertions as well as deletions (see the
caveats below). No indel-burden filter is applied, because a 4,000-indel cutoff would
remove 200 of the 202 MSI-H tumours.

### Pooled insertions / deletions by tract length

Values above 1 mean insertions exceed deletions. Pooled across tumours rather than a
median of per-tumour ratios, since most tumours are zero in most channels.

| Tract | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | R9+ |
|---|---|---|---|---|---|---|---|---|---|
| poly-T, MMR proficient | 0.59 | 0.46 | 0.45 | 0.31 | 0.83 | **2.29** | **1.79** | **1.64** | 0.94 |
| poly-T, MSI-H | 0.46 | 1.48 | 0.61 | 0.39 | 0.35 | 0.60 | 0.33 | 0.26 | 0.15 |
| poly-C, MMR proficient | 0.10 | 0.07 | 0.10 | 0.22 | 0.28 | 0.56 | **2.17** | **1.77** | 0.45 |
| poly-C, MSI-H | 0.10 | 0.13 | 0.31 | 1.01 | 0.25 | 0.25 | 0.51 | 0.47 | 0.19 |

### Coarse bins at the Koh et al. 2025 thresholds

Koh et al. report 1 bp T deletions at poly-T6+ in MMR-deficient lines and 1 bp T
insertions at poly-T5+ in polymerase-dysfunction lines, so these two bins are at
different tract lengths and this ratio is **not** a same-substrate comparison.

| | deletions from 6+ | insertions into 5+ | ins/del |
|---|---|---|---|
| poly-T, MMR proficient | 811,202 | 1,310,978 | **1.62** |
| poly-T, MSI-H | 7,245,469 | 2,015,686 | **0.28** |
| poly-C, MMR proficient | 49,502 | 66,125 | **1.34** |
| poly-C, MSI-H | 401,904 | 183,258 | **0.46** |

### What this shows

1. **The two MMR states diverge with tract length, in opposite directions.** In
   MMR-proficient tumours the ratio climbs with length and crosses 1 at poly-T6,
   peaking near 2.3. In MSI-H tumours it falls monotonically past R2, reaching 0.15 at
   R9+, so deletions outnumber insertions about 7:1 at long tracts. On the coarse bins
   that is a 5.8-fold swing in ins:del between MMR states for poly-T, and 2.9-fold for
   poly-C.

2. **This reproduces Koh et al. 2025 in an independent cohort**, from raw channel
   counts with no signature extraction. Their polymerase-dysfunction insertions at
   poly-T5+ and MMR-deficient deletions at poly-T6+ appear here as an MMR-proficient
   insertion excess starting at exactly poly-T6 and an MSI-H deletion excess growing
   over the same range. It is ID1 versus ID2 falling out of the counts directly.

3. **The human insertion bias is length-specific, not global.** At short tracts (R1 to
   R4) deletions dominate in both MMR groups, which agrees with the yeast and in vitro
   measurements. The insertion excess appears only at R6 to R8. That partly reconciles
   the discordance above: the yeast reporter and genome assays sample short
   homopolymers, exactly the range where our human data agree with them, so the two
   bodies of work may not actually conflict.

### Where InsDel_N_beta sits on this axis

The ID89 signature `InsDel_N_beta` turns out to be the interesting one for the
slippage question. Full analysis with code in
`~/github/zz-liu_2025_draft_release/del-vs-ins-analysis/n_beta_ins_del.qmd`.

**N_beta is the balanced long-poly-T signature.** Taking the insertion and deletion
content of each ID89 signature profile:

| Signature | % insertions | % deletions |
|---|---|---|
| InsDel1a | 93.7 | 6.3 |
| InsDel2a | 25.2 | 74.8 |
| **InsDel_N_beta** | **42.8** | **57.2** |

1a is essentially pure insertion and 2a mostly deletion, the two directional extremes.
N_beta sits between them, and its top channels interleave both directions at the same
long tracts: `A[Del(T):R(8,)]A`, `A[Del(T):R(5,7)]A`, `A[Ins(T):R(5,7)]C`,
`A[Ins(T):R(8,)]A`. That is what an undirected slippage process should look like, a
tract of 8 T's going to 7 about as often as it goes to 9.

**Group ratios follow from mixing.** Restricting to MMR-proficient tumours, so that
burden and MSI status do not drive the contrast, and using the coarse bins
(ins into 5+ over del from 6+):

| Group | n | median indel burden | poly-T ins/del |
|---|---|---|---|
| no N_beta | 4,143 | 494 | 3.05 |
| N_beta + 1 family only | 1,318 | 966 | 1.82 |
| N_beta alone | 1,191 | 1,096 | 1.00 |
| N_beta + 2 family only | 121 | 1,848 | 0.56 |

The third row is the consistency check. Tumours carrying N_beta and neither directional
family land at 1.00, exactly what a 43/57 signature predicts. Adding the insertion
family pushes to 1.82, adding the deletion family pulls to 0.56.

**N_beta is MMR proficient.** Only 10 of 2,640 N_beta tumours are MSI-H (0.4%),
against 4.4% of tumours without it.

**Co-occurrence is lopsided but that is not the interesting part.** Of the 2,640 N_beta
tumours, 1,195 carry neither directional family, 1,319 carry the 1 family only, 124 the
2 family only, and just 2 carry both. The 1 and 2 families are near mutually exclusive.
The 10:1 skew toward the 1 family mostly reflects how common that family is: among the
4,143 MSS tumours **without** N_beta, 3,096 carry the 1 family only, 91 the 2 family
only, 62 both, and 894 neither.

**Why this matters here.** The report above frames the problem as a conflict between
deletion-biased enzymology and insertion-biased human data. N_beta suggests the human
data contain a component that is neither. If N_beta is the underlying slippage process
and 1a and 2a are what happens when something biases it, then the quantity to explain
is not "why are human indels insertion-biased" but "what tips a balanced process one
way or the other". MMR status is one such tipping factor and is visible in the tables
above, but N_beta tumours are almost all MMR proficient, so it cannot be the only one.

### Caveats

- **The R9+ drop, checked.** The ratio falls back below 1 at R9+ in both MMR-proficient
  rows. Two candidate explanations were tested and one survives.

  *Not a binning artefact.* The concern was that `Del:R(9,)` and `Ins:R(9,)` might pool
  different tract-length distributions. They do not. In ICAMS
  `categorize_1_justified_indel()` the repeat count `R` is the tract length **before**
  the mutation for insertions as well as deletions, the insertion path explicitly
  subtracting one because the context it matches against is post-insertion:

  ```r
  indel_str_count_in_ref = nchar(all_repeated_seq) / ins_or_del_seq_len
  if (ins_or_del == "i") {
    indel_str_count_in_ref = indel_str_count_in_ref - 1
    # x_context was the sequence after the insertion. We want
    # indel_str_count_in_ref to reflect the repeat count prior to the insertion
  }
  ```

  So both channels draw on the same set of genomic tracts of length >= 9 and the bins
  are symmetric.

  *Not pipeline-specific.* Splitting the MMR-proficient tumours by cohort gives nearly
  identical curves for HMF (n = 4,038) and PCAWG (n = 2,735), poly-T ins/del:

  | Cohort | R4 | R5 | R6 | R7 | R8 | R9+ |
  |---|---|---|---|---|---|---|
  | HMF | 0.316 | 0.835 | 2.262 | 1.753 | 1.607 | 0.939 |
  | PCAWG | 0.312 | 0.819 | 2.362 | 1.867 | 1.762 | 0.964 |

  Two cohorts with different sequencing and calling pipelines agree to within about 3%
  at every bin, including the drop. The MSI-H split behaves the same way.

  *What is left.* Either the ins/del ratio genuinely declines above tract length 9,
  making the MMR-proficient curve unimodal with a peak at poly-T6, or short-read
  sequencing loses 1 bp insertions faster than 1 bp deletions once the tract exceeds
  what a read can span. The cohort check cannot separate these, since HMF and PCAWG are
  both Illumina short-read and would share any technology-wide limitation. That is the
  same `k` problem as below, concentrated at the longest tracts, and it needs long-read
  or duplex data to settle.
- poly-C at R7 and R8 rests on small counts.
- **Calling sensitivity, and which of our claims it threatens.** Short-read calling
  need not detect 1 bp insertions and 1 bp deletions at long homopolymers equally well.
  Write `k` for the ratio of insertion sensitivity to deletion sensitivity. Every
  ins/del value in the tables above is then the true ratio multiplied by `k`.

  The **between-group comparison survives any value of `k`.** Both groups pass through
  the same pipeline, so `k` cancels when the two ratios are divided:
  1.62 / 0.28 = 5.8 for poly-T regardless of `k`. The 5.8-fold swing between MMR states
  is the robust result.

  The **absolute position of the curve does not survive.** "Insertions exceed deletions
  at poly-T6 in MMR-proficient tumours" compares an observed ratio against the fixed
  threshold of 1, and `k` does not cancel against a constant. At k = 0.5 the true
  poly-T6 ratio would be 4.6 rather than 2.29, at k = 2 it would be 1.15, and the
  tract length at which the curve crosses 1 slides accordingly. So the crossing point
  landing at exactly poly-T6, which is what makes the agreement with Koh et al. 2025
  look tidy, is the weaker of the two claims.

  Two limits on this argument itself. We have not established which direction `k` runs,
  only that there is no reason for it to be 1. And "internally controlled" assumes
  calling behaves the same way in both groups, which is imperfect here, since MSI-H
  tumours carry roughly an order of magnitude more indels and may pass through
  burden-dependent filters differently. Settling this would need a truth set, for
  example long-read or duplex-sequenced homopolymer calls in a few samples from each
  group.

## Bottom line

1. Both +T and -T at poly-T are real and mechanistically distinguished: ID1 (nascent
   strand slippage, insertion) vs ID2 (template strand slippage, deletion), with
   experimental separation in Koh et al. 2025 and de novo germline confirmation in
   Happ et al. 2026.
2. There is **no good evidence that A:T base pairs slip more readily than G:C at
   matched tract length.** The genome-wide poly-A/poly-T dominance is driven by tract
   length and tract abundance. Anyone claiming an intrinsic A/T slippage preference
   should be pointed at Lujan et al. 2015 and Happ et al. 2026.
3. **MMR does not preferentially repair one direction.** Lujan et al. 2015 find the
   del:ins ratio is unchanged between MMR-proficient and MMR-deficient strains, and
   Romanova & Crouse 2013 find MutS-alpha favours insertion loops if anything.
4. **In our cohort the ins:del ratio depends on tract length and MMR state together.**
   Deletions dominate at short tracts in both MMR states, matching the yeast and in
   vitro data. Insertions overtake deletions only at poly-T6 to poly-T8 and only in
   MMR-proficient tumours. So the apparent conflict between the human insertion bias
   and the deletion-biased enzymology is largely a matter of which tract lengths each
   body of work samples.

## Related material in this folder

- The indel signature classification paper (`README.md`, `long-abst.md`, `abstract-150.md`).
  Relevant to that work: the ID1/ID2 nascent-vs-template split, the tract-length
  dependence of insertion vs deletion excess (Koh et al. 2025), and the point that
  A/T dominance is a target-abundance effect, which bears on opportunity
  normalization for any new indel classification.
- [`comparison-with-koh-2025.md`](comparison-with-koh-2025.md)
- The cohort analysis behind the "Our own data" section lives outside this repo, in
  `~/github/zz-liu_2025_draft_release/del-vs-ins-analysis/` (script, augmented
  6,978 x 524 table, and a 4-page PDF of the scatters and ratio summary).
