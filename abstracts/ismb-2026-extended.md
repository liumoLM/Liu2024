## Interpreting indel mutational signatures in 6,975 tumors based on 
more informative indel classification schemes

Mo Liu, Mi Ni Huang, Xue Ming Wu, Qi Zheng, Runtian Yao, Ying Yang, 
Runxi Shen, Steven G. Rozen

### Background

Mutational signatures are distinctive patterns of mutations imprinted 
on genomes by mutagenic processes or exposures. They can be identified 
experimentally—by exposing cells, organoids, or animal models to 
suspected mutagens—or computationally, as latent factors extracted 
from large collections of somatic mutation data. Mutational signatures 
are important for variant interpretation in cancer. They can reveal 
exposures to exogenous carcinogens (e.g., tobacco smoke, aristolochic 
acid) or endogenous mutagenic processes (e.g., defective mismatch 
repair, defective homologous recombination, activation of
APOBEC cytidiine deaminases). They can guide treatment decisions (e.g., identifying 
homologous recombination deficiency). They can help distinguish 
elevated mutation rates from selection by estimating the rates of 
different mutation types.

While single-base substitution (SBS) signatures have been extensively 
studied (as evidenced by 78 SBS signatures in COSMIC v3.5, https://cancer.sanger.ac.uk/signatures/sbs/),
the mutational signatures of small 
insertions and deletions (indels) have been comparatively 
under-studied, with only 25 indel signatures in COSMIC. This gap is 
significant because indel signatures carry important information for 
cancer classification, epidemiology, and understanding DNA repair 
mechanisms. Until recently, only one indel classification scheme was 
widely used. This classification scheme recognizing 83 types of 
indel based on the number of bases 
affected, the identity of the inserted or deleted base, flanking repeat 
context, and the presence of microhomology. Koh, Nanda, and colleagues 
recently proposed two new classification schemes recognizing 89 and 476 
types of indels, respectively (https://doi.org/10.1038/s41588-025-02152-y). 
Importantly, the new classification schemes distinguish 
different single-base T or C indels according to their non-T (or non-C) 
flanking bases, which substantially improves discrimination between 
mutational signatures.

### Results

We analyzed somatic indels in 6,975 whole-genome-sequenced tumors 
spanning 32 cancer types. We performed de novo signature discovery in all three 
classification schemes (83-type, 89-type, and 476-type) using two 
independent algorithms: mSigHdp, based on hierarchical Dirichlet 
processes, and SigProfilerExtractor, based on non-negative matrix 
factorization. We systematically mapped signatures across the three 
classification schemes by collapsing 476-type signatures to 89-type 
counterparts.
Because lossless translation between the new classification schemes
and the more common, 83-type classification is not possible,
we used tumor spectra dominated by individual 
signatures to link 83-type signatures to signatures in the two new 
schemes. We further profiled the genomic topography of each signature, 
including transcription strand asymmetry, replication strand asymmetry, 
genic versus intergenic enrichment, and replication timing effects.

We identified 44 mutational signatures in the 89-type classification 
and, independently, their 44 corresponding signatures in the 476-type 
classification, with high concordance between the collapsed 476-type 
and directly extracted 89-type signatures (mean cosine similarity > 
0.989). In the 83-type classification, we identified 34 signatures, of 
which 20 matched COSMIC reference signatures. Of the 44 signatures in 
the 89-type system, 18 corresponded to previously described signatures.

The new classification schemes proved superior in most cases. Several 
signatures that were merged in the 83-type system could be resolved in 
the 89/476-type systems. For instance, the 83-type ID1 corresponded to 
four distinct 89-type signatures (InsDel1a–d), each displaying 
different cancer-type-specific activity patterns and distinct 
correlation profiles with SBS signatures. Conversely, in two cases the 
83/476-type systems distinguished signature pairs that were conflated 
in the 89-type system, illustrating that the three schemes provide 
complementary information.

A key finding was the identification of a new signature, InsDel_F, 
which more accurately captures the mutational footprint of RNase H2 
deficiency than the previously proposed ID4. InsDel_F has a distinctive 
profile of 2-base deletions with microhomology (predominantly TCT→T 
or TGT→T) and 2-base deletions from paired doublets, while 
conspicuously lacking the 3- and 4-base deletions present in ID4. This 
profile closely matches experimental data from both RNase H2-deficient 
cell lines and mouse tumors.

Correlation analysis between indel and SBS signature activities 
revealed biologically coherent modules, including clusters associated 
with APOBEC activity, homologous recombination deficiency, UV exposure, 
mismatch repair deficiency, tobacco smoking, and gastrointestinal 
reactive oxygen species. Topographic analysis showed that 17 of 33 
signatures exhibited transcription strand asymmetries and 12 showed 
replication strand asymmetries, consistent with the known biology of 
their associated mutagenic processes.

### Implications for variant interpretation

These indel signatures have direct implications for variant 
interpretation in cancer. First, identifying the mutagenic processes 
active in a tumor helps distinguish driver mutations from passengers: 
mutations consistent with an active signature's profile are more likely 
to reflect elevated background mutation rates rather than selection, 
while mutations inconsistent with any active signature warrant closer 
scrutiny as potential drivers. Second, signatures serve as biomarkers 
for clinically actionable mutagenic states. For example, the mismatch 
repair deficiency signatures can flag tumors likely to respond to immune checkpoint 
inhibitors. Similarly, the homologous recombination deficiency signatures 
(C_ID6, InsDel6) identify candidates for PARP inhibitor therapy. Third, 
indel signatures can reveal carcinogen exposures—such as the 
aristolochic acid signatures that show distinctive flanking-base 
preferences in the new classification schemes—informing cancer 
prevention efforts. The finer resolution of the 89-type and 476-type 
schemes enhances these applications by resolving previously merged 
signatures into biologically distinct components with different 
cancer-type distributions, different SBS correlations, and therefore 
different implications for the interpretation of individual variants.

### Availability

We provide a website for exploring indel signatures and their 
relationships across the three classification schemes, including 
signatures' relationships to genomic topography, correlations with SBS 
signatures, and cancer-type-specific activity patterns. 

At the time
of submission of this absract we present offer overview of the indel mutational
signatures repertoire at https://doi.org/10.5281/zenodo.18451842. This
page will point to a more dynamic web site when it is ready as well as
repositories of the analysis code.



