Subject: Methods and code for the topography inputs and the TCGA exomes

Hi Mo,

We are documenting the topography analysis so that a reader with data access could rerun it. A few inputs and steps are not yet in the code repo or the methods. Could you help with these?

Per-indel signature probabilities (`all.indel.ID83.partial.credit.1026.txt.gz`):

1. The script that made the 1026 files, and its inputs (which annotated VCFs, which exposures, which signature profiles).
2. One or two sentences for the methods on how a probability is computed. For example, for an indel of class c in tumor t, is the probability for signature s equal to exposure(t, s) × profile(s, c) divided by the sum over all signatures?
3. Are the 83-type probabilities the 89-type ones summed over the 89-to-83 mapping, as for the exposures?

Simulated indels (`all_simulated_with_annotation.txt`): you sent the SigProfilerSimulator command on 2026-10-07. We still need:

4. The input folder you passed to it (which VCFs).
5. The seed. As far as we know, SigProfilerSimulator uses its default seed file unless one is given.
6. The script that turned its output into `all_simulated_with_annotation.txt` and annotated it with ICAMS.

Strand reference:

7. The script that made `single.indel.aggregate.strand.added.matrix.txt` (the +/-/Q strand of each real indel), or who made it.

TCGA exomes:

8. The spectra include 3 TCGA whole exomes (`WES::TCGA-EW-A2FV-01A-11D-A17D-09`, `WES::TCGA-D8-A27V-01A-12D-A17D-09`, and `WES::TCGA-13-0889-01A-01W-0420-08`). They have no signature assignments, and the manuscript does not mention them. Where did their indel calls come from, and were they used in signature extraction? The data repo README says that `TCGA-EW-A2FV` is the only tumor with its signature. Which signature is that? We need a sentence on the exomes in the data sources and the extraction methods.

Best,
Steve
