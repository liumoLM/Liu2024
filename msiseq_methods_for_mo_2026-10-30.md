# MSIseq methods: please check (2026-10-30)

Mo, this is the proposed methods text for the MSI-H calls, now in the
section "MSI/MSS status and high/low TMB status" of `ms.qmd`. Is it
correct?

> MSIseq (v1.0.0) metrics were computed with `Compute.input.variables()`
> on the SBSs and indels of all 6,975 tumors (`uniform.seq.len` = 3,200
> Mb), where S.ind is the number of indels in simple sequence repeats per
> Mb and T.ind is the total number of indels per Mb. We first called
> tumors MSI-H when S.ind > 0.4 (more than about 1,280 indels in simple
> repeats). Inspection of the spectra of tumors near this threshold
> revealed likely false positives: their indels in simple repeats were
> mostly single-T insertions and deletions in long poly-T tracts
> (InsDel_N_beta and InsDel1a), their total indel burden was low, and they
> showed no SBS signatures typical of mismatch repair deficiency. We
> therefore additionally required T.ind > 1.2 (more than about 3,840
> indels), which gave 204 MSI-H tumors.

## Questions

1. Which SBS signatures did you look for in the suspected false positives
   (for example SBS6, SBS15, SBS21, SBS26, SBS44)? Should the text name
   them?
2. How was 1.2 chosen? For example, was it set to exclude the false
   positives you found in the spot checks?
3. The sentence about poly-T tracts comes from our own check, not from your
   spot checks (see below). Does it match what you saw?

## Our check of the tumors removed by T.ind > 1.2

S.ind > 0.4 alone gives 271 MSI-H tumors. T.ind > 1.2 removes 67 of them,
leaving 204. The 67 removed tumors:

- have T.ind from 0.69 to 1.20 (about 2,200 to 3,800 indels). The 204 kept
  MSI-H tumors have a median T.ind of 12.3.
- have InsDel_N_beta as their largest signature in 53 of 67 (InsDel2a in
  10, InsDel14 in 4). InsDel_N_beta plus InsDel1a make up a median of 76%
  of their indels.
- are mostly colon (29) and esophagus (27) tumors.

The rule is applied by `msiseq/call_msi_h.R` in the code repo, from your
code of 2026-10-30.
