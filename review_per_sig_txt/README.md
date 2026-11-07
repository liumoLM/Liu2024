# Review of the per-signature notes (2026-11-06)

Scripts used to check the notes in `unified_indels_data/per_sig_txt/*.md`
(shown on the vignette pages) against the current signatures. The findings
are in `paper_to_do.md` under "Open: Steve".

Run `1_collapse_signatures.R` first. It uses the vignette's collapse
functions in `~/github/Liu2026_code_and_data` and writes `sigs.rds`
(untracked). The other scripts read `sigs.rds` and print their checks:

- `2_83type_features.R`: single-base T and C indel proportions, microhomology
  and repeat deletions, for the 83-type signatures and the 476-type
  signatures collapsed to 83 types.
- `3_89_476_features.R`: the same features in the 89-type signatures and in
  the 476-type signatures collapsed to 89 types, cosine similarities between
  signature pairs, and 476-type T and C indels by tract length.
- `4_cosmic_koh_assignments.R`: comparisons with COSMIC ID9 and ID4 and with
  Koh et al. signatures, ID_N as a mix of C_ID1 and C_ID2, the InsDel15 and
  InsDel16 exemplar spectra, and co-occurrence with MSI-H in the assignments.
- `5_raw_runs_id_n.R`: whether ID_N, InsDel_N_alpha, and InsDel_N_beta appear
  in the raw SigProfiler and mSigHdp runs.
