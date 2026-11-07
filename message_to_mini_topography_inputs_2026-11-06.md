Subject: ID83 topography input files

Hi Mini,

Thank you again for the ID83 scripts and `append_bothstrand.R`. They are now in the code repo.

For our records, could you send us copies of the topography inputs that we do not have here?

1. `all_simulated_with_annotation.txt`, and `all_simulated_with_annotation.with_bothstrand.txt` with the summary and log files that `append_bothstrand.R` wrote next to it
2. `single.indel.aggregate.strand.added.matrix.txt` (the strand reference that `r.03_mergeInf.R` reads)
3. `all.meta.0528.rds`

These have per-indel data of controlled-access tumors, so please send them by a private link (for example MEGA or Duke Box), not by email. Compressed files are fine.

The reason is that the inputs name the Hartwig tumors by Hartwig internal IDs (CPCT..., DRUP..., WIDE...), which cannot appear in a publication. We plan to convert every input to HMF IDs here, take cancer types from Table S16 instead of `all.meta.0528.rds`, and rerun `normalize.sh` to r.04 to check that the r.04 matrices match your rerun of 2026-10-30. So you do not need to change anything on your side.

Best,
Steve
