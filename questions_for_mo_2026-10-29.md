# Questions for Mo (and Mini), video meeting 2026-10-29

## 1. MSIseq cutoff

Mini's script `run_MSIseq_overlap_6975.R` calls `MSIseq.classify()` and
writes `MSI_classified_overlap_6975.tsv`. The script does not contain the
rule S.ind > 0.4 and T.ind > 1.2. That rule exists only as the column
`MSI_status_Sind_Tind_rule`, which was added by hand to
`mutationNum_overlap_6975.tsv`.

Steve ran `MSIseq.classify()` on the metrics in `mutationNum_overlap_6975.tsv`:

| | classify MSI-H | classify Non-MSI-H |
|---|---|---|
| rule MSI-H | 204 | 0 |
| rule Non-MSI-H | 70 | 6,701 |

So `MSIseq.classify()` gives 274 MSI-H. The rule keeps 204 of them and drops
70 borderline tumors (median T.ind 1.03, median S.ind 0.45).

`MSIseq.classify()` has no threshold argument. Its built-in classifier
(MSIseq 1.0.0, retrained from the package's training data to print it) is a
single split, S.ind > 0.394737, and on our tumors it matches that split
exactly. So the rule is MSIseq's own S.ind cutoff rounded to 0.4, plus an
added T.ind > 1.2 that MSIseq does not use. The added T.ind cutoff removes
most of the 70 tumors.

- Why add T.ind > 1.2, and how was 1.2 chosen?
- Why not use the `MSIseq.classify()` calls (274 MSI-H)?
- If we keep the rule, what do we say in the methods to justify it?

## 2. Topography rerun (Mini)

- Does the 1026 partial-credit file use the new assignments?
- When will Mini rerun r.01 to r.04 with it?
- Which folders (`results/` or `results0322/`) and columns (`ID.class` or
  `COSMIC_83`, `ID_O`/`ID_P`) made the current r.04 matrices?
- Is Mo redoing the simulations to drop the indels in genes on both strands?

## 3. Figure S1 with the new assignments

- 122 dots were added (923 to 1,029), 86 of them from 1 to 3 tumors (for
  example ID_L / InsDel_L in 8 cancer types). Are these expected from the
  mSigAct refits?
- Stomach InsDel1a rose more than tenfold (0.13 to 1.60 mean mutations/Mb)
  in the same 14 tumors. Why? (Answered 2026-10-29, no need to ask: one
  hypermutated tumor, CPCT02450014T, where InsDel1a replaced InsDel_D in
  the mSigAct refit. The new fit is better, and Steve kept it.)

## 4. Methods details

- mSigAct refit: mSigAct2 1.0.0 (commit b3e11fc) and mSigTools 1.0.8?
- Can Mo send the code that merged the mSigAct refits into the assignment
  tables?
- Hartwig data download date or release?

## 5. Age and sex analyses (Tables S15 and S16, Figure S6)

- Who writes the methods paragraph and the Table S16 legend?
- (Other decisions, if there is time: see "Open: Steve with Mo" in
  `paper_to_do.md`.)
