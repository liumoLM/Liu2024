# Age and sex models for Tables S15 and S16

Oct 22, 2026 · @Steve

Within cancer types, C\_ID1/InsDel1a, C\_ID5/InsDel5a and C\_ID9/InsDel9 are clearly age-associated, and the effect survives adjustment for sex, cohort and MSI-H status. InsDel5b, which looked age-associated when all tumors were pooled, is not. Age is a real but small contributor: cancer type explains far more of the variation between tumors.

## The problem with the pooled correlation

The "Overall" row of Table S16 was one Spearman correlation of signature activity with age over all tumors pooled. Cancer types differ both in patient age and in signature activity, so a pooled correlation can appear with no age effect inside any cancer type.

InsDel5b is the example. The manuscript said it was age-associated overall but in no single cancer type. Combining the per-type correlations from Table S16 (Fisher z, weighted by n − 3) removes it:

| Signature | Pooled rho (FDR) | Within-type rho (p) | Types with rho > 0 |  |
| --- | --- | --- | --- | --- |
| InsDel5b | 0.47 (5e-4) | 0.19 (0.26) | 3 of 5 |  |
| InsDel5c | 0.33 (6e-4) | 0.26 (0.01) | 5 of 9 |  |
| C\_ID4/InsDel4 | 0.23 (5e-7) | 0.25 (3e-8) | 12 of 18 |  |
| C\_ID5 | 0.22 (3e-11) | 0.34 (5e-29) | 17 of 20 |  |
| InsDel5a | 0.16 (8e-6) | 0.34 (7e-25) | 15 of 20 |  |
| C\_ID1 | 0.21 (4e-28) | 0.18 (7e-22) | 15 of 20 |  |
| InsDel1a | 0.20 (2e-20) | 0.14 (5e-12) | 15 of 20 |  |

Only 59 tumors have InsDel5b, mostly liver, prostate and pancreas, and the per-type correlations point in different directions. The sex analysis (Table S15) is per cancer type and has no pooled row, but it is not adjusted for age.

## Already changed (committed 2026-10-22)

These use the new signature assignments adopted the same day.

- `build_fig_s6/build_fig_s6_stats.R` no longer computes the pooled Overall correlation. Its new `age_glm_table()` fits exposure \~ age + cancer type and presence \~ age + cancer type (no sex, cohort or MSI-H yet) and writes `age.glm.adjusted.for.cancer.type.txt`.
- Table S16 now has two sheets: "By cancer type" (the per-type Spearman correlations, as before) and "Adjusted for cancer type" (the new models).
- Figure S6B no longer has the Overall row.
- `ms.qmd`: the InsDel5b sentences were deleted. A new sentence (green, for review) says that with cancer type as a covariate the exposures of C\_ID1 and InsDel1a rise 8% per decade of age and those of C\_ID5 and InsDel5a 19%, and that older patients are more likely to have them. The numbers come from `ms_numbers.R`.

The joint models with sex, cohort and MSI-H below are a prototype, not yet used by any figure or table. The script is build\_fig\_s6/age\_sex\_glm\_prototype.R in the code repo, and its results are in build\_fig\_s6/age\_sex\_glm\_prototype.tsv.

## The joint models

Two models per signature, fit to the 4,960 HMF and PCAWG tumors with known age and sex (CNS excluded, as in the existing age analysis). Covariates: age in decades, sex, cohort (PCAWG vs HMF), MSI-H by MSIseq, and cancer type as a fixed effect.

1. **Exposure** (the main test of clock-like behavior): negative binomial regression of the attributed indel count on the covariates, in tumors where the signature is present, in cancer types with at least 3 such tumors.
2. **Presence**: logistic regression of presence on the covariates plus log10 total indels, in all tumors of cancer types where the signature is present in some but not all tumors. Total indels is included because the chance of detecting a signature depends on how many indels a tumor has.

Details:

- Present = more than 1% of the tumor's attributed indels, as in the sex analysis.
- Each covariate is tested by dropping it from the full model (likelihood-ratio test). P values are BH-adjusted across signatures for each model and covariate.
- A model is fit only if at least 30 tumors have the signature.
- Total indels is not in the exposure model: burden rises with age partly because of these signatures, so adjusting for it would remove the effect being tested.
- Each 83-type/89-type pair with identical exposures (for example C\_ID9/InsDel9) gives identical results.

On SHAP: for regression models SHAP values restate the coefficients and give no test, and for a machine-learning model they would mostly report that cancer type explains most of the exposure. Sex is nearly determined by cancer type for breast, prostate, ovary and uterus, so SHAP's split of credit between them would be arbitrary. The drop-one deviance below answers the variable-importance question for these models.

## Age effects after adjustment

The four signatures the manuscript calls clock-like, plus C\_ID9/InsDel9, are age-associated in both models. Exposure effects are fold change per decade; presence effects are odds ratios per decade; FDR in parentheses.

| Signature | Exposure fold per decade (FDR) | Presence OR per decade (FDR) | Reading |
| --- | --- | --- | --- |
| C\_ID5 | 1.19 (2e-23) | 1.29 (7e-12) | Clock-like |
| InsDel5a | 1.19 (4e-19) | 1.41 (3e-16) | Clock-like |
| C\_ID1 | 1.10 (2e-15) | 1.17 (2e-7) | Clock-like |
| InsDel1a | 1.10 (4e-12) | 1.19 (7e-9) | Clock-like |
| C\_ID9/InsDel9 | 1.07 (2e-3) | 1.18 (8e-6) | Clock-like |
| C\_ID4/InsDel4 | 1.09 (3e-4) | 1.11 (0.04) | Weaker, both models |
| C\_ID10/InsDel10 | 1.26 (2e-5) | 0.95 (0.7) | Exposure only, 128 tumors |
| ID\_A | 1.18 (2e-5) | 0.95 (0.5) | Exposure only |
| InsDel5c | 1.21 (2e-4) | 1.00 (1.0) | Exposure only |
| InsDel\_A\_alpha | 1.21 (3e-4) | 0.90 (0.2) | Exposure only |
| InsDel1c | 1.09 (5e-3) | 0.99 (0.9) | Exposure only |
| ID\_N | 1.08 (2e-5) | 0.90 (6e-3) | Opposite directions |
| InsDel\_N\_beta | 1.07 (5e-5) | 0.89 (4e-3) | Opposite directions |
| InsDel5b | 1.19 (0.05) | 0.84 (0.2) | Not age-associated |

The ID\_N exposure effect survives adjustment for MSI-H status, but presence decreases with age. In the earlier age + cancer type model, ID\_N also had a strong age × cancer type interaction (p = 3e-19), so its age effect differs between cancer types. Following it up would mean per-type age estimates from the interaction model, before any claim in the text.

Across the middle half of patients (ages 51 to 70), a 19% rise per decade means about 1.4-fold more C\_ID5 exposure, and from age 30 to 80 about 2.4-fold.

## How much age explains

Age is highly significant but explains about 1–4% of the variation in exposure, while cancer type explains 15–60%. The tiny p values come from the number of tumors (hundreds to thousands per model), not from a large share of the variation.

Two measures:

- **Likelihood-ratio share**: the chi-square for dropping age (1 df) divided by the chi-square of the full model against a model with no covariates. Shares of different covariates do not add to 1, because covariates overlap and each share is that covariate's unique contribution.
- **Deviance R²**: 1 − deviance(full) / deviance(no covariates). For the negative binomial, all models are refit with the full model's dispersion so deviances are comparable. A covariate's unique R² is how much R² drops when it alone is removed. For presence this equals McFadden's R². For exposure McFadden's R² is 0.03–0.07 because count-level noise dominates the negative binomial log-likelihood, so it is not informative there.

| Model | Signature | Tumors | Age χ² (p) | Age share of model χ² | R² full model | Unique R²: age | Cancer type | MSI-H | Cohort |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Exposure | C\_ID5 | 1,073 | 107 (4e-25) | 0.11 | 0.58 | 0.042 | 0.52 | 0.026 | 0.013 |
| Exposure | InsDel5a | 883 | 86 (2e-20) | 0.09 | 0.63 | 0.036 | 0.60 | 0.000 | 0.012 |
| Exposure | C\_ID1 | 3,006 | 69 (9e-17) | 0.03 | 0.47 | 0.011 | 0.15 | 0.17 | 0.028 |
| Exposure | InsDel1a | 2,387 | 53 (3e-13) | 0.05 | 0.37 | 0.013 | 0.19 | 0.069 | 0.028 |
| Exposure | C\_ID9 | 742 | 11 (7e-4) | 0.03 | 0.39 | 0.009 | 0.24 | 0.005 | 0.050 |
| Exposure | C\_ID10 | 128 | 21 (4e-6) | 0.26 | 0.45 | 0.093 | 0.28 | 0.000 | 0.000 |
| Exposure | ID\_N | 1,623 | 23 (2e-6) | 0.02 | 0.55 | 0.006 | 0.28 | 0.080 | 0.020 |
| Presence | C\_ID5 | 4,600 | 54 (2e-13) | 0.05 | 0.24 | 0.011 | 0.15 | 0.002 | 0.009 |
| Presence | InsDel5a | 4,600 | 75 (4e-18) | 0.07 | 0.25 | 0.017 | 0.18 | 0.003 | 0.003 |
| Presence | C\_ID1 | 4,600 | 33 (1e-8) | 0.04 | 0.15 | 0.006 | 0.14 | 0.000 | 0.000 |
| Presence | InsDel1a | 4,600 | 39 (3e-10) | 0.05 | 0.14 | 0.006 | 0.09 | 0.008 | 0.018 |
| Presence | C\_ID9 | 4,600 | 24 (1e-6) | 0.06 | 0.11 | 0.006 | 0.07 | 0.001 | 0.006 |
| Presence | C\_ID10 | 4,135 | 0.4 (0.5) | 0.005 | 0.08 | 0.000 | 0.06 | 0.000 | 0.002 |
| Presence | ID\_N | 4,600 | 11 (9e-4) | 0.006 | 0.31 | 0.002 | 0.21 | 0.060 | 0.007 |

Across all signatures, the full model's deviance R² is 0.11–0.84 for exposure (median 0.56) and 0.04–0.72 for presence (median 0.19). Whether a tumor has a signature is therefore mostly not predicted by these covariates.

For the manuscript: "C\_ID5 exposure rises about 19% per decade within cancer types" is accurate. "Age is a major determinant of C\_ID5 activity" would not be.

## Sex, cohort and MSI-H

In the joint model, sex has one effect shared by all cancer types, so it does not answer the per-type question of Table S15. Only a few signatures show a sex effect (male vs female, FDR):

- C\_ID5 exposure is higher in men: 1.19 (6e-3).
- C\_ID3 and InsDel3b presence is more common in men: OR 1.68 (7e-3) and 1.59 (0.03), consistent with the lung/tobacco result.
- C\_ID12 presence is more common in women: OR 0.23 (7e-3), consistent with the new pancreas result in Table S15.

Cohort and MSI-H have large effects and belong in the models. Many signatures differ strongly between HMF metastases and PCAWG primaries (for example InsDel1c presence OR 12 for PCAWG vs HMF). MSI-H explains more of C\_ID1 and ID\_N exposure than age does.

**Separation.** Some presence estimates are 0.00 or very large: cohort for C\_ID19 and InsDel19a/b, and MSI-H for ID\_N and InsDel\_N\_beta. There the signature is present in only one group, so the logistic estimate has no finite value. A published version would need Firth logistic regression (the `logistf` R package) or similar.

**Table S15 is not adjusted for age.** The next section refits it within each cancer type with age, cohort and MSI-H as covariates.

## Table S15 adjusted for age, cohort and MSI-H

All the Table S15 sex associations at q < 0.2 stay below 0.2 after adjustment, with nearly the same odds ratios, so none is an artifact of age. Two new borderline associations appear: CNS C\_ID19/InsDel19b and Lymphoid C\_ID1.

Method: each of the 204 cancer-type × signature groups tested in Table S15 is refit as present \~ sex + age + cohort + MSI-H, by Firth logistic regression (the `logistf` R package), in tumors with known age and sex. Firth regression gives finite estimates when a signature is all present or all absent in one sex. Odds ratios are female vs male, as in Table S15 (> 1 = female-enriched). Covariates with one value in a group are dropped. P values are profile penalized-likelihood tests, BH-adjusted within the 83-type and the 89-type signatures separately, as in Table S15. A variant adds log10 total indels.

| Cancer type | Signature | Tumors (present) | S15 Fisher OR (q) | Adjusted OR (q) | Also adjusted for total indels: OR (q) |
| --- | --- | --- | --- | --- | --- |
| Pancreas | C\_ID12 / InsDel12 | 377 (21) | 8.05 (0.014) | 6.59 (0.021) | 6.53 (0.023) |
| Lung | InsDel3b | 255 (197) | 0.35 (0.031) | 0.37 (0.060) | 0.33 (0.14) |
| Lung | C\_ID3 | 255 (199) | 0.38 (0.081) | 0.40 (0.135) | 0.39 (0.25) |
| Esophagus | InsDel1a | 192 (48) | 3.40 (0.105) | 3.16 (0.157) | 3.13 (0.14) |
| Esophagus | C\_ID1 | 192 (71) | 2.96 (0.182) | 3.08 (0.135) | 3.03 (0.22) |
| Liver | InsDel5a | 355 (150) | 0.50 (0.163) | 0.50 (0.157) | 0.51 (0.14) |
| Liver | C\_ID5 | 355 (171) | 0.53 (0.182) | 0.53 (0.157) | 0.53 (0.22) |
| Liver | InsDel8 / C\_ID8 | 355 (39) | 0.27 (0.16–0.18) | 0.31 (0.157–0.188) | 0.31 (0.18–0.22) |
| CNS | C\_ID19 | 338 (54) | 1.01 (1.0) | 5.34 (0.139) | 4.87 (0.22) |
| CNS | InsDel19b | 338 (53) | 0.96 (1.0) | 4.76 (0.191) | 4.26 (0.24) |
| Lymphoid | C\_ID1 | 215 (64) | 2.02 (0.33) | 2.19 (0.177) | 1.73 (0.61) |
| Bone/soft tissue | C\_ID4 / InsDel4 | 110 (10) | 1.40 (0.81–0.85) | 6.34 (0.19–0.21) | 6.47 (0.20–0.23) |

- **CNS C\_ID19/InsDel19b** appear only after adjustment (OR 1.0 → about 5), so sex is entangled with age or cohort in CNS. This needs a look before it is believed.
- **Bone/soft tissue C\_ID4/InsDel4** is not trustworthy: 99 of 213 tumors have no age, which leaves 10 with the signature, and the Fisher test used a different set of tumors.
- **Lymphoid C\_ID1** is borderline and disappears with total indels added.
- **Adjusting for total indels** weakens Lung most. Smoking raises the indel count, so that adjustment removes part of the effect being tested. The model without total indels should be the main analysis.

**CNS and lymphoid tumors are very heterogeneous.** Each lumps together subtypes that can differ in patient age and sex, for example pediatric and adult brain tumors, or lymphoma subtypes. A sex or age association in these groups could reflect which subtypes the tumors come from rather than an effect within any one subtype, and adjusting for age and cohort does not remove that.

### Subtypes behind the CNS, lymphoid and bone/soft tissue associations

The subtypes explain all three new associations, so none of them should be treated as real. Subtypes come from the PCAWG histology codes (`all.meta.txt`, from `build_all_meta.R` in the data repo) and the HMF `cancerSubtype` (`HMF_metadata.tsv`, controlled access, not in any repo).

**CNS, C\_ID19.** C\_ID19 is in 55 of 63 HMF glioblastomas but none of 41 PCAWG glioblastomas, which have the same median age (60) and a similar sex ratio. So in CNS it marks the HMF cohort, possibly prior treatment. CNS also mixes pediatric tumors that are sex-balanced with adult glioblastomas that are male-skewed. The apparent sex effect after adjustment comes from that mixture.

| CNS subtype | Tumors | Female / male | Median age | With C\_ID19 |
| --- | --- | --- | --- | --- |
| PCAWG medulloblastoma | 146 | 67 / 79 | 9 | 3 |
| PCAWG pilocytic astrocytoma | 89 | 47 / 42 | 8 | 5 |
| HMF glioblastoma | 63 | 15 / 29 | 60 | 55 |
| PCAWG glioblastoma | 41 | 13 / 28 | 60 | 0 |
| PCAWG oligodendroglioma | 18 | 9 / 9 | 40.5 | 0 |

**Lymphoid, C\_ID1.** CLL is two-thirds male and has much less C\_ID1 than B-cell non-Hodgkin lymphoma. Mixing the two gives an apparent female enrichment of C\_ID1.

| Lymphoid subtype | Tumors | Female / male | Median age | With C\_ID1 |
| --- | --- | --- | --- | --- |
| PCAWG B-cell non-Hodgkin lymphoma | 107 | 51 / 56 | 57 | 44 (41%) |
| PCAWG chronic lymphocytic leukemia | 95 | 31 / 64 | 62 | 15 (16%) |
| HMF non-Hodgkin | 21 | 2 / 8 | 71 | 9 |

**Bone/soft tissue, C\_ID4.**

- About 30 subtypes, including GIST, leiomyosarcoma, osteosarcoma, liposarcoma and benign bone tumors.
- None of its 137 PCAWG tumors has an age, so the adjusted model used HMF tumors only.
- The 10 tumors with C\_ID4 are scattered across 9 subtypes.
- Leiomyosarcoma is 18 female / 5 male, which by itself can tilt the sex ratio.

Other coarse groups mix subtypes too, for example lung adenocarcinoma and squamous cell carcinoma. In lung, though, the male excess and the tobacco signature are both caused by smoking, which is a real biological link.

To decide: drop sex and age associations in CNS, lymphoid and bone/soft tissue; use subtype as a covariate where it is known (it is missing from the shared data for HMF); or report them with this caveat.

Script: `build_fig_s6/sex_glm_prototype.R` in the code repo. Results for all 204 groups: `build_fig_s6/sex_glm_prototype.tsv`.

## Questions to decide

- [ ] Use the joint models (age, sex, cohort, MSI-H, cancer type) for Table S16, or keep the age + cancer type models now in the code repo?
- [ ] Exposure as the main test and presence as secondary, or report both equally?
- [ ] Replace the Table S15 Fisher tests with per-type logistic models adjusted for age and cohort?
- [ ] Fit presence models with Firth regression to handle separation?
- [ ] Follow up ID\_N and InsDel\_N\_beta (strong age × cancer type interaction, exposure and presence in opposite directions), or leave them out of the text?
- [ ] Is a 30-tumor minimum per model and 3 present tumors per cancer type right?
- [ ] Report deviance R² or likelihood-ratio shares in the paper, to make clear that age explains little of the variation?
- [ ] Methods paragraph and Table S16 legend: who writes them once the analysis is settled?
