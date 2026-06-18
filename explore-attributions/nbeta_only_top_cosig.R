# For the "N_beta only" tumors (defined as samples with InsDel_N_beta
# attribution and no attribution to any of InsDel1a/1b/1c or InsDel2a/2b/2c),
# rank the next most common signatures (other than InsDel_N_beta itself)
# by number of N_beta-only samples in which they are present (attribution > 0).
# Then for the top-3 such signatures, produce a count table with the same
# column structure as nbeta_id1_id2_counts.xlsx: Major Cancer Types
# (MSS non-hypermutated samples only), MSI-H, MSS hypermutated, Row total.
#
# The output xlsx contains TWO sheets:
#   "ranked_signatures": full ranking of signatures among the 1195 N_beta-only samples
#   "top3_counts": cancer-type x stratum count table for the top 3 co-signatures

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tidyr)
})

this_dir <- "/home/steve/github/0Liu2024/explore-attributions"
sup_dir  <- file.path(this_dir, "..", "Sup Tables")

HYPER_CUT <- 5000

meta <- as.data.frame(read_excel(
  file.path(sup_dir, "Table S17 metadata of 6975 samples.xlsx"),
  sheet = "Sheet1"
))
meta$major <- meta[["Major Cancer Type"]]
meta$msi <- ifelse(
  meta$MSI_status %in% c("MSI-H", "MSI"), "MSI",
  ifelse(meta$MSI_status == "MSS", "MSS", NA_character_)
)

assign_mat <- read.delim(
  file.path(sup_dir, "Table S10 83-type and 89-type signature assignment.tsv"),
  row.names = 1,
  check.names = FALSE
)
samples <- colnames(assign_mat)
totals  <- colSums(assign_mat)

ix <- match(samples, meta$Patient)
samp <- data.frame(
  sample = samples,
  major  = meta$major[ix],
  msi    = meta$msi[ix],
  total  = totals,
  stringsAsFactors = FALSE
)
samp$stratum <- with(samp, ifelse(
  msi == "MSI", "MSI-H",
  ifelse(total >= HYPER_CUT, "MSS hypermutated", "MSS non-hypermutated")
))

id1_sigs <- c("C_ID1/InsDel1a", "C_ID1/InsDel1b", "C_ID1/InsDel1c")
id2_sigs <- c("C_ID2/InsDel2a", "C_ID2/InsDel2b", "C_ID2/InsDel2c")
n_beta   <- "ID_N/InsDel_N_beta"

pres <- as.matrix(assign_mat) > 0
has_nbeta <- pres[n_beta, samp$sample]
n_id1 <- colSums(pres[id1_sigs, samp$sample, drop = FALSE])
n_id2 <- colSums(pres[id2_sigs, samp$sample, drop = FALSE])

nbeta_only_idx <- which(has_nbeta & n_id1 == 0 & n_id2 == 0
                        & !is.na(samp$major) & !is.na(samp$stratum))
cat("N_beta-only sample count:", length(nbeta_only_idx), "\n\n")

# Rank all signatures except N_beta and the InsDel1/InsDel2 children by
# number of N_beta-only samples in which they are present
exclude <- c(n_beta, id1_sigs, id2_sigs)
candidate_sigs <- setdiff(rownames(pres), exclude)
pres_sub <- pres[candidate_sigs, samp$sample[nbeta_only_idx], drop = FALSE]
counts <- rowSums(pres_sub)
ranked <- data.frame(
  signature = names(counts),
  n_nbeta_only_samples_present = unname(counts),
  stringsAsFactors = FALSE
)
ranked$pct_of_nbeta_only <- 100 * ranked$n_nbeta_only_samples_present / length(nbeta_only_idx)
ranked <- ranked[order(-ranked$n_nbeta_only_samples_present), ]
cat("Top 10 signatures co-present in N_beta-only samples:\n")
print(head(ranked, 10), row.names = FALSE, digits = 4)

top3 <- ranked$signature[1:3]
cat("\nTop 3 co-signatures (rank 2, 3, 4 overall):\n", paste(top3, collapse = ", "), "\n")

# Build cancer-type x stratum count table for each top-3 co-signature
type_order <- sort(unique(samp$major[!is.na(samp$major)]))

build_row <- function(sig) {
  ok <- nbeta_only_idx[pres[sig, samp$sample[nbeta_only_idx]]]
  sub <- samp[ok, ]
  row <- data.frame(category = paste0("N_beta only AND ", sig, " present"),
                    stringsAsFactors = FALSE)
  mss_nh_counts <- table(sub$major[sub$stratum == "MSS non-hypermutated"])
  for (ct in type_order) {
    row[[ct]] <- as.integer(mss_nh_counts[ct])
    if (is.na(row[[ct]])) row[[ct]] <- 0L
  }
  row$`MSI-H` <- sum(sub$stratum == "MSI-H")
  row$`MSS hypermutated` <- sum(sub$stratum == "MSS hypermutated")
  row$`Row total` <- sum(row[, -1, drop = FALSE])
  row
}

# Also include the original "N_beta only" baseline row (= 1195 total)
baseline_sub <- samp[nbeta_only_idx, ]
baseline <- data.frame(category = "N_beta only (any other sig allowed except InsDel1/2)",
                       stringsAsFactors = FALSE)
mss_nh_counts <- table(baseline_sub$major[baseline_sub$stratum == "MSS non-hypermutated"])
for (ct in type_order) {
  baseline[[ct]] <- as.integer(mss_nh_counts[ct])
  if (is.na(baseline[[ct]])) baseline[[ct]] <- 0L
}
baseline$`MSI-H` <- sum(baseline_sub$stratum == "MSI-H")
baseline$`MSS hypermutated` <- sum(baseline_sub$stratum == "MSS hypermutated")
baseline$`Row total` <- sum(baseline[, -1, drop = FALSE])

top3_table <- rbind(baseline, build_row(top3[1]), build_row(top3[2]), build_row(top3[3]))

cat("\n=== Top-3 co-signature counts among N_beta-only samples ===\n")
print(top3_table, row.names = FALSE)

out_xlsx <- file.path(this_dir, "nbeta_only_top_cosig.xlsx")
write_xlsx(list(ranked_signatures = ranked,
                top3_counts       = top3_table),
           out_xlsx)
message("Wrote ", out_xlsx)
