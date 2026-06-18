# Same as nbeta_id1_id2_counts.R but with presence threshold = 0
# (any nonzero attribution counts as "present"). Writes an .xlsx file.

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

# Threshold 0: any nonzero attribution counts as present
pres <- as.matrix(assign_mat) > 0
samp$has_nbeta <- pres[n_beta, samp$sample]
samp$n_id1 <- colSums(pres[id1_sigs, samp$sample, drop = FALSE])
samp$n_id2 <- colSums(pres[id2_sigs, samp$sample, drop = FALSE])

samp$category <- with(samp, ifelse(
  !has_nbeta, NA_character_,
  ifelse(n_id1 == 0 & n_id2 == 0, "N_beta only",
  ifelse(n_id1 >= 1 & n_id2 == 0, "N_beta + >=1 InsDel1, 0 InsDel2",
  ifelse(n_id1 == 0 & n_id2 >= 1, "N_beta + 0 InsDel1, >=1 InsDel2",
                                  "N_beta + >=1 InsDel1, >=1 InsDel2")))
))

rows_order <- c("N_beta only",
                "N_beta + >=1 InsDel1, 0 InsDel2",
                "N_beta + 0 InsDel1, >=1 InsDel2",
                "N_beta + >=1 InsDel1, >=1 InsDel2")

type_order <- sort(unique(samp$major[!is.na(samp$major)]))

with_cat <- samp[!is.na(samp$category) & !is.na(samp$major), ]
mss_nh   <- with_cat[with_cat$stratum == "MSS non-hypermutated", ]
ct_tab <- mss_nh %>%
  group_by(category, major) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  pivot_wider(names_from = major, values_from = n, values_fill = 0)
msi_col <- with_cat %>%
  dplyr::filter(stratum == "MSI-H") %>%
  group_by(category) %>%
  summarise(`MSI-H` = dplyr::n(), .groups = "drop")
hyper_col <- with_cat %>%
  dplyr::filter(stratum == "MSS hypermutated") %>%
  group_by(category) %>%
  summarise(`MSS hypermutated` = dplyr::n(), .groups = "drop")

full <- data.frame(category = rows_order, stringsAsFactors = FALSE)
for (ct in type_order) {
  if (ct %in% colnames(ct_tab)) {
    full[[ct]] <- ct_tab[[ct]][match(rows_order, ct_tab$category)]
  } else {
    full[[ct]] <- 0
  }
}
full[is.na(full)] <- 0
full$`MSI-H` <- msi_col$`MSI-H`[match(rows_order, msi_col$category)]
full$`MSS hypermutated` <- hyper_col$`MSS hypermutated`[match(rows_order, hyper_col$category)]
full[is.na(full)] <- 0
full$`Row total` <- rowSums(full[, -1, drop = FALSE])

out_xlsx <- file.path(this_dir, "nbeta_id1_id2_counts_thr0.xlsx")
write_xlsx(full, out_xlsx)
message("Wrote ", out_xlsx)
cat("\nThreshold = 0 (any nonzero attribution counts as present)\n\n")
print(full, row.names = FALSE)
