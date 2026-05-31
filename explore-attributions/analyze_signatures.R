# Analyses #1 (MSI vs MSS per major cancer type) and #2 (signature pairwise
# co-occurrence) on Table S10 + Table S17.
#
# Inputs (same folder as this script):
#   Table S10 83-type and 89-type signature assignment.tsv
#   Table S17 metadata of 6975 samples.xlsx
#
# Outputs (same folder):
#   analysis1_msi_vs_mss.csv       per-signature, per-Major-Cancer-Type MSI vs MSS Fisher test
#   analysis2_cooccurrence_global.csv  pairwise co-occurrence across all samples
#   analysis2_cooccurrence_global_phi.pdf  heatmap of phi coefficients (global)
#   analysis2_cooccurrence_phi_by_tissue.pdf  phi heatmaps faceted by Major Cancer Type
#
# A sample is treated as "exposed" to a signature if that signature's share
# of the sample's total assigned indels is >= PRESENCE_FRACTION (default 0.05).
# This is more robust to noise than exposure > 0.

library(readxl)
library(ggplot2)

# Resolve script directory robustly
this_dir <- local({
  f <- NULL
  try(f <- sys.frame(1)$ofile, silent = TRUE)
  if (is.null(f)) {
    a <- commandArgs(trailingOnly = FALSE)
    m <- a[grep("^--file=", a)]
    if (length(m) > 0) f <- sub("^--file=", "", m)
  }
  if (is.null(f)) return(getwd())
  f <- gsub("~\\+~", " ", f)
  dirname(normalizePath(f))
})

sup_dir <- file.path(this_dir, "..", "Sup Tables")

PRESENCE_FRACTION <- 0.05
MIN_GROUP_N <- 5  # min positives + negatives in each contingency cell row total

# ---- Load data ----
message("Loading metadata")
meta <- as.data.frame(read_excel(
  file.path(sup_dir, "Table S17 metadata of 6975 samples.xlsx"),
  sheet = "Sheet1"
))
stopifnot(all(c("Patient", "Major Cancer Type", "MSI_status") %in% colnames(meta)))
meta$major <- meta[["Major Cancer Type"]]
# Collapse MSI_status: MSI-H/MSI -> MSI, MSS -> MSS, else NA
meta$msi <- ifelse(
  meta$MSI_status %in% c("MSI-H", "MSI"), "MSI",
  ifelse(meta$MSI_status == "MSS", "MSS", NA_character_)
)

message("Loading assignment matrix")
assign_mat <- read.delim(
  file.path(sup_dir, "Table S10 83-type and 89-type signature assignment.tsv"),
  row.names = 1,
  check.names = FALSE
)
sig_names <- rownames(assign_mat)
samples <- colnames(assign_mat)

# Build per-sample fraction matrix (signature share of that sample's indels)
sample_totals <- colSums(assign_mat)
frac_mat <- sweep(as.matrix(assign_mat), 2, sample_totals, FUN = "/")
frac_mat[is.na(frac_mat)] <- 0
# Presence matrix: TRUE if signature share >= PRESENCE_FRACTION
presence <- frac_mat >= PRESENCE_FRACTION  # rows=sigs, cols=samples

# Map sample -> major cancer type, msi
samp_meta <- meta[match(samples, meta$Patient), c("major", "msi")]
rownames(samp_meta) <- samples

# =============================================================================
# Analysis #1: MSI vs MSS Fisher within each Major Cancer Type
# =============================================================================
message("Running analysis #1: MSI vs MSS per Major Cancer Type")

majors <- sort(unique(na.omit(samp_meta$major)))
res1 <- list()
for (mj in majors) {
  in_major <- which(samp_meta$major == mj & samp_meta$msi %in% c("MSI", "MSS"))
  if (length(in_major) < MIN_GROUP_N) next
  msi_flag <- samp_meta$msi[in_major] == "MSI"
  n_msi <- sum(msi_flag); n_mss <- sum(!msi_flag)
  if (n_msi < 2 || n_mss < 2) next
  for (sig in sig_names) {
    pres <- presence[sig, in_major]
    a <- sum(pres & msi_flag)     # MSI present
    b <- sum(!pres & msi_flag)    # MSI absent
    c_ <- sum(pres & !msi_flag)   # MSS present
    d <- sum(!pres & !msi_flag)   # MSS absent
    if (a + c_ == 0) next  # no positives at all
    tab <- matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE)
    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) next
    res1[[length(res1) + 1]] <- data.frame(
      signature = sig,
      major_cancer_type = mj,
      n_msi = n_msi,
      n_mss = n_mss,
      msi_pos = a,
      msi_neg = b,
      mss_pos = c_,
      mss_neg = d,
      msi_pos_rate = a / n_msi,
      mss_pos_rate = c_ / n_mss,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}
res1_df <- do.call(rbind, res1)
res1_df$p_adj <- p.adjust(res1_df$p_value, method = "BH")
res1_df <- res1_df[order(res1_df$p_adj, res1_df$p_value), ]

out1 <- file.path(this_dir, "analysis1_msi_vs_mss.csv")
write.csv(res1_df, out1, row.names = FALSE)
message("Wrote: ", out1, " (", nrow(res1_df), " rows)")

# Print top hits
top <- head(res1_df[res1_df$p_adj < 0.05, ], 20)
if (nrow(top) > 0) {
  message("Top MSI vs MSS hits (q < 0.05):")
  print(top[, c("signature", "major_cancer_type", "msi_pos_rate",
                "mss_pos_rate", "odds_ratio", "p_adj")], row.names = FALSE)
}

# =============================================================================
# Analysis #2: Pairwise signature co-occurrence
# =============================================================================
message("Running analysis #2: pairwise signature co-occurrence (global)")

cooc <- function(pres_mat, sigs) {
  # pres_mat: rows = sigs, cols = samples (logical)
  out <- list()
  n_samp <- ncol(pres_mat)
  for (i in seq_along(sigs)) {
    for (j in seq_along(sigs)) {
      if (j <= i) next
      x <- pres_mat[sigs[i], ]
      y <- pres_mat[sigs[j], ]
      a <- as.numeric(sum(x & y))
      b <- as.numeric(sum(x & !y))
      c_ <- as.numeric(sum(!x & y))
      d <- as.numeric(sum(!x & !y))
      n <- a + b + c_ + d
      # phi coefficient
      denom <- sqrt((a + b) * (c_ + d) * (a + c_) * (b + d))
      phi <- if (denom > 0) (a * d - b * c_) / denom else NA_real_
      # Haldane-corrected log OR
      log_or <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c_ + 0.5)))
      ft <- tryCatch(
        fisher.test(matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE)),
        error = function(e) NULL
      )
      p_val <- if (is.null(ft)) NA_real_ else ft$p.value
      out[[length(out) + 1]] <- data.frame(
        sig1 = sigs[i], sig2 = sigs[j],
        n11 = a, n10 = b, n01 = c_, n00 = d, n = n,
        phi = phi, log_OR = log_or, p_value = p_val,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

cooc_global <- cooc(presence, sig_names)
cooc_global$p_adj <- p.adjust(cooc_global$p_value, method = "BH")
cooc_global <- cooc_global[order(cooc_global$p_adj, -abs(cooc_global$phi)), ]

out2 <- file.path(this_dir, "analysis2_cooccurrence_global.csv")
write.csv(cooc_global, out2, row.names = FALSE)
message("Wrote: ", out2, " (", nrow(cooc_global), " pairs)")

# ---- Heatmap of global phi ----
phi_mat <- matrix(0, length(sig_names), length(sig_names),
                  dimnames = list(sig_names, sig_names))
for (k in seq_len(nrow(cooc_global))) {
  i <- cooc_global$sig1[k]; j <- cooc_global$sig2[k]
  phi_mat[i, j] <- cooc_global$phi[k]
  phi_mat[j, i] <- cooc_global$phi[k]
}
diag(phi_mat) <- 1

# Hierarchical cluster by phi distance
d <- as.dist(1 - phi_mat)
hc <- hclust(d, method = "average")
ord <- hc$order
phi_df <- data.frame(
  sig1 = rep(sig_names[ord], each = length(sig_names)),
  sig2 = rep(sig_names[ord], times = length(sig_names)),
  phi = as.vector(phi_mat[ord, ord])
)
phi_df$sig1 <- factor(phi_df$sig1, levels = sig_names[ord])
phi_df$sig2 <- factor(phi_df$sig2, levels = sig_names[ord])

p_global <- ggplot(phi_df, aes(sig1, sig2, fill = phi)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1),
    name = expression(phi)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0(
      "Signature co-occurrence (phi), global; presence = share >= ",
      PRESENCE_FRACTION
    ),
    x = NULL, y = NULL
  )

out2_pdf <- file.path(this_dir, "analysis2_cooccurrence_global_phi.pdf")
pdf(out2_pdf, width = 11, height = 11)
print(p_global)
dev.off()
message("Wrote: ", out2_pdf)

# =============================================================================
# Analysis #2b: phi heatmaps per Major Cancer Type
# =============================================================================
message("Running per-tissue co-occurrence heatmaps")

per_tissue_pdf <- file.path(this_dir, "analysis2_cooccurrence_phi_by_tissue.pdf")
pdf(per_tissue_pdf, width = 11, height = 11)
for (mj in majors) {
  ix <- which(samp_meta$major == mj)
  if (length(ix) < 20) next
  pres_t <- presence[, ix, drop = FALSE]
  # Drop signatures that never appear in this tissue
  sig_keep <- rownames(pres_t)[rowSums(pres_t) > 0]
  if (length(sig_keep) < 3) next
  pres_t <- pres_t[sig_keep, , drop = FALSE]
  ct <- cooc(pres_t, sig_keep)
  if (nrow(ct) == 0) next
  ct$p_adj <- p.adjust(ct$p_value, method = "BH")
  phi_mt <- matrix(0, length(sig_keep), length(sig_keep),
                   dimnames = list(sig_keep, sig_keep))
  for (k in seq_len(nrow(ct))) {
    phi_mt[ct$sig1[k], ct$sig2[k]] <- ct$phi[k]
    phi_mt[ct$sig2[k], ct$sig1[k]] <- ct$phi[k]
  }
  phi_mt[is.na(phi_mt)] <- 0
  diag(phi_mt) <- 1
  hc_t <- tryCatch(
    hclust(as.dist(1 - phi_mt), method = "average"),
    error = function(e) NULL
  )
  if (is.null(hc_t)) next
  ord_t <- hc_t$order
  df_t <- data.frame(
    sig1 = rep(sig_keep[ord_t], each = length(sig_keep)),
    sig2 = rep(sig_keep[ord_t], times = length(sig_keep)),
    phi = as.vector(phi_mt[ord_t, ord_t])
  )
  df_t$sig1 <- factor(df_t$sig1, levels = sig_keep[ord_t])
  df_t$sig2 <- factor(df_t$sig2, levels = sig_keep[ord_t])
  p_t <- ggplot(df_t, aes(sig1, sig2, fill = phi)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1),
      name = expression(phi)
    ) +
    coord_fixed() +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0(mj, " (n=", length(ix), ")"),
      x = NULL, y = NULL
    )
  print(p_t)
}
dev.off()
message("Wrote: ", per_tissue_pdf)

message("Done.")
