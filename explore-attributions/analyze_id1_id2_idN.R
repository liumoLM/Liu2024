# Focused analysis of the InsDel1a/1b/1c, InsDel2a/2b/2c, and
# InsDel_Nalpha/InsDel_Nbeta splits in relation to:
#   - MSI vs MSS
#   - Hypermutator (total indels >= 5000) vs non-hypermutator
#   - The MMR block (C_ID7, InsDel_D, InsDel_J)
#
# Outputs (all written to this folder):
#   id1_id2_idN_cooccurrence_focused.pdf  (panel A: phi heatmap, focused subset)
#   id1_id2_idN_stripplot_by_stratum.pdf  (panel B: per-sample share by stratum)
#   id1_id2_idN_summary.csv               (per-signature x stratum summary stats)
#   idN_alpha_vs_beta_investigation.csv   (extra: N alpha vs N beta details)
#   idN_alpha_vs_beta_panels.pdf          (extra: N alpha vs N beta visuals)
#
# Inputs:
#   ../Sup Tables/Table S10 83-type and 89-type signature assignment.tsv
#   ../Sup Tables/Table S17 metadata of 6975 samples.xlsx

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

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

PRESENCE_FRAC <- 0.05
HYPER_CUT     <- 5000  # total assigned indels per sample

# ---- Load ----
message("Loading metadata")
meta <- as.data.frame(read_excel(
  file.path(sup_dir, "Table S17 metadata of 6975 samples.xlsx"),
  sheet = "Sheet1"
))
meta$major <- meta[["Major Cancer Type"]]
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
samples   <- colnames(assign_mat)

sample_totals <- colSums(assign_mat)
frac_mat <- sweep(as.matrix(assign_mat), 2, sample_totals, FUN = "/")
frac_mat[is.na(frac_mat)] <- 0
presence <- frac_mat >= PRESENCE_FRAC

# Per-sample meta aligned to sample order
ix <- match(samples, meta$Patient)
samp <- data.frame(
  sample = samples,
  major  = meta$major[ix],
  msi    = meta$msi[ix],
  total  = sample_totals,
  stringsAsFactors = FALSE
)
samp$stratum <- with(samp, ifelse(
  msi == "MSI", "MSI-H",
  ifelse(total >= HYPER_CUT, "MSS hypermutated", "MSS non-hypermutated")
))
samp$stratum <- factor(samp$stratum,
                       levels = c("MSS non-hypermutated", "MSS hypermutated", "MSI-H"))
message("Stratum counts:")
print(table(samp$stratum, useNA = "ifany"))

# ---- Focused signature subset ----
focus_sigs <- c(
  "C_ID1/InsDel1a", "C_ID1/InsDel1b", "C_ID1/InsDel1c",
  "C_ID2/InsDel2a", "C_ID2/InsDel2b", "C_ID2/InsDel2c",
  "ID_N/InsDel_N_alpha", "ID_N/InsDel_N_beta",
  "C_ID7/InsDel7", "ID_D/InsDel_D", "ID_J/InsDel_J"
)
focus_sigs <- focus_sigs[focus_sigs %in% sig_names]

# =============================================================================
# Panel A: focused co-occurrence heatmap (phi)
# =============================================================================
message("Building focused co-occurrence heatmap")

cooc_pair <- function(x, y) {
  a <- as.numeric(sum(x & y)); b <- as.numeric(sum(x & !y))
  c_ <- as.numeric(sum(!x & y)); d <- as.numeric(sum(!x & !y))
  denom <- sqrt((a + b) * (c_ + d) * (a + c_) * (b + d))
  phi <- if (denom > 0) (a * d - b * c_) / denom else NA_real_
  p_val <- tryCatch(
    fisher.test(matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE))$p.value,
    error = function(e) NA_real_
  )
  list(n11 = a, n10 = b, n01 = c_, n00 = d,
       phi = phi, p_value = p_val)
}

pres_focus <- presence[focus_sigs, , drop = FALSE]
rows <- list()
for (i in seq_along(focus_sigs)) {
  for (j in seq_along(focus_sigs)) {
    if (j <= i) next
    r <- cooc_pair(pres_focus[i, ], pres_focus[j, ])
    rows[[length(rows) + 1]] <- data.frame(
      sig1 = focus_sigs[i], sig2 = focus_sigs[j],
      n11 = r$n11, n10 = r$n10, n01 = r$n01, n00 = r$n00,
      phi = r$phi, p_value = r$p_value, stringsAsFactors = FALSE
    )
  }
}
cooc_df <- do.call(rbind, rows)
cooc_df$p_adj <- p.adjust(cooc_df$p_value, method = "BH")
write.csv(cooc_df,
          file.path(this_dir, "id1_id2_idN_cooccurrence_focused.csv"),
          row.names = FALSE)

phi_mat <- matrix(0, length(focus_sigs), length(focus_sigs),
                  dimnames = list(focus_sigs, focus_sigs))
for (k in seq_len(nrow(cooc_df))) {
  phi_mat[cooc_df$sig1[k], cooc_df$sig2[k]] <- cooc_df$phi[k]
  phi_mat[cooc_df$sig2[k], cooc_df$sig1[k]] <- cooc_df$phi[k]
}
diag(phi_mat) <- 1

# Keep manually-defined order (so 1a/1b/1c, 2a/2b/2c, Na/Nb, MMR group)
ord <- focus_sigs
phi_long <- as.data.frame(as.table(phi_mat[ord, ord]))
colnames(phi_long) <- c("sig1", "sig2", "phi")
phi_long$sig1 <- factor(phi_long$sig1, levels = ord)
phi_long$sig2 <- factor(phi_long$sig2, levels = rev(ord))

# Significance asterisks
sig_lookup <- cooc_df
sig_lookup$key <- paste(sig_lookup$sig1, sig_lookup$sig2, sep = "||")
phi_long$key1 <- paste(phi_long$sig1, phi_long$sig2, sep = "||")
phi_long$key2 <- paste(phi_long$sig2, phi_long$sig1, sep = "||")
phi_long$p_adj <- sig_lookup$p_adj[match(phi_long$key1, sig_lookup$key)]
ix2 <- is.na(phi_long$p_adj)
phi_long$p_adj[ix2] <- sig_lookup$p_adj[match(phi_long$key2[ix2], sig_lookup$key)]
phi_long$star <- ifelse(is.na(phi_long$p_adj), "",
                 ifelse(phi_long$p_adj < 1e-20, "***",
                 ifelse(phi_long$p_adj < 1e-10, "**",
                 ifelse(phi_long$p_adj < 0.01,  "*", ""))))

p_heat <- ggplot(phi_long, aes(sig1, sig2, fill = phi)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sprintf("%.2f", phi)),
            size = 2.6, color = "black") +
  geom_text(aes(label = star), nudge_y = -0.28, size = 2.6, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = expression(phi)) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank()) +
  labs(title = paste0("Co-occurrence (phi) of ID1 / ID2 / ID_N children and MMR block, ",
                      "presence = share >= ", PRESENCE_FRAC),
       subtitle = paste0("n = ", ncol(presence), " samples. ",
                         "*: q<0.01  **: q<1e-10  ***: q<1e-20"),
       x = NULL, y = NULL)

pdf(file.path(this_dir, "id1_id2_idN_cooccurrence_focused.pdf"),
    width = 8, height = 7.5)
print(p_heat)
dev.off()
message("Wrote id1_id2_idN_cooccurrence_focused.pdf")

# =============================================================================
# Panel B: per-sample strip plot, focused signatures by stratum
# =============================================================================
message("Building per-sample strip plot")

frac_df <- as.data.frame(t(frac_mat[focus_sigs, , drop = FALSE]),
                         check.names = FALSE)
frac_df$sample <- rownames(frac_df)
frac_long <- pivot_longer(frac_df, -sample,
                          names_to = "signature", values_to = "share")
frac_long <- merge(frac_long, samp[, c("sample", "stratum")], by = "sample")
frac_long$signature <- factor(frac_long$signature, levels = focus_sigs)
frac_long <- frac_long[!is.na(frac_long$stratum), ]

p_strip <- ggplot(frac_long, aes(stratum, share, color = stratum)) +
  geom_jitter(width = 0.30, height = 0,
              size = 0.45, alpha = 0.35) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.55, color = "black", fatten = 1.4) +
  geom_hline(yintercept = PRESENCE_FRAC,
             linetype = "dashed", color = "grey30", linewidth = 0.3) +
  facet_wrap(~ signature, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                                "MSS hypermutated"     = "#ff7f00",
                                "MSI-H"                = "#e31a1c")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1),
        strip.text  = element_text(face = "bold")) +
  labs(title = "Per-sample signature share by MSI / hypermutator stratum",
       subtitle = "Dashed line = 5% presence threshold. Black bar = median.",
       x = NULL, y = "Signature share of sample's assigned indels")

pdf(file.path(this_dir, "id1_id2_idN_stripplot_by_stratum.pdf"),
    width = 9.5, height = 9)
print(p_strip)
dev.off()
message("Wrote id1_id2_idN_stripplot_by_stratum.pdf")

# =============================================================================
# Per-signature x stratum summary statistics
# =============================================================================
summary_df <- frac_long %>%
  group_by(signature, stratum) %>%
  summarise(
    n_samples       = n(),
    n_present_5pct  = sum(share >= PRESENCE_FRAC),
    pct_present     = 100 * n_present_5pct / n_samples,
    median_share    = median(share),
    mean_share      = mean(share),
    max_share       = max(share),
    .groups = "drop"
  ) %>%
  arrange(signature, stratum)

# Also add: total absolute indels attributed to each signature by stratum
counts_focus <- t(assign_mat[focus_sigs, , drop = FALSE])
counts_focus <- as.data.frame(counts_focus, check.names = FALSE)
counts_focus$sample <- rownames(counts_focus)
cf_long <- pivot_longer(counts_focus, -sample,
                        names_to = "signature", values_to = "count")
cf_long <- merge(cf_long, samp[, c("sample", "stratum")], by = "sample")
cf_long <- cf_long[!is.na(cf_long$stratum), ]
totals_df <- cf_long %>%
  group_by(signature, stratum) %>%
  summarise(total_indels_attributed = sum(count), .groups = "drop")

summary_df <- left_join(summary_df, totals_df, by = c("signature", "stratum"))
write.csv(summary_df,
          file.path(this_dir, "id1_id2_idN_summary.csv"),
          row.names = FALSE)
message("Wrote id1_id2_idN_summary.csv")
cat("\n=== Per-stratum summary (focused signatures) ===\n")
print(as.data.frame(summary_df), row.names = FALSE, digits = 3)

# =============================================================================
# InsDel_N alpha vs InsDel_N beta investigation
# =============================================================================
message("Running N_alpha / N_beta investigation")

n_alpha <- "ID_N/InsDel_N_alpha"
n_beta  <- "ID_N/InsDel_N_beta"

# Counts of presence and total attribution, by stratum
n_inv <- summary_df %>%
  dplyr::filter(signature %in% c(n_alpha, n_beta))

# Among samples with N_alpha >= 5%: how many also have N_beta >= 5%?
both_alpha <- presence[n_alpha, ]
both_beta  <- presence[n_beta,  ]
overlap <- data.frame(
  metric = c(
    "samples with N_alpha share >= 5%",
    "samples with N_beta  share >= 5%",
    "samples with both >= 5%",
    "samples with neither >= 5%"
  ),
  count = c(
    sum(both_alpha),
    sum(both_beta),
    sum(both_alpha & both_beta),
    sum(!both_alpha & !both_beta)
  )
)

# By stratum
by_stratum <- data.frame(
  stratum = c("MSS non-hyper", "MSS hyper", "MSI-H"),
  n_alpha_pos = c(
    sum(both_alpha & samp$stratum == "MSS non-hypermutated"),
    sum(both_alpha & samp$stratum == "MSS hypermutated"),
    sum(both_alpha & samp$stratum == "MSI-H")
  ),
  n_beta_pos = c(
    sum(both_beta & samp$stratum == "MSS non-hypermutated"),
    sum(both_beta & samp$stratum == "MSS hypermutated"),
    sum(both_beta & samp$stratum == "MSI-H")
  ),
  n_stratum = c(
    sum(samp$stratum == "MSS non-hypermutated"),
    sum(samp$stratum == "MSS hypermutated"),
    sum(samp$stratum == "MSI-H")
  )
)
by_stratum$pct_alpha <- 100 * by_stratum$n_alpha_pos / by_stratum$n_stratum
by_stratum$pct_beta  <- 100 * by_stratum$n_beta_pos  / by_stratum$n_stratum

# How does the per-sample N_alpha+N_beta total compare to C_ID2 children?
# Specifically: in samples with high N_beta, what is happening to C_ID2 attribution?
nb_frac  <- frac_mat[n_beta, ]
na_frac  <- frac_mat[n_alpha, ]
id2a_frac <- frac_mat["C_ID2/InsDel2a", ]
id2b_frac <- frac_mat["C_ID2/InsDel2b", ]
id2c_frac <- frac_mat["C_ID2/InsDel2c", ]
id2_total_frac <- id2a_frac + id2b_frac + id2c_frac

# Spearman correlations across all samples and within strata
spearman_table <- function(x, y, label) {
  c1 <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
  data.frame(pair = label, rho_spearman = c1)
}
cor_rows <- rbind(
  spearman_table(na_frac, id2_total_frac, "N_alpha vs sum(InsDel2a/b/c)"),
  spearman_table(nb_frac, id2_total_frac, "N_beta  vs sum(InsDel2a/b/c)"),
  spearman_table(na_frac, nb_frac,        "N_alpha vs N_beta"),
  spearman_table(na_frac, frac_mat["C_ID7/InsDel7", ], "N_alpha vs C_ID7"),
  spearman_table(nb_frac, frac_mat["C_ID1/InsDel1a", ], "N_beta  vs InsDel1a")
)

cat("\n=== N_alpha vs N_beta overlap (all samples) ===\n")
print(overlap, row.names = FALSE)
cat("\n=== N_alpha and N_beta presence by stratum ===\n")
print(by_stratum, row.names = FALSE, digits = 3)
cat("\n=== Spearman correlations across all samples ===\n")
print(cor_rows, row.names = FALSE, digits = 3)

write.csv(rbind(
    data.frame(section = "summary_by_stratum",
               item    = paste(n_inv$signature, n_inv$stratum, sep = " | "),
               value   = paste0("present ", n_inv$n_present_5pct, "/",
                                n_inv$n_samples, " (",
                                round(n_inv$pct_present, 1), "%), median share ",
                                signif(n_inv$median_share, 3),
                                ", total indels ",
                                n_inv$total_indels_attributed)),
    data.frame(section = "presence_overlap",
               item    = overlap$metric,
               value   = as.character(overlap$count)),
    data.frame(section = "presence_by_stratum",
               item    = by_stratum$stratum,
               value   = paste0("N_alpha ", by_stratum$n_alpha_pos, "/",
                                by_stratum$n_stratum, " (",
                                round(by_stratum$pct_alpha, 1), "%); ",
                                "N_beta ",  by_stratum$n_beta_pos,  "/",
                                by_stratum$n_stratum, " (",
                                round(by_stratum$pct_beta,  1), "%)")),
    data.frame(section = "spearman_correlations",
               item    = cor_rows$pair,
               value   = signif(cor_rows$rho_spearman, 3))
  ),
  file.path(this_dir, "idN_alpha_vs_beta_investigation.csv"),
  row.names = FALSE
)
message("Wrote idN_alpha_vs_beta_investigation.csv")

# Visuals: 3-panel PDF
# (i) N_alpha share vs N_beta share, color = stratum
df_ab <- data.frame(
  sample = samples,
  N_alpha = na_frac,
  N_beta  = nb_frac,
  ID2_sum = id2_total_frac,
  C_ID7   = frac_mat["C_ID7/InsDel7", ],
  InsDel1a = frac_mat["C_ID1/InsDel1a", ],
  stratum = samp$stratum
)
df_ab <- df_ab[!is.na(df_ab$stratum), ]

p_ab1 <- ggplot(df_ab, aes(N_alpha, N_beta, color = stratum)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = PRESENCE_FRAC, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = PRESENCE_FRAC, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                                "MSS hypermutated"     = "#ff7f00",
                                "MSI-H"                = "#e31a1c")) +
  theme_bw(base_size = 10) +
  labs(title = "InsDel_N_alpha vs InsDel_N_beta share, per sample",
       x = "InsDel_N_alpha share", y = "InsDel_N_beta share")

p_ab2 <- ggplot(df_ab, aes(ID2_sum, N_beta, color = stratum)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                                "MSS hypermutated"     = "#ff7f00",
                                "MSI-H"                = "#e31a1c")) +
  theme_bw(base_size = 10) +
  labs(title = "InsDel_N_beta share vs sum(InsDel2a/b/c) share, per sample",
       x = "Sum of InsDel2a + InsDel2b + InsDel2c share",
       y = "InsDel_N_beta share")

p_ab3 <- ggplot(df_ab, aes(C_ID7, N_alpha, color = stratum)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                                "MSS hypermutated"     = "#ff7f00",
                                "MSI-H"                = "#e31a1c")) +
  theme_bw(base_size = 10) +
  labs(title = "InsDel_N_alpha share vs C_ID7 share, per sample",
       x = "C_ID7/InsDel7 share", y = "InsDel_N_alpha share")

pdf(file.path(this_dir, "idN_alpha_vs_beta_panels.pdf"),
    width = 8, height = 10)
print(p_ab1); print(p_ab2); print(p_ab3)
dev.off()
message("Wrote idN_alpha_vs_beta_panels.pdf")

message("Done.")
