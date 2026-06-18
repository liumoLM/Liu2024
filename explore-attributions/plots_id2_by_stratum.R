# Violin + hamburger plots for InsDel2a, InsDel2b, InsDel2c.
# All plots in one PDF (4 pages):
#   page 1: violin plot of 2a/2b/2c per-sample share by Major Cancer Type x stratum
#   page 2: hamburger for InsDel2a, 3 rows = strata, x = Major Cancer Type (alpha)
#   page 3: same for InsDel2b
#   page 4: same for InsDel2c
# Strata (no double counting):
#   MSS non-hypermutated, MSI-H, MSS hypermutated (>= 5000 total indels).

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(scales)
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

HYPER_CUT <- 5000
STRATUM_LEVELS <- c("MSS non-hypermutated", "MSI-H", "MSS hypermutated")
STRATUM_COLORS <- c("MSS non-hypermutated" = "#1f78b4",
                    "MSI-H"                = "#e31a1c",
                    "MSS hypermutated"     = "#ff7f00")

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
frac    <- sweep(as.matrix(assign_mat), 2, totals, FUN = "/")
frac[is.na(frac)] <- 0

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
samp$stratum <- factor(samp$stratum, levels = STRATUM_LEVELS)
samp <- samp[!is.na(samp$major) & !is.na(samp$stratum), ]
all_types <- sort(unique(samp$major))
strat_n   <- table(samp$stratum)

sigs <- c("C_ID2/InsDel2a", "C_ID2/InsDel2b", "C_ID2/InsDel2c")
tags <- c("InsDel2a", "InsDel2b", "InsDel2c")

# ---- Page 1: violin ----
df_sig <- as.data.frame(t(frac[sigs, , drop = FALSE]), check.names = FALSE)
df_sig$sample <- rownames(df_sig)
df_long <- pivot_longer(df_sig, -sample,
                        names_to = "signature", values_to = "share")
df_long <- merge(df_long, samp[, c("sample", "major", "stratum")], by = "sample")
df_long$signature <- factor(df_long$signature, levels = sigs)
df_long$major <- factor(df_long$major, levels = all_types)

p_violin <- ggplot(df_long, aes(major, share, fill = stratum, color = stratum)) +
  geom_violin(position = position_dodge(width = 0.85),
              width = 0.85, scale = "width",
              alpha = 0.55, linewidth = 0.25,
              trim = TRUE, na.rm = TRUE) +
  stat_summary(aes(group = stratum),
               fun = median, geom = "point",
               position = position_dodge(width = 0.85),
               shape = 95, size = 4, color = "black") +
  facet_grid(signature ~ ., switch = "y") +
  scale_fill_manual(values = STRATUM_COLORS, name = NULL,
                    breaks = STRATUM_LEVELS) +
  scale_color_manual(values = STRATUM_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 9) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  labs(title = "InsDel2a / 2b / 2c per-sample share by Major Cancer Type and stratum",
       subtitle = paste0("Each sample appears in exactly one stratum. ",
                         "Hypermutator cutoff = ", HYPER_CUT, " total indels. ",
                         "Black bar = median."),
       x = NULL, y = "Signature share of sample's assigned indels")

# ---- Hamburger helper ----
plot_hamburger <- function(sig_name) {
  vals <- as.numeric(assign_mat[sig_name, samp$sample])
  df <- data.frame(
    sample      = samp$sample,
    cancer_type = factor(samp$major, levels = all_types),
    stratum     = samp$stratum,
    mutations   = vals
  )
  counts_df <- df %>%
    group_by(stratum, cancer_type, .drop = FALSE) %>%
    summarise(total_samples   = dplyr::n(),
              nonzero_samples = sum(mutations > 0),
              .groups = "drop")
  df_pos <- df[df$mutations > 0, ]
  medians_df <- df_pos %>%
    group_by(stratum, cancer_type, .drop = FALSE) %>%
    summarise(median_mutations = ifelse(dplyr::n() > 0, median(mutations), NA_real_),
              .groups = "drop")
  df_pos$x_pos <- as.numeric(df_pos$cancer_type)
  df_pos <- df_pos[order(df_pos$stratum, df_pos$cancer_type, df_pos$mutations), ]
  df_pos$within_rank <- ave(seq_len(nrow(df_pos)),
                            paste(df_pos$stratum, df_pos$cancer_type),
                            FUN = seq_along)
  df_pos$group_size <- ave(seq_len(nrow(df_pos)),
                           paste(df_pos$stratum, df_pos$cancer_type),
                           FUN = length)
  df_pos$x_plot <- df_pos$x_pos +
    0.8 * (df_pos$within_rank - 1) / pmax(df_pos$group_size - 1, 1) - 0.4
  df_pos$x_plot[df_pos$group_size == 1] <- df_pos$x_pos[df_pos$group_size == 1]
  medians_df$x_pos <- as.numeric(medians_df$cancer_type)
  counts_df$x_pos  <- as.numeric(counts_df$cancer_type)
  n_types <- length(all_types)
  bg_df <- data.frame(
    xmin = seq(0.5, n_types - 0.5, by = 1),
    xmax = seq(1.5, n_types + 0.5, by = 1),
    fill = rep(c("#EDF8B1", "#2C7FB8"), length.out = n_types)
  )
  if (nrow(df_pos) > 0) {
    y_label_pos <- max(min(df_pos$mutations) / 10, 0.1)
  } else {
    y_label_pos <- 0.1
  }
  strat_labels <- setNames(
    paste0(STRATUM_LEVELS, " (n=", strat_n[STRATUM_LEVELS], ")"),
    STRATUM_LEVELS
  )
  ggplot(df_pos, aes(x = x_plot, y = mutations)) +
    geom_rect(data = bg_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.2) +
    scale_fill_identity() +
    geom_point(size = 1.0, alpha = 0.55, color = "black", shape = 16) +
    geom_segment(data = medians_df[!is.na(medians_df$median_mutations), ],
                 aes(x = x_pos - 0.4, xend = x_pos + 0.4,
                     y = median_mutations, yend = median_mutations),
                 color = "red", linewidth = 0.8, linetype = "dashed",
                 inherit.aes = FALSE) +
    geom_text(data = counts_df,
              aes(x = x_pos, y = y_label_pos,
                  label = paste0(nonzero_samples, "/", total_samples)),
              vjust = 2.0, size = 2.4, lineheight = 0.9,
              inherit.aes = FALSE) +
    facet_wrap(~ stratum, ncol = 1, scales = "free_y",
               labeller = labeller(stratum = strat_labels)) +
    scale_x_continuous(breaks = seq_along(all_types),
                       labels = all_types,
                       position = "top",
                       expand = expansion(add = 0.5)) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(sig_name, " mutations attributed per sample"),
         subtitle = "Each dot = a sample. Red dashed line = per-type median over nonzero samples. Counts beneath bars: nonzero / total samples in that stratum.",
         x = NULL, y = "Mutations attributed (log scale)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 0.0, vjust = -0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          plot.margin = margin(t = 5, r = 5, b = 60, l = 5, unit = "pt"))
}

out_pdf <- file.path(this_dir, "plots_InsDel2_by_stratum.pdf")
pdf(out_pdf, width = 13, height = 11)
print(p_violin)
for (s in sigs) {
  message("Plotting ", s)
  suppressWarnings(print(plot_hamburger(s)))
}
dev.off()
message("Wrote ", out_pdf)
