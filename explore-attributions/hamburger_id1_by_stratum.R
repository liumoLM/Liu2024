# Hamburger plots for InsDel1a, InsDel1b, InsDel1c.
# One PDF per signature (3 PDFs), each with 3 stacked rows:
#   row 1: MSS non-hypermutated
#   row 2: MSI-H
#   row 3: MSS hypermutated
# Each row shows per-sample mutation counts attributed to the signature
# across all Major Cancer Types in alphabetical order. Each dot = a sample,
# red dashed line = per-type median, y axis log scale. No double counting:
# each sample appears in exactly one stratum.

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(scales)
  library(dplyr)
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
  sample  = samples,
  major   = meta$major[ix],
  msi     = meta$msi[ix],
  total   = totals,
  stringsAsFactors = FALSE
)
samp$stratum <- with(samp, ifelse(
  msi == "MSI", "MSI-H",
  ifelse(total >= HYPER_CUT, "MSS hypermutated", "MSS non-hypermutated")
))
samp$stratum <- factor(samp$stratum, levels = STRATUM_LEVELS)
samp <- samp[!is.na(samp$major) & !is.na(samp$stratum), ]

all_types <- sort(unique(samp$major))

# stratum sample sizes (constant across signatures)
strat_n <- table(samp$stratum)

plot_one <- function(sig_row, sig_name) {
  vals <- as.numeric(sig_row[samp$sample])
  df <- data.frame(
    sample      = samp$sample,
    cancer_type = factor(samp$major, levels = all_types),
    stratum     = samp$stratum,
    mutations   = vals
  )

  # Per (stratum x type) totals and nonzero counts BEFORE excluding zero
  counts_df <- df %>%
    group_by(stratum, cancer_type, .drop = FALSE) %>%
    summarise(total_samples   = dplyr::n(),
              nonzero_samples = sum(mutations > 0),
              .groups = "drop")

  df_pos <- df[df$mutations > 0, ]

  # Medians per (stratum x type) using nonzero samples
  medians_df <- df_pos %>%
    group_by(stratum, cancer_type, .drop = FALSE) %>%
    summarise(median_mutations = ifelse(dplyr::n() > 0, median(mutations), NA_real_),
              .groups = "drop")

  # x-position per type (alphabetical)
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

  # Background alternating bands
  n_types <- length(all_types)
  bg_df <- data.frame(
    xmin = seq(0.5, n_types - 0.5, by = 1),
    xmax = seq(1.5, n_types + 0.5, by = 1),
    fill = rep(c("#EDF8B1", "#2C7FB8"), length.out = n_types)
  )

  # Y label position (below smallest positive)
  if (nrow(df_pos) > 0) {
    y_min <- min(df_pos$mutations)
    y_label_pos <- max(y_min / 10, 0.1)
  } else {
    y_label_pos <- 0.1
  }

  # Strata facet labels with n
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
         subtitle = "Each dot = a sample. Red dashed line = per-type median. Counts beneath bars: nonzero / total samples in that stratum.",
         x = NULL, y = "Mutations attributed (log scale)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 0.0, vjust = -0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          plot.margin = margin(t = 5, r = 5, b = 60, l = 5, unit = "pt"))
}

sigs <- c("C_ID1/InsDel1a", "C_ID1/InsDel1b", "C_ID1/InsDel1c")
file_tags <- c("InsDel1a", "InsDel1b", "InsDel1c")

for (i in seq_along(sigs)) {
  sig <- sigs[i]; tag <- file_tags[i]
  message("Plotting ", sig)
  p <- plot_one(assign_mat[sig, , drop = TRUE], sig)
  out_pdf <- file.path(this_dir, sprintf("hamburger_%s_by_stratum.pdf", tag))
  pdf(out_pdf, width = 13, height = 11)
  print(p)
  dev.off()
  message("Wrote ", out_pdf)
}
message("Done.")
