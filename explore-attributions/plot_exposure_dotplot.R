# Reproduce the Figure 4a/4b style "active-proportion x median-exposure"
# dot plot from Table S10 + Table S17, and additionally split samples into
# (a) MSI, (b) MSS-hypermutated, (c) MSS-non-hypermutated subsets so that
# rare exposures aren't drowned out by MSI/hypermutator samples.
#
# Inputs (same folder):
#   Table S10 83-type and 89-type signature assignment.tsv
#   Table S17 metadata of 6975 samples.xlsx
#
# Output (same folder):
#   Figure4_dotplot_panels.pdf
#     Page 1: All samples
#     Page 2: MSI samples only
#     Page 3: MSS hypermutated samples only (total indels >= 5000, excludes MSI)
#     Page 4: MSS non-hypermutated samples (the "background")
#
# Visual encoding (matches Figure 4):
#   row = signature (rows include both 83-type and 89-type rows, with the
#         83-type total computed by summing the 89-type joint rows that share
#         the same 83-type prefix)
#   column = Major Cancer Type
#   size = proportion of samples in that cancer type with the signature
#          present (exposure > 0)
#   color = median mutations per Mb among samples with the signature present,
#           log10 color scale.

library(readxl)
library(ggplot2)
library(scales)

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

GENOME_MB <- 3000
HYPERMUT_INDEL_THRESHOLD <- 5000  # matches the high-TMB cutoff in Methods

# ---- Load ----
meta <- as.data.frame(read_excel(
  file.path(sup_dir, "Table S17 metadata of 6975 samples.xlsx"),
  sheet = "Sheet1"
))
meta$major <- meta[["Major Cancer Type"]]
meta$msi <- ifelse(meta$MSI_status %in% c("MSI", "MSI-H"), "MSI",
                   ifelse(meta$MSI_status == "MSS", "MSS", NA_character_))

assign_mat <- as.matrix(read.delim(
  file.path(sup_dir, "Table S10 83-type and 89-type signature assignment.tsv"),
  row.names = 1,
  check.names = FALSE
))
joint_sigs <- rownames(assign_mat)

# Build 83-type aggregated rows
groups <- split(joint_sigs, sub("/.*", "", joint_sigs))
mat83 <- t(sapply(groups, function(rs) colSums(assign_mat[rs, , drop = FALSE])))
rownames(mat83) <- names(groups)

# Build 89-type-only rows (strip the "X/" prefix)
mat89 <- assign_mat
rownames(mat89) <- sub("^.*/", "", joint_sigs)

# Combine for plotting. Interleave 83 row followed by its 89 children, to
# match the layout in Figure 4.
ordered_rows <- list()
for (g in names(groups)) {
  ordered_rows[[length(ordered_rows) + 1]] <- list(
    name = g, vec = mat83[g, ]
  )
  children <- sub("^.*/", "", groups[[g]])
  for (ch in children) {
    ordered_rows[[length(ordered_rows) + 1]] <- list(
      name = ch, vec = mat89[ch, ]
    )
  }
}
row_names <- vapply(ordered_rows, function(x) x$name, character(1))
plot_mat <- do.call(rbind, lapply(ordered_rows, function(x) x$vec))
rownames(plot_mat) <- row_names
samples <- colnames(plot_mat)

# Map sample -> major, msi, hypermut
samp_meta <- meta[match(samples, meta$Patient), c("major", "msi")]
rownames(samp_meta) <- samples
samp_meta$total <- colSums(assign_mat)
samp_meta$hyper <- samp_meta$total >= HYPERMUT_INDEL_THRESHOLD

# ---- Plotting helper ----
build_dotplot_df <- function(plot_mat, samp_meta, sample_subset) {
  pm <- plot_mat[, sample_subset, drop = FALSE]
  sm <- samp_meta[sample_subset, , drop = FALSE]
  majors <- sort(unique(na.omit(sm$major)))
  rows <- list()
  for (sig in rownames(pm)) {
    for (mj in majors) {
      ix <- which(sm$major == mj)
      if (length(ix) == 0) next
      vals <- pm[sig, ix]
      n_total <- length(vals)
      pos <- vals > 0
      n_pos <- sum(pos)
      if (n_pos == 0) {
        rows[[length(rows) + 1]] <- data.frame(
          signature = sig, major = mj,
          active_prop = 0, median_per_mb = NA_real_,
          n_pos = 0, n_total = n_total,
          stringsAsFactors = FALSE
        )
        next
      }
      med <- median(vals[pos]) / GENOME_MB
      rows[[length(rows) + 1]] <- data.frame(
        signature = sig, major = mj,
        active_prop = n_pos / n_total,
        median_per_mb = med,
        n_pos = n_pos, n_total = n_total,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)
  df$signature <- factor(df$signature, levels = rev(rownames(plot_mat)))
  df$major <- factor(df$major, levels = majors)
  df
}

dotplot <- function(df, title) {
  ggplot(df[df$active_prop > 0, ],
         aes(x = major, y = signature,
             size = active_prop, color = median_per_mb)) +
    geom_point() +
    scale_size_area(
      max_size = 6,
      breaks = c(0.05, 0.25, 0.5, 0.75, 1.0),
      labels = scales::percent,
      name = "Active proportion"
    ) +
    scale_color_gradientn(
      colors = c("#3B0F70", "#8C2981", "#DE4968",
                 "#FE9F6D", "#FCFDBF"),
      trans = "log10",
      labels = scales::label_number(),
      name = "Median mut/Mb\n(in present samples)"
    ) +
    scale_x_discrete(position = "top") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
      panel.grid = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

# ---- Subsets ----
all_samples <- samples
msi_samples <- samples[!is.na(samp_meta$msi) & samp_meta$msi == "MSI"]
mss_hyper_samples <- samples[
  !is.na(samp_meta$msi) & samp_meta$msi == "MSS" & samp_meta$hyper
]
mss_lo_samples <- samples[
  !is.na(samp_meta$msi) & samp_meta$msi == "MSS" & !samp_meta$hyper
]

subsets <- list(
  list(name = "All samples", samp = all_samples),
  list(name = sprintf("MSI samples (n=%d)", length(msi_samples)),
       samp = msi_samples),
  list(name = sprintf(
    "MSS hypermutated samples, total indels >= %d (n=%d)",
    HYPERMUT_INDEL_THRESHOLD, length(mss_hyper_samples)
  ), samp = mss_hyper_samples),
  list(name = sprintf(
    "MSS non-hypermutated samples (n=%d)", length(mss_lo_samples)
  ), samp = mss_lo_samples)
)

out_pdf <- file.path(this_dir, "Figure4_dotplot_panels.pdf")
pdf(out_pdf, width = 11, height = 13)
for (s in subsets) {
  message("Building: ", s$name)
  if (length(s$samp) < 5) {
    message("  skipping, too few samples")
    next
  }
  df <- build_dotplot_df(plot_mat, samp_meta, s$samp)
  print(dotplot(df, s$name))
}
dev.off()
message("Wrote: ", out_pdf)
