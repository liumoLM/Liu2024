# Violin plots of InsDel1a, InsDel1b, InsDel1c per-sample share, by
# Major Cancer Type and by MSI / hypermutator stratum. No sample is
# double-counted: each sample falls into exactly one of three strata:
#   MSS non-hypermutated, MSS hypermutated (>= 5000 total indels), MSI-H.
# All cancer types shown, alphabetical order, all on one PDF page.

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
samp$stratum <- factor(samp$stratum,
                       levels = c("MSS non-hypermutated",
                                  "MSS hypermutated",
                                  "MSI-H"))

sig_targets <- c("C_ID1/InsDel1a", "C_ID1/InsDel1b", "C_ID1/InsDel1c")
df_sig <- as.data.frame(t(frac[sig_targets, , drop = FALSE]),
                        check.names = FALSE)
df_sig$sample <- rownames(df_sig)
df_long <- pivot_longer(df_sig, -sample,
                        names_to = "signature", values_to = "share")
df_long <- merge(df_long, samp[, c("sample", "major", "stratum")], by = "sample")
df_long <- df_long[!is.na(df_long$major) & !is.na(df_long$stratum), ]

df_long$signature <- factor(df_long$signature, levels = sig_targets)
df_long$major <- factor(df_long$major, levels = sort(unique(df_long$major)))

# Per-violin n labels for top facet
n_labels <- df_long %>%
  dplyr::filter(signature == "C_ID1/InsDel1a") %>%
  group_by(major, stratum) %>%
  summarise(n = dplyr::n(), .groups = "drop")

p <- ggplot(df_long, aes(major, share, fill = stratum, color = stratum)) +
  geom_violin(position = position_dodge(width = 0.85),
              width = 0.85, scale = "width",
              alpha = 0.55, linewidth = 0.25,
              trim = TRUE, na.rm = TRUE) +
  stat_summary(aes(group = stratum),
               fun = median, geom = "point",
               position = position_dodge(width = 0.85),
               shape = 95, size = 4, color = "black") +
  facet_grid(signature ~ ., switch = "y") +
  scale_fill_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                               "MSS hypermutated"     = "#ff7f00",
                               "MSI-H"                = "#e31a1c"),
                    name = NULL) +
  scale_color_manual(values = c("MSS non-hypermutated" = "#1f78b4",
                                "MSS hypermutated"     = "#ff7f00",
                                "MSI-H"                = "#e31a1c"),
                     guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 9) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  labs(title = "InsDel1a / 1b / 1c per-sample share by Major Cancer Type and stratum",
       subtitle = paste0("Each sample appears in exactly one stratum. ",
                         "Hypermutator cutoff = ", HYPER_CUT, " total indels. ",
                         "Black bar = median."),
       x = NULL, y = "Signature share of sample's assigned indels")

out_pdf <- file.path(this_dir, "violins_id1_by_type_and_stratum.pdf")
pdf(out_pdf, width = 13, height = 8)
print(p)
dev.off()
message("Wrote ", out_pdf)

# Also save the per-group n table so the user can see sample counts.
n_full <- df_long %>%
  dplyr::filter(signature == "C_ID1/InsDel1a") %>%
  group_by(major, stratum) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  pivot_wider(names_from = stratum, values_from = n, values_fill = 0)
write.csv(n_full,
          file.path(this_dir, "violins_id1_sample_counts.csv"),
          row.names = FALSE)
message("Wrote violins_id1_sample_counts.csv")
print(as.data.frame(n_full))
