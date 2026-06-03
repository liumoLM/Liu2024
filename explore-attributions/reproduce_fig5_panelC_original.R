#!/usr/bin/env Rscript
# Reproduce Figure 5 panel C in its original form: cancer-type columns
# plus MSI-H, with groups NOT mutually exclusive. A sample is in its
# cancer-type column whether or not it is also MSI-H.
#
# For each (89-type signature, group), Fisher's exact test compares the
# proportion of samples in the group with exposure > 0 against the
# proportion in all remaining samples in the 6975-sample universe.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

setwd("/home/steve/github/0Liu2024")

meta <- read_excel("Sup Tables/Table S17 metadata of 6975 samples.xlsx")
exp  <- read.delim(
  "Sup Tables/Table S10 83-type and 89-type signature assignment.tsv",
  check.names = FALSE
)

sig_names <- sub(".*/", "", exp$Signature)
M <- as.matrix(exp[, -1])
rownames(M) <- sig_names

stopifnot(all(meta$Patient %in% colnames(M)))
M <- M[, meta$Patient]

CANCERS <- c("Bladder", "Colon", "Esophagus", "Kidney", "Liver",
             "Lung", "Ovary", "Prostate", "Skin")
GROUPS  <- c(CANCERS, "MSI-H")

# Build a sample-by-group membership matrix. A sample can belong to
# more than one group (e.g. an MSI Bladder sample is in both Bladder
# and MSI-H).
membership <- sapply(GROUPS, function(g) {
  if (g == "MSI-H") meta$MSI_status == "MSI"
  else              meta$`Major Cancer Type` == g
})
rownames(membership) <- meta$Patient

cat("Group sizes (non-mutually-exclusive):\n")
print(colSums(membership))

fisher_one <- function(present, in_group) {
  tab <- table(factor(present,  levels = c(FALSE, TRUE)),
               factor(in_group, levels = c(FALSE, TRUE)))
  ft <- fisher.test(tab)
  data.frame(OR = unname(ft$estimate), p = ft$p.value)
}

results <- expand.grid(signature = rownames(M), group = GROUPS,
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(
    {
      present  <- M[signature, ] > 0
      in_group <- membership[, group]
      r <- fisher_one(present, in_group)
      tibble(OR = r$OR, p = r$p)
    }
  ) |>
  ungroup() |>
  mutate(FDR = p.adjust(p, method = "BH"))

write.csv(
  results,
  "explore-attributions/fig5_panelC_original_fisher.csv",
  row.names = FALSE
)

make_plot <- function(cutoff) {
  sig_keep <- results |>
    dplyr::filter(FDR < cutoff, OR > 1) |>
    pull(signature) |>
    unique()
  sig_keep <- rownames(M)[rownames(M) %in% sig_keep]
  cat(sprintf("FDR < %g: %d signatures\n", cutoff, length(sig_keep)))

  plot_df <- results |>
    dplyr::filter(signature %in% sig_keep, OR > 1, FDR < cutoff) |>
    mutate(
      signature  = factor(signature, levels = rev(sig_keep)),
      group      = factor(group,     levels = GROUPS),
      log10_OR   = pmin(log10(pmin(OR, 1e6)), 4),
      neglog_FDR = pmin(-log10(FDR), 100)
    )

  ggplot(plot_df, aes(x = group, y = signature)) +
    geom_point(aes(size = log10_OR, color = neglog_FDR)) +
    scale_size_continuous(name = "log10(OR)", range = c(1, 8),
                          limits = c(0, 4)) +
    scale_color_gradient(name = "-log10(FDR)",
                         low = "lightblue", high = "darkblue",
                         limits = c(0, 5), oob = scales::squish) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank()) +
    labs(x = NULL, y = NULL,
         title = sprintf(
           "Figure 5C reproduction (non-mutually-exclusive groups, FDR < %g)",
           cutoff),
         subtitle = "Cancer types + MSI-H, as in the original figure")
}

cutoffs <- c(0.001, 0.01, 0.05, 0.10, 0.25)
for (cut in cutoffs) {
  tag <- sub("\\.", "p", format(cut, nsmall = 0, scientific = FALSE))
  p <- make_plot(cut)
  ggsave(sprintf("explore-attributions/fig5_panelC_original_FDR%s.pdf", tag),
         p, width = 7, height = 6)
  ggsave(sprintf("explore-attributions/fig5_panelC_original_FDR%s.png", tag),
         p, width = 7, height = 6, dpi = 150)
}

# Also keep the default (FDR < 0.05) at the legacy filename.
p_default <- make_plot(0.05)
ggsave("explore-attributions/fig5_panelC_original.pdf",
       p_default, width = 7, height = 6)
ggsave("explore-attributions/fig5_panelC_original.png",
       p_default, width = 7, height = 6, dpi = 150)

cat("Wrote fisher CSV + one PDF/PNG per FDR cutoff in",
    paste(cutoffs, collapse = ", "), "\n")
