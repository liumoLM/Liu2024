#!/usr/bin/env Rscript
# Reproduce Figure 5 panel C, adding an "MSS-hyper" column for MSS samples
# with > 5,000 total indels. Samples are assigned to groups mutually
# exclusively in this order:
#   1. MSI-H        (MSI_status == "MSI")
#   2. MSS-hyper    (MSS_status == "MSS" and total indels > 5000)
#   3. one of the cancer types in CANCERS (MSS, <= 5000 indels)
# Samples not matching any group are excluded from the test.
#
# For each (89-type signature, group), Fisher's exact test compares the
# proportion of samples in the group with exposure > 0 against the
# proportion in all remaining samples (the full 6975-sample universe).
# Samples not assigned to any displayed group still contribute to the
# "remaining" pool, the way the original panel C does.
# Dot size = log10(odds ratio), dot color = -log10(FDR), as in the paper.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

setwd("/home/steve/github/0Liu2024")

meta <- read_excel("Sup Tables/Table S17 metadata of 6975 samples.xlsx")
exp  <- read.delim(
  "Sup Tables/Table S10 83-type and 89-type signature assignment.tsv",
  check.names = FALSE
)

# Strip 83-type prefix, keep only the InsDel* name on the right of "/".
sig_names <- sub(".*/", "", exp$Signature)
M <- as.matrix(exp[, -1])
rownames(M) <- sig_names

# Sanity: samples line up.
stopifnot(all(meta$Patient %in% colnames(M)))
M <- M[, meta$Patient]

total_indels <- colSums(M)

CANCERS <- c("Bladder", "Colon", "Esophagus", "Kidney", "Liver",
             "Lung", "Ovary", "Prostate", "Skin")

group <- rep(NA_character_, nrow(meta))
group[meta$MSI_status == "MSI"] <- "MSI-H"
mss_hyper <- meta$MSI_status == "MSS" & total_indels > 5000 & is.na(group)
group[mss_hyper] <- "MSS-hyper"
mss_norm <- meta$MSI_status == "MSS" & total_indels <= 5000 & is.na(group)
in_cancer <- mss_norm & meta$`Major Cancer Type` %in% CANCERS
group[in_cancer] <- meta$`Major Cancer Type`[in_cancer]

cat("Group sizes:\n")
print(table(group, useNA = "ifany"))

# Keep all samples. A sample with NA group is in "remaining" for every
# Fisher test.
M_k     <- M
group_k <- ifelse(is.na(group), "_other_", group)

GROUPS <- c(CANCERS, "MSI-H", "MSS-hyper")

fisher_one <- function(present, in_group) {
  tab <- table(factor(present,  levels = c(FALSE, TRUE)),
               factor(in_group, levels = c(FALSE, TRUE)))
  ft <- fisher.test(tab)
  data.frame(OR = unname(ft$estimate), p = ft$p.value)
}

results <- expand.grid(signature = rownames(M_k), group = GROUPS,
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(
    {
      present  <- M_k[signature, ] > 0
      in_group <- group_k == group
      r <- fisher_one(present, in_group)
      tibble(OR = r$OR, p = r$p)
    }
  ) |>
  ungroup() |>
  mutate(FDR = p.adjust(p, method = "BH"))

# Keep signatures with at least one significant, enriched cell.
sig_keep <- results |>
  dplyr::filter(FDR < 0.05, OR > 1) |>
  pull(signature) |>
  unique()

cat("Signatures shown:", length(sig_keep), "\n")

# Preserve the row order from the exposure table.
sig_keep <- rownames(M_k)[rownames(M_k) %in% sig_keep]

plot_df <- results |>
  dplyr::filter(signature %in% sig_keep, OR > 1, FDR < 0.05) |>
  mutate(
    signature  = factor(signature, levels = rev(sig_keep)),
    group      = factor(group,     levels = GROUPS),
    log10_OR   = log10(pmin(OR, 1e6)),
    neglog_FDR = pmin(-log10(FDR), 100)
  )

write.csv(
  results,
  "explore-attributions/fig5_panelC_with_MSS_hyper_fisher.csv",
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x = group, y = signature)) +
  geom_point(aes(size = log10_OR, color = neglog_FDR)) +
  scale_size_continuous(name = "log10(OR)", range = c(1, 8)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "lightblue", high = "darkblue") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Figure 5C reproduction with mutually-exclusive groups (FDR < 0.05)",
       subtitle = "MSS-hyper = MSS samples with > 5,000 total indels")

ggsave("explore-attributions/fig5_panelC_with_MSS_hyper.pdf",
       p, width = 7, height = 6)
ggsave("explore-attributions/fig5_panelC_with_MSS_hyper.png",
       p, width = 7, height = 6, dpi = 150)

cat("Wrote:\n",
    " explore-attributions/fig5_panelC_with_MSS_hyper_fisher.csv\n",
    " explore-attributions/fig5_panelC_with_MSS_hyper.pdf\n",
    " explore-attributions/fig5_panelC_with_MSS_hyper.png\n", sep = "")
