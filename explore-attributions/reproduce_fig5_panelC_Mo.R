#!/usr/bin/env Rscript
# Tight reproduction of Mo's Figure 5 panel C.
#
# A signature is shown only if some (signature, group) cell satisfies:
#   OR >= 5  AND  FDR < 0.001  AND  >= MIN_N samples in the group
#   actually carry the signature (exposure > 0).
# With MIN_N = 25, this picks exactly the 16 signatures Mo shows.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

setwd("/home/steve/github/0Liu2024")

OR_MIN  <- 5
FDR_MAX <- 1e-3
N_MIN   <- 25

meta <- read_excel("Sup Tables/Table S17 metadata of 6975 samples.xlsx")
exp  <- read.delim(
  "Sup Tables/Table S10 83-type and 89-type signature assignment.tsv",
  check.names = FALSE
)
sig_names <- sub(".*/", "", exp$Signature)
M <- as.matrix(exp[, -1]); rownames(M) <- sig_names
M <- M[, meta$Patient]

CANCERS <- c("Bladder", "Colon", "Esophagus", "Kidney", "Liver",
             "Lung", "Ovary", "Prostate", "Skin")
GROUPS  <- c(CANCERS, "MSI-H")

membership <- sapply(GROUPS, function(g) {
  if (g == "MSI-H") meta$MSI_status == "MSI"
  else              meta$`Major Cancer Type` == g
})

fisher_one <- function(present, in_group) {
  tab <- table(factor(present,  levels = c(FALSE, TRUE)),
               factor(in_group, levels = c(FALSE, TRUE)))
  ft <- fisher.test(tab)
  data.frame(OR = unname(ft$estimate), p = ft$p.value,
             n_in_grp_present = sum(present & in_group))
}

results <- expand.grid(signature = rownames(M), group = GROUPS,
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(
    {
      present  <- M[signature, ] > 0
      in_group <- membership[, group]
      r <- fisher_one(present, in_group)
      tibble(OR = r$OR, p = r$p, n_in_grp_present = r$n_in_grp_present)
    }
  ) |>
  ungroup() |>
  mutate(FDR = p.adjust(p, method = "BH"))

sig_keep <- results |>
  dplyr::filter(OR >= OR_MIN, FDR < FDR_MAX, n_in_grp_present >= N_MIN) |>
  pull(signature) |> unique()
sig_keep <- rownames(M)[rownames(M) %in% sig_keep]
cat("Selected signatures (", length(sig_keep), "):\n", sep = "")
print(sig_keep)

plot_df <- results |>
  dplyr::filter(signature %in% sig_keep, OR > 1, FDR < 0.05) |>
  mutate(
    signature  = factor(signature, levels = rev(sig_keep)),
    group      = factor(group,     levels = GROUPS),
    log10_OR   = pmin(log10(pmin(OR, 1e6)), 4),
    neglog_FDR = pmin(-log10(FDR), 100)
  )

write.csv(results,
          "explore-attributions/fig5_panelC_Mo_fisher.csv",
          row.names = FALSE)

p <- ggplot(plot_df, aes(x = group, y = signature)) +
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
         "Figure 5C reproduction (Mo-style: OR>=%g, FDR<%g, n>=%d)",
         OR_MIN, FDR_MAX, N_MIN),
       subtitle = "Non-mutually-exclusive groups")

ggsave("explore-attributions/fig5_panelC_Mo.pdf",
       p, width = 7, height = 5)
ggsave("explore-attributions/fig5_panelC_Mo.png",
       p, width = 7, height = 5, dpi = 150)

cat("Wrote explore-attributions/fig5_panelC_Mo.{pdf,png,csv}\n")
