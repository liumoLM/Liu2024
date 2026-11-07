x <- readRDS("~/github/0Liu2024/review_per_sig_txt/sigs.rds"); list2env(x, globalenv())
d <- "~/github/unified_indels_data/raw_signature_from_all_hdp_SP_runs/"
cos <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))
for (f in c("SigProfiler.COSMIC83","mSigHdp.COSMIC83","SigProfiler.Koh89","mSigHdp.Koh89")) {
  r <- read.delim(paste0(d,f,".signatures.all.runs.txt"), check.names=FALSE)
  S <- if (grepl("83",f)) s83 else s89
  # sanity: max cos to every one of our sigs
  mx <- sapply(S, function(s) max(sapply(r, cos, b=s)))
  cat(f, " median over our sigs of best cos:", round(median(mx),3), "\n")
  q <- if (grepl("83",f)) "ID_N" else c("InsDel_N_alpha","InsDel_N_beta")
  for (s in q) { cs <- sapply(r, cos, b=S[[s]]); cat("  ", s, " best", round(max(cs),3), " n>0.9:", sum(cs>0.9), " of", length(cs), "\n")}
}
