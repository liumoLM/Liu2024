x <- readRDS("~/github/0Liu2024/review_per_sig_txt/sigs.rds"); list2env(x, globalenv())
colnames(m83) <- sub("_converted","",colnames(m83)); colnames(m89) <- sub("_converted","",colnames(m89))
f83 <- function(v){ r <- names(v) <- rownames(s83); c(
 DelT0=v["DEL:T:1:0"], DelT1=v["DEL:T:1:1"], DelT3=v["DEL:T:1:3"], DelT4=v["DEL:T:1:4"], DelT5p=v["DEL:T:1:5+"],
 InsT0=v["INS:T:1:0"], InsT3=v["INS:T:1:3"],InsT4=v["INS:T:1:4"], InsT5p=v["INS:T:1:5+"],
 DelC=sum(v[grep("DEL:C:1",r)]), InsC=sum(v[grep("INS:C:1",r)]),
 MH=sum(v[grep("MH",r)]), DelRep=sum(v[grep("DEL:repeats",r)]), InsRep=sum(v[grep("INS:repeats",r)]))}
tab <- function(M, nm) { out <- t(sapply(nm, function(n) f83(M[,n]))); colnames(out) <- sub("\\..*","",colnames(out)); round(out,3)}
ct <- read.delim("~/github/unified_indels_data/connection_table.tsv")
print(tab(s83, colnames(s83)))
cat("\n--- collapsed 476->83\n"); print(tab(m83, colnames(m83)))
