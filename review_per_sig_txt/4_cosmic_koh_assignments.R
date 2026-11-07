x <- readRDS("~/github/0Liu2024/review_per_sig_txt/sigs.rds"); list2env(x, globalenv())
d <- "~/github/unified_indels_data/"
cos <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))
cm <- read.delim(paste0(d,"external_data/COSMIC_v3.6_ID_GRCh37.txt"), row.names=1)
# map COSMIC row names to 83 names via mSigSpectra? inspect order
cat(head(rownames(cm),8), "\n"); 
# convert "1:Del:C:0" -> "DEL:C:1:0"; "2:Del:R:0" -> DEL:repeats:2:0 ; "5:Del:M:1" -> DEL:MH:5+:1
conv <- function(r){ p <- strsplit(r,":")[[1]]; op <- toupper(p[2]); t <- p[3]; t <- c(C="C",T="T",R="repeats",M="MH")[t]
  len <- ifelse(p[1]=="5","5+",p[1]); n <- ifelse(p[4]=="5","5+",p[4]); paste(op,t,len,n,sep=":")}
rownames(cm) <- sapply(rownames(cm), conv); stopifnot(all(rownames(s83) %in% rownames(cm))); cm <- cm[rownames(s83),]
cat("ID9 DelT5p", cm["DEL:T:1:5+","ID9"], " ID9 >1bp:", sum(cm[grep("repeats:[2-5]|MH",rownames(cm)),"ID9"]), " C_ID9 >1bp:", sum(s83[grep("repeats:[2-5]|MH",rownames(s83)),"C_ID9"]),"\n")
cat("ID9 InsT5p", cm["INS:T:1:5+","ID9"], " C_ID9 InsT5p", s83["INS:T:1:5+","C_ID9"], "\n")
cat("cos C_ID9 v ID9", cos(s83$C_ID9, cm$ID9), "; C_ID4 v ID4", cos(s83$C_ID4, cm$ID4), "; ID_F v ID4", cos(s83$ID_F, cm$ID4), "\n")
k <- read.delim(paste0(d,"external_data/Koh_signatures.tsv"), row.names=1); k <- k[rownames(s89),]
cat("cos InsDel9 vs Koh: "); print(round(sort(sapply(k, cos, b=s89$InsDel9), decreasing=TRUE)[1:3],3))
g <- function(v,p) sum(v[grepl(p, rownames(s89), fixed=TRUE)])
for (nm in c("InD9b")) cat(nm, "InsT R57+R8p:", g(setNames(k[[nm]],rownames(s89)),"Ins(T):R(5,7)")+g(setNames(k[[nm]],rownames(s89)),"Ins(T):R(8,)"), "\n")
v9 <- setNames(s89$InsDel9, rownames(s89)); cat("InsDel9 89 InsT R5+:", g(v9,"Ins(T):R(5,7)")+g(v9,"Ins(T):R(8,)"), "\n")
v9 <- setNames(m89 <- NULL, NULL)
# ID_N NNLS
A <- as.matrix(s83[,c("C_ID1","C_ID2")]); f <- nnls::nnls(A, s83$ID_N); cat("ID_N ~ C_ID1+C_ID2 coef", round(f$x,3), " cos", cos(A%*%f$x, s83$ID_N), "\n")
# spectra exemplars
sp89 <- read.delim(paste0(d,"liu_et_al_89_spectra.tsv"), row.names=1, check.names=FALSE)
for (p in list(c("InsDel15","WES::TCGA-EW-A2FV-01A-11D-A17D-09"), c("InsDel16","WES::TCGA-13-0889-01A-01W-0420-08"))) {
  sp <- sp89[[p[2]]]; cat(p[1], "max abs diff 89 sig vs normalized spectrum:", max(abs(s89[[p[1]]] - sp/sum(sp))), " n indels", sum(sp), "\n")}
# assignments: ID_J and C_ID7
a83 <- read.delim(paste0(d,"liu_et_al_83_assignment.tsv"), row.names=1, check.names=FALSE)
si <- read.delim(paste0(d,"sample_info.tsv")); msi <- si$Patient[si$MSIseq_MSI.H %in% TRUE]
j <- a83["ID_J",] > 0; c7 <- a83["C_ID7",] > 0; ism <- colnames(a83) %in% msi
cat("ID_J>0:", sum(j), " of which MSI-H:", sum(j & ism), " with C_ID7:", sum(j & c7), " MSI-H & C_ID7:", sum(j&c7&ism), "\n")
cat("MSI-H tumors:", sum(ism), " MSI-H with C_ID7:", sum(c7&ism), " MSI-H with ID_J:", sum(j&ism), "\n")
a89 <- read.delim(paste0(d,"liu_et_al_89_assignment.tsv"), row.names=1, check.names=FALSE)
for (s in c("InsDel_D","InsDel1a","InsDel_J","InsDel7")) { y <- a89[s,]>0; cat(s, " n>0", sum(y), " MSI-H", sum(y & colnames(a89) %in% msi), "\n")}
