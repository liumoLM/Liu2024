x <- readRDS("~/github/0Liu2024/review_per_sig_txt/sigs.rds"); list2env(x, globalenv())
colnames(m89) <- sub("_converted","",colnames(m89))
cos <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))
r89 <- rownames(s89); r476 <- rownames(s476)
g <- function(M, pat, fixed=TRUE) colSums(M[grep(pat, rownames(M), fixed=fixed),,drop=FALSE])
f89 <- function(M) round(rbind(
 DelT_R14=g(M,"Del(T):R(1,4)"), DelT_R57=g(M,"Del(T):R(5,7)"), DelT_R8p=g(M,"Del(T):R(8,)"),
 InsT_R04=g(M,"Ins(T):R(0,4)"), InsT_R57=g(M,"Ins(T):R(5,7)"), InsT_R8p=g(M,"Ins(T):R(8,)"),
 DelC=g(M,"Del(C)"), InsC=g(M,"Ins(C)"), MH=g(M,":M"), Ins2p=g(M,"^Ins\\(",FALSE), Del2p=g(M,"^Del\\(",FALSE)),3)
sel <- c("InsDel1a","InsDel1b","InsDel1c","InsDel2a","InsDel2b","InsDel2c","InsDel3a","InsDel3b","InsDel5a","InsDel5b","InsDel5c","InsDel9","InsDel19a","InsDel19b","InsDel_C","InsDel13","InsDel_D","InsDel_K","InsDel_M","InsDel_N_alpha","InsDel_N_beta")
cat("89-type extracted\n"); print(t(f89(s89[,sel])))
cat("476->89\n"); print(t(f89(m89[,sel])))
cat("\ncos89 InsDel11 v 16:", cos(s89$InsDel11,s89$InsDel16), " 476:", cos(s476$InsDel11,s476$InsDel16), "\n")
cat("cos89 InsDel1a v D:", cos(s89$InsDel1a,s89$InsDel_D), " 476:", cos(s476$InsDel1a,s476$InsDel_D), "\n")
cat("cos89 InsDel13 v C:", cos(s89$InsDel13,s89$InsDel_C), " 476:", cos(s476$InsDel13,s476$InsDel_C), "\n")
# 476 single-base detail
rr <- function(sig, op, base){ v <- s476[,sig]; nm <- r476; i <- grepl(paste0("[",op,"(",base,"):R"), nm, fixed=TRUE)
  R <- sub(".*:R([^]]*)\\].*","\\1",nm[i]); tapply(v[i], factor(R, levels=unique(R)), sum)}
for (s in c("InsDel11","InsDel16")) {cat(s,"InsC by R:\n"); print(round(rr(s,"Ins","C"),3)); cat(s,"InsT by R:\n"); print(round(rr(s,"Ins","T"),3))}
for (s in c("InsDel13","InsDel_C")) {cat(s,"DelT by R:\n"); print(round(rr(s,"Del","T"),3)); cat(s,"InsT by R:\n"); print(round(rr(s,"Ins","T"),3))}
for (s in c("InsDel_N_alpha","InsDel_N_beta")) {cat(s,"DelT by R:\n"); print(round(rr(s,"Del","T"),3)); cat(s,"InsT by R:\n"); print(round(rr(s,"Ins","T"),3))}
cat("89 Del(T) followed by A (R1-4):\n"); print(round(sapply(c("InsDel13","InsDel_C"), function(s) c(xTA=sum(s89[grep("\\[Del\\(T\\):R\\(1,4\\)\\]A",r89),s]), allR14=sum(s89[grep("Del(T):R(1,4)",r89,fixed=TRUE),s]))),3))
for (s in c("InsDel_A_alpha","InsDel_A_beta")) { v <- setNames(s476[,s], r476)
 cat(s, " TCA>TA", v["T[Del(C):R1]A"], " TCT>TT", v["T[Del(C):R1]T"], " GCT>GT", v["G[Del(C):R1]T"], " GCCA>GCA", v["G[Del(C):R2]A"],
  " Del2 total", sum(v[grep("^Del2:", r476)]), " all del C", sum(v[grepl("Del(C)",r476,fixed=TRUE)]), "\n")}
for (s in c("InsDel_F","InsDel4")) { v <- setNames(s476[,s], r476); d <- v[grep("^Del",r476)]; cat(s,"\n"); print(round(sort(d,decreasing=TRUE)[1:10],3)); cat("Del2 all:", sum(v[grep("^Del2",r476)]), "\n")}
