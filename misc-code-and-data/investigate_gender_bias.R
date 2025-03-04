library(tidyr)
library(dplyr)
rm(list = ls())

signature.assignment <- readRDS("./signature.assignment.df.rds") ##signature assignment result
(signature.assignment > 0.05) %>%
t() %>% as_tibble() %>%
  `$<-`(new_sample_ID, colnames(signature.assignment)) %>%
  relocate(new_sample_ID) -> sa
# sa looks like: 
# A tibble: 6,176 × 34
# new_sample_ID     C_ID1 C_ID2 C_ID3 C_ID4 C_ID5 C_ID6 C_ID7 C_ID8 C_ID9 C_ID10 C_ID11 C_ID12 C_ID13 C_ID14 C_ID17 C_ID18 C_ID19 C_ID23 H_ID24 H_ID25 H_ID26
# <chr>             <lgl> <lgl> <lgl> <lgl> <lgl> <lgl> <lgl> <lgl> <lgl> <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl>  <lgl> 
# 1 Biliary::CPCT020… TRUE  TRUE  FALSE FALSE TRUE  TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  TRUE   FALSE  FALSE  FALSE  TRUE  
# 2 Biliary::CPCT020… FALSE TRUE  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  TRUE   FALSE  FALSE  TRUE   FALSE 
# 3 Biliary::CPCT020… TRUE  TRUE  FALSE TRUE  FALSE FALSE FALSE FALSE TRUE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  TRUE   FALSE  FALSE  FALSE  FALSE 
# 4 Biliary::CPCT020… TRUE  TRUE  FALSE FALSE TRUE  TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  TRUE   FALSE  FALSE  FALSE  FALSE 
# 5 Biliary::CPCT020… TRUE  TRUE  FALSE TRUE  FAL
stopifnot(ncol(signature.assignment) == nrow(sa))
any(duplicated(sa$new_sample_ID))

metainfo <- readRDS("./metainfo.rds") ## meta info for PCAWG and HMF
## I only used the samples with reconstruction cos.sim > 0.9, so some samples don't have assignment result
metainfo <- metainfo[metainfo$new_sample_ID %in% colnames(signature.assignment),]
stopifnot(nrow(metainfo) == nrow(sa))
nrow(filter(metainfo, is.na(Gender)))
nrow(filter(metainfo, Gender == ""))

metainfo %>%
  filter(!is.na(Gender)) %>%
  filter(Gender != "") -> m2
  stopifnot(nrow(m2)== 4357) 

m2 %>% 
  dplyr::left_join(sa) -> tt
stopifnot(nrow(m2) == nrow(tt))
stopifnot(!any(duplicated(tt$new_sample_ID)))

tt %>%
  filter(!(cancertype %in% c("Prostate","Uterus","Breast","Ovary"))) %>%
  group_by(cancertype, Gender) %>% select(-c(2, 5:12)) ->
  tt2
  
stopifnot(nrow(tt2) == 3067)

summarize(tt2, 
          across(
            .cols = 2:34,
            .fns = list(
              Tcount = ~ sum(.),
              Fcount = ~ sum(!.)
            ),
            .names = "{.col}.{.fn}")) -> s1
# s1 looks like:
# # A tibble: 36 × 68
# Groups:   cancertype [18]
# cancertype      Gender C_ID1.Tcount C_ID1.Fcount C_ID2.Tcount C_ID2.Fcount C_ID3.Tcount C_ID3.Fcount C_ID4.Tcount C_ID4.Fcount C_ID5.Tcount C_ID5.Fcount C_ID6.Tcount C_ID6.Fcount C_ID7.Tcount C_ID7.Fcount C_ID8.Tcount C_ID8.Fcount
# <chr>           <chr>         <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>        <int>
# 1 Biliary         Female           38            1           38            1            0           39           11           28            7           32            5           34            1           38            3           36
# 2 Biliary         Male             42            1           38            5            0           43           13           30           10           33            4           39            2           41            2           41
# 3 Bladder         Female           25            2           27            0            0           27            7           20            4           23            3           24            8           19            8           19
# 4 Bladder         Male             70            1           71            0            0           71           19           52           20           51           11           60           14           57           11           60

s1 %>% 
  pivot_longer(
    # pivot all columns except the grouping variables
    cols = -c(cancertype, Gender),
    # specify that we want to create two output variables (columns)
    # from each column name in the input row, e.g.
    # C_ID1.Tcount generates a value for new column sig ("C_ID1)
    # takes the value in C_ID1.Tcount and puts it column Tcount
    names_to = c("sig", ".value"),
    names_sep = "\\."
  ) %>% 
  arrange(cancertype, sig, Gender) %>%
  relocate(sig, .after = cancertype) -> s2

# s2 looks like:
# A tibble: 1,188 × 5
# Groups:   cancertype [18]
# cancertype sig    Gender Tcount Fcount
# <chr>      <chr>  <chr>   <int>  <int>
# 1 Biliary    C_ID1  Female     38      1
# 2 Biliary    C_ID1  Male       42      1
# 3 Biliary    C_ID10 Female      8     31
# 4 Biliary    C_ID10 Male       11     32
# ...

# Remove cancertype, sig groups in which the
# sig (signature) is rare or almost always present.
s2 %>%
  group_by(cancertype, sig) %>%
  filter(sum(Tcount) > 20) %>%
  filter(sum(Fcount) > 20) ->
  s3

fishtest = function(tfmatrix) {
  out1 = fisher.test(tfmatrix[[1]])
  # browser() # out1$estimate is odds ratio
  return(out1$p.value)
}

# Did not figure out how to return 2 columns 
# as output from one summarizing funtion...
oddsratio = function(tfmatrix) {
  out1 = fisher.test(tfmatrix[[1]])
  # browser() # out1$estimate is odds ratio
  return(out1$estimate)
}

make_matrix = function(Gender,tcount, fcount) {
  m1 = cbind(tcount, fcount)
  rownames(m1) = Gender
  return(m1)
}

s3 %>% 
  mutate(tfmatrix = list(make_matrix(Gender, Tcount, Fcount)),
         p.value = fishtest(tfmatrix),
         odds_ratio = oddsratio(tfmatrix)) ->
  s4

# Male and Female rows contain the same tfmatrix and p.value values.
# We want only one before computing FDRs
s5 = filter(s4, Gender == "Female")

s5$p.adj = p.adjust(s5$p.value)
possible = filter(s5, p.adj < 0.20)

possible
# A tibble: 3 × 8
# Groups:   cancertype, sig [3]
#  cancertype sig    Gender Tcount Fcount tfmatrix          p.value odds_ratio     p.adj
# <chr>      <chr>  <chr>   <int>  <int> <list>              <dbl>      <dbl>     <dbl>
# 1 Lung       C_ID3  Female     74     32 <int [2 × 2]> 0.000252        0.281  0.0321   
# 2 Other      C_ID10 Female     29     30 <int [2 × 2]> 0.000456        4.25   0.0574   
# 3 Other      C_ID19 Female     36     23 <int [2 × 2]> 0.000000106     0.0272 0.0000136
# 4 Skin       C_ID13 Female     62     34 <int [2 × 2]> 0.00145         0.383  0.181    
# 
# Speculation: Differences in signatures in Other represent difference in 
# cancer subtype prevalence between men and women.

possible[1, ]$tfmatrix[[1]]
#        tcount fcount
# Female     74     32
# Male      124     15
# In Lung, ID3 is enriched in men, probably reflects more
# lung cancer cases in men due to smoking than 
# lung cancer cases in women.

possible[4, ]$tfmatrix[[1]]
#        tcount fcount
# Female     62     34
# Male      134     28
# In Skin, C_ID13 is enriched in men. Why? Different types of UV exposure?

library(BonEV)
s5$bonev = Bon_EV(s5$p.value, 0.1)$Bon_EV_adjp

lowodds = head(arrange(s5, odds_ratio), 15)
lowodds
# Liver, H_ID25 has low p.value and strong odds ratio
filter(lowodds, cancertype == "Liver", sig == "H_ID25")$tfmatrix
#        tcount fcount
# Female     44     45
# Male      157     83
# --> Any hints on origin of H_ID25? More prevalent in men.
# 
# Kidney, H_ID19 has low p.value and strong odds ratio
filter(lowodds, cancertype == "Kidney", sig == "C_ID19")$tfmatrix
#        tcount fcount
# Female     12     52
# Male       47     92
# --> Any hints on origin of H_ID19? More prevalent in men.

junk = apply(filter(lowodds, p.value < 0.05),
             MARGIN = 1, 
             function(rr) { 
               print(paste(rr$cancertype, rr$sig, rr$odds_ratio, rr$p.value))
               print(rr$tfmatrix)
               cat("\n")}
)


highodds = head(arrange(s5, desc(odds_ratio)), 10)
highodds

junk = apply(
  filter(highodds, p.value < 0.05),
  MARGIN = 1,
  function(rr) { 
    print(paste(rr$cancertype, rr$sig, rr$odds_ratio, rr$p.value))
    print(rr$tfmatrix)
    cat("\n")}
)
