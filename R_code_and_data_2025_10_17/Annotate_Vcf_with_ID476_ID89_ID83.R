library(ICAMS)
library(indelsig.tools.lib)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load the data
# the exmaple data has five samples: SP1003, CPCT02020389T, SP1677, SP111066, SP54363
example.vcfs <- readRDS("./example.indel.vcfs.rds")

temp <- example.vcfs %>% filter(Sample=="SP1003")

temp$mutID <- paste(temp[,2],temp[,3],temp[,4],temp[,5],sep="-")
## indel476, indel89 annd indel83 (ICAMS) need the same input: CHROM, POS, REF, ALT, Sample
## annotate VCF with indel476
mutations_indel476 <- indelsig.tools.lib::indel_classifier_full(temp, genome.v = "hg19")
## annotate VCF with indel89

mutations_indel89 <- indelsig.tools.lib::indel_classifier89(temp, genome.v = "hg19")

temp.ICAMS <- temp
colnames(temp.ICAMS)[2:5] <- c("CHROM", "POS", "REF", "ALT")
mutations_indel83 <- ICAMS::VCFsToIDCatalogs(list(temp.ICAMS), ref.genome = "hg19", return.annotated.vcfs = TRUE)
mutations_indel83 <- mutations_indel83$annotated.vcfs[[1]]

# Create the full catalog
## this gives ID476 catalog
catalog <- indelsig.tools.lib::gen_fullcatalogue(mutations_indel476, sample_col = 1)

# Combine results
temp <- cbind(temp,
              "indel476.class" = mutations_indel476$type_4[match(temp$mutID, mutations_indel476$mutID)], ##type_4 column shows the classification
              "indel89.class" = mutations_indel89$type_4[match(temp$mutID, mutations_indel89$mutID)],##type_4 column shows the classification
              "indel83.class" = mutations_indel83$ID.class[match(temp$mutID, mutations_indel83$mutID)],
              "seq.context" = mutations_indel83$seq.context[match(temp$mutID, mutations_indel83$mutID)])

