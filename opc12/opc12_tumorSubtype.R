# Subtype opc12

library(pheatmap)
library(RColorBrewer)

source("./functions/tumorSubtyping.R")
source("./functions/normalizeTPM.R")

f.opc12 <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/rsem_opc12.csv"
f.opc12.info <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/info_opc12.csv"

opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)

opc12_tpm <- normalizeTPM(rsem = opc12,index_counts = 10:ncol(opc12))
opc12_tpm_log <- opc12_tpm
row.names(opc12_tpm_log) <- opc12_tpm_log$symbol
opc12_tpm_log <- log2(opc12_tpm_log[,10:ncol(opc12_tpm_log)] + 1)

opc12_tpm_log_subtype <- opc12_tpm_log[conversion$MGI.symbol,]
opc12_tpm_log_subtype <- opc12_tpm_log_subtype[complete.cases(opc12_tpm_log_subtype),]

opc12_tpm_log_subtype_pooled <- opc12_tpm_log_subtype[,opc12.info$type == "pooled"]

mouseTumorGenes <- unlist(lapply(tumorSubtype$GeneSymbol,function(x) conversion$MGI.symbol[which(conversion$HGNC.symbol %in% x)]))


pheatmap(mat = as.matrix(opc12_tpm_log_subtype_order),
         color = rev(brewer.pal(11,"RdBu")),
         cluster_rows = F,
         scale = "none")

subtypeHM <- tumorSubtype[,2:4]
row.names(subtypeHM) <- tumorSubtype$GeneSymbol
pheatmap(mat = as.matrix(subtypeHM),
         color = rev(brewer.pal(11,"RdBu")),
         cluster_rows = F,cluster_cols = F)
plot(ecdf(sort(rank(opc12_tpm_log_subtype_order[1,],ties.method = "min")/ncol(opc12_tpm_log_subtype_order))))

# apply(opc12_tpm_log_subtype_order
opc12_tpm_10cell <- opc12_tpm[,9+which(opc12.info$type == "ten-cell")]
row.names(opc12_tpm_10cell) <- opc12_tpm$symbol
opc12_tpm_10cell_subset <- opc12_tpm_10cell[mouseTumorGenes,]
opc12_tpm_10cell_subset <- opc12_tpm_10cell_subset[complete.cases(opc12_tpm_10cell_subset),]
row.names(opc12_tpm_10cell_subset) <- conversion$HGNC.symbol[conversion$MGI.symbol %in% row.names(opc12_tpm_10cell_subset)]


opc12_tpm_subtype <- opc12_tpm[,c(3,1,10:ncol(opc12_tpm))]
names(opc12_tpm_subtype)[1:2] <- c("NAME","Description")
a <- opc12_tpm_subtype[opc12_tpm_subtype$NAME %in% mouseTumorGenes,]

a$NAME %in% mouseTumorGenes
conversion$HGNC.symbol[conversion$MGI.symbol %in% a$NAME]
a$NAME

!duplicated(unlist(lapply(a$NAME,function(x) conversion$HGNC.symbol[which(conversion$MGI.symbol %in% x)])))
