# Subtype opc12
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)

# source("./functions/tumorSubtyping.R")
source("./functions/collapseIsoforms.R")
source("./functions/normalizeTPM.R")

load("./temp/speciesConversion.RData")

f.opc12 <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/rsem_opc12.csv"
f.opc12.info <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/info_opc12.csv"

opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)

conversion <- convertMouseToHuman(opc12$symbol,human,mouse)

opc12$symbol <- conversion$HGNC.symbol[match(opc12$symbol,conversion$MGI.symbol)]
opc12 <- opc12[complete.cases(opc12),]

opc12 <- collapseIsoforms(opc12,10:ncol(opc12))

opc12_tpm <- normalizeTPM(rsem = opc12,index_counts = 10:ncol(opc12))

opc12_tpm_subtype <- opc12_tpm[,c(3,1,10:ncol(opc12_tpm))]
names(opc12_tpm_subtype)[1:2] <- c("NAME","Description")

write.table(x = "#1.2",file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(x = paste(nrow(opc12_tpm_subtype),ncol(opc12_tpm_subtype)-2,sep = "\t"),file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",append = T,row.names = F,col.names = F)
write.table(x = opc12_tpm_subtype,file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",row.names = F,append = T)

source("./ssgsea.GBM.classification/R/msig.library.12.R")
source("./ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R")

runSsGSEAwithPermutation(profile_data_file = "./ssgsea.GBM.classification/opc12.gct",number_perms = 100)

p_opc12 <- read.table("./ssgsea.GBM.classification/p_result_opc12.gct.txt",header = T, row.names = 1)
p_opc12 <- p_opc12[,4:6]
p_opc12[,1:3] <- matrix(p.adjust(p = as.vector(as.matrix(p_opc12)),method = "fdr"),nrow = 96,ncol = 3)
p_opc12 <- p_opc12 < 0.05

phm_annotation <- data.frame(row.names = paste0("opc12.",sprintf("%02d",1:96)),
                             sampleType = opc12.info$type,
                             sex = opc12.info$sex,
                             mouseID = opc12.info$mouse)
phm_annotation_colors <- list(sampleType = c(pooled = "#e41a1c",`ten-cell` = "#377eb8"),
                              sex = c(female = "#006d2c",male = "#54278f"),
                              mouseID = c(F8519 = "#99d8c9",F8520 = "#2ca25f",
                                          M8516 = "#bcbddc",M8518 = "#756bb1"))
pdf(file = "~/subtype_opc12.pdf",width = 4, height = 8)
pheatmap(1-p_opc12,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
         annotation_colors = phm_annotation_colors,show_rownames = F)
dev.off()






