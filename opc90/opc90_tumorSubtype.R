# Subtype opc90
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)

# source("./functions/tumorSubtyping.R")
source("./functions/collapseIsoforms.R")
source("./functions/normalizeTPM.R")

load("./temp/speciesConversion.RData")

f.opc90 <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/rsem_opc90.csv"
f.opc90.info <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

conversion <- convertMouseToHuman(opc90$symbol,human,mouse)

opc90$symbol <- conversion$HGNC.symbol[match(opc90$symbol,conversion$MGI.symbol)]
opc90 <- opc90[complete.cases(opc90),]

opc90 <- collapseIsoforms(opc90,10:ncol(opc90))

opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))

opc90_tpm_subtype <- opc90_tpm[,c(3,1,10:ncol(opc90_tpm))]
names(opc90_tpm_subtype)[1:2] <- c("NAME","Description")

write.table(x = "#1.2",file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(x = paste(nrow(opc90_tpm_subtype),ncol(opc90_tpm_subtype)-2,sep = "\t"),file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",append = T,row.names = F,col.names = F)
write.table(x = opc90_tpm_subtype,file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",row.names = F,append = T)

source("./ssgsea.GBM.classification/R/msig.library.12.R")
source("./ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R")

runSsGSEAwithPermutation(profile_data_file = "./ssgsea.GBM.classification/opc90.gct",number_perms = 100)

p_opc90 <- read.table("./ssgsea.GBM.classification/p_result_opc90.gct.txt",header = T, row.names = 1)
p_opc90 <- p_opc90[,4:6]
# p_opc90[,1:3] <- matrix(p.adjust(p = as.vector(as.matrix(p_opc90)),method = "fdr"),nrow = 96,ncol = 3)
p_opc90 <- p_opc90 < 0.05

phm_annotation <- data.frame(row.names = c(paste0("control.",sprintf("%02d",1:40)),paste0("sample.",sprintf("%02d",1:56))),
                             sampleType = opc90.info$type,
                             sex = opc90.info$sex,
                             mouseID = opc90.info$mouse)
phm_annotation_colors <- list(sampleType = c(pooled = "#e41a1c",`ten-cell` = "#377eb8"),
                              sex = c(female = "#006d2c",male = "#54278f"),
                              mouseID = c(F7460 = "#99d8c9",F6340 = "#2ca25f",
                                          M8170_1 = "#bcbddc",M8170_2 = "#756bb1"))
pdf(file = "~/subtype_opc90.pdf",width = 4, height = 8)
pheatmap(1-p_opc90,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
         annotation_colors = phm_annotation_colors,show_rownames = F,cluster_rows = F,cluster_cols = F)
dev.off()






