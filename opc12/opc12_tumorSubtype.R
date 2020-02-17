# Subtype opc12
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)
library(umap)
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
# p_opc12[,1:3] <- matrix(p.adjust(p = as.vector(as.matrix(p_opc12)),method = "fdr"),nrow = 96,ncol = 3)
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
         annotation_colors = phm_annotation_colors,show_rownames = F,cluster_rows = F,cluster_cols = F)
dev.off()

opc12_tpm_gsc <- opc12_tpm
opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)] <- log2(x = opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)]/100 + 1)

pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
sum(pc$Gene %in% opc12_tpm_gsc$symbol)
opc12_tpm_gsc <- opc12_tpm_gsc[match(pc$Gene,opc12_tpm_gsc$symbol),]

pc1_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
pc2_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))

pc1_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
pc2_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))

pc1_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
pc2_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))

pc1_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
pc2_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))


png("./plots/pc_opc12.png",width = 1200,height = 1000,res = 250)
par(mar=c(4.1,4.1,1,0.5))
# plot(pc1_opc12,pc2_opc12,xlim = c(-60,60),ylim = c(0,120),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 16)
# points(pc1_opc12_control,pc2_opc12_control)
plot(pc1_opc12F_control,pc2_opc12F_control,xlim = c(-60,60),ylim = c(0,140),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 1,lwd = 1.5,col = "#de2d26")
points(pc1_opc12F_sample,pc2_opc12F_sample,pch = 16,col = "#de2d26")
points(pc1_opc12M_control,pc2_opc12M_control,pch = 1,lwd = 1.5,col = "#3182bd")
points(pc1_opc12M_sample,pc2_opc12M_sample,pch = 16,col = "#3182bd")
legend("topright",legend = c("split-pool","ten-cell","female","male"),pch=c(1,16,15,15),col = c("#000000","#000000","#de2d26","#3182bd"))
dev.off()

# plot(c(pc1_opc12F_sample,pc1_opc12M_sample),c(pc2_opc12F_sample,pc2_opc12M_sample),xlim = c(-60,60),ylim = c(0,140),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 1,lwd = 1.5,col = "#de2d26")

png("./plots/hist_pc1_opc12.png",width = 1000,height = 1000,res = 250)
par(mar=c(4.1,4.1,1,0.5))
hist(c(pc1_opc12F_sample,pc1_opc12M_sample),breaks = seq(-40,40,5),ylim = c(0,12),col = "#666666",
     xlab = "PC1",ylab = "Number of observations",las = 1,
     main = "opc12 ten-cell samples")
dev.off()

# Full PC loadings - log10
{
  opc12_tpm_gsc <- opc12_tpm
  opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)] <- log10(x = opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)]/100 + 1)
  
  pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
  sum(pc$Gene %in% opc12_tpm_gsc$symbol)
  opc12_tpm_gsc <- opc12_tpm_gsc[match(pc$Gene,opc12_tpm_gsc$symbol),]
  
  pc1_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  pc2_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  
  pc1_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  pc2_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  
  pc1_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  pc2_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  
  pc1_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  pc2_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  
  
  png("./plots/pc_opc12_log10.png",width = 1200,height = 1000,res = 250)
  par(mar=c(4.1,4.1,1,0.5))
  # plot(pc1_opc12,pc2_opc12,xlim = c(-60,60),ylim = c(0,120),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 16)
  # points(pc1_opc12_control,pc2_opc12_control)
  plot(pc1_opc12F_control,pc2_opc12F_control,xlim = c(-30,30),ylim = c(0,50),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq - Log10",pch = 1,lwd = 1.5,col = "#de2d26")
  points(pc1_opc12F_sample,pc2_opc12F_sample,pch = 16,col = "#de2d26")
  points(pc1_opc12M_control,pc2_opc12M_control,pch = 1,lwd = 1.5,col = "#3182bd")
  points(pc1_opc12M_sample,pc2_opc12M_sample,pch = 16,col = "#3182bd")
  legend("topright",legend = c("split-pool","ten-cell","female","male"),pch=c(1,16,15,15),col = c("#000000","#000000","#de2d26","#3182bd"))
  dev.off()
  
  png("./plots/hist_pc1_opc12_log10.png",width = 1000,height = 1000,res = 250)
  par(mar=c(4.1,4.1,1,0.5))
  hist(c(pc1_opc12F_sample,pc1_opc12M_sample),breaks = seq(-15,15,2),ylim = c(0,12),col = "#666666",
       xlab = "PC1",ylab = "Number of observations",las = 1,
       main = "12dpi ten-cell samples - Log10")
  dev.off()
}

# UMAP
{
  opc12_umap <- umap(d = t(as.matrix(log2(opc12_tpm[,10:ncol(opc12_tpm)]+1))))
  plot(opc12_umap$layout)
  points(opc12_umap$layout[opc12.info$type == "pooled" & opc12.info$sex == "male",],col = "#3182bd",pch = 1,lwd = 1.5)
  points(opc12_umap$layout[opc12.info$type == "pooled" & opc12.info$sex == "female",],col = "#de2d26",pch = 1,lwd = 1.5)
  points(opc12_umap$layout[opc12.info$type == "ten-cell" & opc12.info$sex == "male",],col = "#3182bd",pch = 16)
  points(opc12_umap$layout[opc12.info$type == "ten-cell" & opc12.info$sex == "female",],col = "#de2d26",pch = 16)
}
