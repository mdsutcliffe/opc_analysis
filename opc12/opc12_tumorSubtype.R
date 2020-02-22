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
# 
# opc12_tpm_subtype <- opc12_tpm[,c(3,1,10:ncol(opc12_tpm))]
# names(opc12_tpm_subtype)[1:2] <- c("NAME","Description")
# 
# write.table(x = "#1.2",file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",row.names = F,col.names = F)
# write.table(x = paste(nrow(opc12_tpm_subtype),ncol(opc12_tpm_subtype)-2,sep = "\t"),file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",append = T,row.names = F,col.names = F)
# write.table(x = opc12_tpm_subtype,file = "./ssgsea.GBM.classification/opc12.gct",quote = F,sep = "\t",row.names = F,append = T)
# 
# source("./ssgsea.GBM.classification/R/msig.library.12.R")
# source("./ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R")
# 
# runSsGSEAwithPermutation(profile_data_file = "./ssgsea.GBM.classification/opc12.gct",number_perms = 100)
# 
# p_opc12 <- read.table("./ssgsea.GBM.classification/p_result_opc12.gct.txt",header = T, row.names = 1)
# p_opc12 <- p_opc12[,4:6]
# # p_opc12[,1:3] <- matrix(p.adjust(p = as.vector(as.matrix(p_opc12)),method = "fdr"),nrow = 96,ncol = 3)
# p_opc12 <- p_opc12 < 0.05
# 
# phm_annotation <- data.frame(row.names = paste0("opc12.",sprintf("%02d",1:96)),
#                              sampleType = opc12.info$type,
#                              sex = opc12.info$sex,
#                              mouseID = opc12.info$mouse)
# phm_annotation_colors <- list(sampleType = c(pooled = "#e41a1c",`ten-cell` = "#377eb8"),
#                               sex = c(female = "#006d2c",male = "#54278f"),
#                               mouseID = c(F8519 = "#99d8c9",F8520 = "#2ca25f",
#                                           M8516 = "#bcbddc",M8518 = "#756bb1"))
# pdf(file = "~/subtype_opc12.pdf",width = 4, height = 8)
# pheatmap(1-p_opc12,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
#          annotation_colors = phm_annotation_colors,show_rownames = F,cluster_rows = F,cluster_cols = F)
# dev.off()

opc12_tpm_gsc <- opc12_tpm
opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)] <- log2(x = opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)]/100 + 1)

pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
pc$PC1 <- pc$PC1 / sqrt(sum(pc$PC1^2))
pc$PC2 <- pc$PC1 / sqrt(sum(pc$PC2^2))
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
  # opc12_tpm_gsc <- opc12_tpm
  # opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)] <- log10(x = opc12_tpm_gsc[,10:ncol(opc12_tpm_gsc)]/100 + 1)
  # 
  # pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
  # sum(pc$Gene %in% opc12_tpm_gsc$symbol)
  # opc12_tpm_gsc <- opc12_tpm_gsc[match(pc$Gene,opc12_tpm_gsc$symbol),]
  # 
  # pc1_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  # pc2_opc12F_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  # 
  # pc1_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  # pc2_opc12M_sample <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  # 
  # pc1_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  # pc2_opc12F_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "female"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  # 
  # pc1_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC1 * opc12_tpm_gsc[,x],na.rm = T))
  # pc2_opc12M_control <- sapply(9+which(opc12.info$type == "pooled" & opc12.info$sex == "male"),function(x) sum(pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])),na.rm = T))
  # 
  # 
  # png("./plots/pc_opc12_log10.png",width = 1200,height = 1000,res = 250)
  # par(mar=c(4.1,4.1,1,0.5))
  # # plot(pc1_opc12,pc2_opc12,xlim = c(-60,60),ylim = c(0,120),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 16)
  # # points(pc1_opc12_control,pc2_opc12_control)
  # plot(pc1_opc12F_control,pc2_opc12F_control,xlim = c(-30,30),ylim = c(0,50),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq - Log10",pch = 1,lwd = 1.5,col = "#de2d26")
  # points(pc1_opc12F_sample,pc2_opc12F_sample,pch = 16,col = "#de2d26")
  # points(pc1_opc12M_control,pc2_opc12M_control,pch = 1,lwd = 1.5,col = "#3182bd")
  # points(pc1_opc12M_sample,pc2_opc12M_sample,pch = 16,col = "#3182bd")
  # legend("topright",legend = c("split-pool","ten-cell","female","male"),pch=c(1,16,15,15),col = c("#000000","#000000","#de2d26","#3182bd"))
  # dev.off()
  # 
  # png("./plots/hist_pc1_opc12_log10.png",width = 1000,height = 1000,res = 250)
  # par(mar=c(4.1,4.1,1,0.5))
  # hist(c(pc1_opc12F_sample,pc1_opc12M_sample),breaks = seq(-15,15,2),ylim = c(0,12),col = "#666666",
  #      xlab = "PC1",ylab = "Number of observations",las = 1,
  #      main = "12dpi ten-cell samples - Log10")
  # dev.off()
}

# UMAP
{

  custom.settings <- umap.defaults
  custom.settings$n_epochs <- 1000
  custom.settings$n_neighbors <- 7
  opc12_umap <- umap(d = t(as.matrix(log2(opc12_tpm[,9+which(opc12.info$type == "ten-cell")]+1))),config = custom.settings)
  
  png(filename = "opc12_umap_wholeTranscriptome.png",width = 1000,height = 1000,res = 300)
  par(mar=c(4.1,4.1,1,0.5))
  plot(opc12_umap$layout,xlab = "UMAP 1",ylab = "UMAP 2")
  points(opc12_umap$layout[opc12.info$sex[opc12.info$type == "ten-cell"] == "male",],col = "#3182bd",pch = 16)
  points(opc12_umap$layout[opc12.info$sex[opc12.info$type == "ten-cell"] == "female",],col = "#de2d26",pch = 16)
  dev.off()
}

# Clustergram of PC transcripts
{
  


opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)

opc12 <- collapseIsoforms(opc12,10:ncol(opc12))

opc12_tpm <- normalizeTPM(rsem = opc12,index_counts = 10:ncol(opc12))

pc_mouse <- convertHumanToMouse(pc$Gene,human,mouse)

x <- pc_mouse$MGI.symbol[match(pc_mouse$HGNC.symbol,pc$Gene)]
x <- pc_mouse$MGI.symbol[match(pc$Gene,pc_mouse$HGNC.symbol)]
x <- x[!is.na(x)]
x <- unique(x)

opc12_pc <- opc12_tpm[match(x,opc12_tpm$symbol),]
opc12_pc <- opc12_pc[complete.cases(opc12_pc),]

row.names(opc12_pc) <- opc12_pc$symbol
opc12_pc <- opc12_pc[,10:ncol(opc12_pc)]

library(pheatmap)

backConvert <- convertMouseToHuman(row.names(opc12_pc),human,mouse)
bC2 <- backConvert[!duplicated(backConvert$MGI.symbol) & !duplicated(backConvert$HGNC.symbol),]
bC2$HGNC.symbol[match(row.names(opc12_pc),bC2$MGI.symbol)]
pc_bc <- pc[pc$Gene %in% bC2$HGNC.symbol[match(row.names(opc12_pc),bC2$MGI.symbol)],2]
temp <- as.character(pc[pc$Gene %in% bC2$HGNC.symbol[match(row.names(opc12_pc),bC2$MGI.symbol)],1])
tempConvert <- convertHumanToMouse(temp,human,mouse)
nrow(tempConvert[match(temp,tempConvert$HGNC.symbol),])
unique(tempConvert[match(temp,tempConvert$HGNC.symbol),2]) == row.names(opc12_pc)
pc[x,1:2]
ann <- data.frame(row.names = names(opc12_pc[,opc12.info$type == "ten-cell"]),
                  sex = opc12.info$sex[opc12.info$type == "ten-cell"])
png(filename = "./plots/opc12_pc_clustergram.png",width = 1200,height = 1200,res = 200)
pheatmap(mat = as.matrix(log2(opc12_pc[,opc12.info$type == "ten-cell"]+1)),scale = "column",cluster_rows = F,color = rev(brewer.pal(11,"RdBu")),show_rownames = F)
dev.off()

opc12_pc_subset <- opc12_pc[c(1:50,(nrow(opc12_pc)-49):nrow(opc12_pc)),]
png(filename = "./plots/opc12_pc_clustergram_subset.png",width = 1000,height = 1200,res = 150)
pheatmap(mat = as.matrix(log2(opc12_pc_subset[,opc12.info$type == "ten-cell"]+1)),scale = "column",cluster_rows = F,color = rev(brewer.pal(11,"RdBu")),show_rownames = F,gaps_row = rep(50,3),border_color = NA)
dev.off()

png(filename = "./plots/opc12_pc_clustergram_noScale.png",width = 1200,height = 1200,res = 200)
pheatmap(mat = as.matrix(log2(opc12_pc[,opc12.info$type == "ten-cell"]+1)),scale = "none",cluster_rows = F,color = brewer.pal(9,"Reds"),show_rownames = F)
dev.off()

opc12_pc_subset <- opc12_pc[c(1:50,(nrow(opc12_pc)-49):nrow(opc12_pc)),]
png(filename = "./plots/opc12_pc_clustergram_noScale_subset.png",width = 1000,height = 1200,res = 150)
pheatmap(mat = as.matrix(log2(opc12_pc_subset[,opc12.info$type == "ten-cell"]+1)),scale = "none",cluster_rows = F,color = brewer.pal(9,"Reds"),show_rownames = F,gaps_row = rep(50,3),border_color = NA)
dev.off()
}

pc1_opc12F_sample_load <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) pc$PC1 * opc12_tpm_gsc[,x])
pc2_opc12F_sample_load <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),function(x) pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])))

pc1_opc12M_sample_load <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) pc$PC1 * opc12_tpm_gsc[,x])
pc2_opc12M_sample_load <- sapply(9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male"),function(x) pc$PC2 * (opc12_tpm_gsc[,x]-(pc$PC1 * opc12_tpm_gsc[,x])))

row.names(pc1_opc12F_sample_load) <- pc$Gene
row.names(pc2_opc12F_sample_load) <- pc$Gene
row.names(pc1_opc12M_sample_load) <- pc$Gene
row.names(pc2_opc12M_sample_load) <- pc$Gene

colnames(pc1_opc12F_sample_load) <- names(opc12_tpm_gsc)[9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female")]
colnames(pc2_opc12F_sample_load) <- names(opc12_tpm_gsc)[9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "female")]
colnames(pc1_opc12M_sample_load) <- names(opc12_tpm_gsc)[9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male")]
colnames(pc2_opc12M_sample_load) <- names(opc12_tpm_gsc)[9+which(opc12.info$type == "ten-cell" & opc12.info$sex == "male")]

m <- cbind(pc1_opc12F_sample_load,pc1_opc12M_sample_load)
m <- m[complete.cases(m),]
ann <- data.frame(row.names = c(colnames(pc1_opc12F_sample_load),colnames(pc1_opc12M_sample_load)),
                  sex = c(rep("F",28),rep("M",28)),gaps_row = rep(28,2))

png(filename = "./plots/opc12_pc_clustergram_load.png",width = 1400,height = 1000,res = 150)
pheatmap(mat = t(m[nrow(m):1,]),cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),breaks = seq(-4,4,length.out = 12),show_colnames = F,main = "12 dpi - 10cRNA-seq projections - PC1",
         annotation_row = ann,cluster_rows = F,gaps_row = rep(28,2))
dev.off()
png(filename = "./plots/opc12_pc_clustergram_load_subset.png",width = 1400,height = 1000,res = 100)
pheatmap(mat = t(m[c(nrow(m):(nrow(m)-49),50:1),]),cluster_cols = F,gaps_col = rep(50,3),color = rev(brewer.pal(11,"RdBu")),main = "12 dpi - 10cRNA-seq projections - PC1",border_color = NA,breaks = seq(-4,4,length.out = 12),annotation_row = ann,cluster_rows = F,gaps_row = rep(28,2))
dev.off()
