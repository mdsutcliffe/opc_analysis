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
# 
# opc90_tpm_subtype <- opc90_tpm[,c(3,1,10:ncol(opc90_tpm))]
# names(opc90_tpm_subtype)[1:2] <- c("NAME","Description")
# 
# write.table(x = "#1.2",file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",row.names = F,col.names = F)
# write.table(x = paste(nrow(opc90_tpm_subtype),ncol(opc90_tpm_subtype)-2,sep = "\t"),file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",append = T,row.names = F,col.names = F)
# write.table(x = opc90_tpm_subtype,file = "./ssgsea.GBM.classification/opc90.gct",quote = F,sep = "\t",row.names = F,append = T)
# 
# source("./ssgsea.GBM.classification/R/msig.library.12.R")
# source("./ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R")
# 
# runSsGSEAwithPermutation(profile_data_file = "./ssgsea.GBM.classification/opc90.gct",number_perms = 100)
# 
# p_opc90 <- read.table("./ssgsea.GBM.classification/p_result_opc90.gct.txt",header = T, row.names = 1)
# p_opc90 <- p_opc90[,4:6]
# # p_opc90[,1:3] <- matrix(p.adjust(p = as.vector(as.matrix(p_opc90)),method = "fdr"),nrow = 96,ncol = 3)
# p_opc90 <- p_opc90 < 0.05
# 
# phm_annotation <- data.frame(row.names = c(paste0("control.",sprintf("%02d",1:40)),paste0("sample.",sprintf("%02d",1:56))),
#                              sampleType = opc90.info$type,
#                              sex = opc90.info$sex,
#                              mouseID = opc90.info$mouse)
# phm_annotation_colors <- list(sampleType = c(pooled = "#e41a1c",`ten-cell` = "#377eb8"),
#                               sex = c(female = "#006d2c",male = "#54278f"),
#                               mouseID = c(F7460 = "#99d8c9",F6340 = "#2ca25f",
#                                           M8170_1 = "#bcbddc",M8170_2 = "#756bb1"))
# pdf(file = "~/subtype_opc90.pdf",width = 4, height = 8)
# pheatmap(1-p_opc90,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
#          annotation_colors = phm_annotation_colors,show_rownames = F,cluster_rows = F,cluster_cols = F)
# dev.off()
# 

opc90_tpm_gsc <- opc90_tpm
opc90_tpm_gsc[,10:ncol(opc90_tpm_gsc)] <- log2(x = opc90_tpm_gsc[,10:ncol(opc90_tpm_gsc)]/100 + 1)

pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
sum(pc$Gene %in% opc90_tpm_gsc$symbol)
opc90_tpm_gsc <- opc90_tpm_gsc[match(pc$Gene,opc90_tpm_gsc$symbol),]

pc1_opc90 <- sapply(9+which(opc90.info$type == "ten-cell"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90 <- sapply(9+which(opc90.info$type == "ten-cell"),function(x) sum(pc$PC2 * opc90_tpm_gsc[,x],na.rm = T))

pc1_opc90_control <- sapply(9+which(opc90.info$type == "pooled"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90_control <- sapply(9+which(opc90.info$type == "pooled"),function(x) sum(pc$PC2 * opc90_tpm_gsc[,x],na.rm = T))

pc1_opc90F_sample <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90F_sample <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female"),function(x) sum(pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])),na.rm = T))

pc1_opc90M_sample <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90M_sample <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male"),function(x) sum(pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])),na.rm = T))

pc1_opc90F_control <- sapply(9+which(opc90.info$type == "pooled" & opc90.info$sex == "female"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90F_control <- sapply(9+which(opc90.info$type == "pooled" & opc90.info$sex == "female"),function(x) sum(pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])),na.rm = T))

pc1_opc90M_control <- sapply(9+which(opc90.info$type == "pooled" & opc90.info$sex == "male"),function(x) sum(pc$PC1 * opc90_tpm_gsc[,x],na.rm = T))
pc2_opc90M_control <- sapply(9+which(opc90.info$type == "pooled" & opc90.info$sex == "male"),function(x) sum(pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])),na.rm = T))

png("./plots/pc_opc90.png",width = 1200,height = 1000,res = 250)
par(mar=c(4.1,4.1,1,0.5))
# plot(pc1_opc90,pc2_opc90,xlim = c(-60,60),ylim = c(0,120),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "12 dpi 10cRNA-seq",pch = 16)
# points(pc1_opc90_control,pc2_opc90_control)
plot(pc1_opc90F_control,pc2_opc90F_control,xlim = c(-60,60),ylim = c(0,140),xlab = "PC1",ylab = "PC2",las = 1,frame = F,main = "90 dpi 10cRNA-seq",pch = 1,lwd = 1.5,col = "#de2d26")
points(pc1_opc90F_sample,pc2_opc90F_sample,pch = 16,col = "#de2d26")
points(pc1_opc90M_control,pc2_opc90M_control,pch = 1,lwd = 1.5,col = "#3182bd")
points(pc1_opc90M_sample,pc2_opc90M_sample,pch = 16,col = "#3182bd")
legend("topright",legend = c("split-pool","ten-cell","female","male"),pch=c(1,16,15,15),col = c("#000000","#000000","#de2d26","#3182bd"))
dev.off()



# Clustergram of PC transcripts

f.opc90 <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/rsem_opc90.csv"
f.opc90.info <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

opc90 <- collapseIsoforms(opc90,10:ncol(opc90))

opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))

pc_mouse <- convertHumanToMouse(pc$Gene,human,mouse)

x <- pc_mouse$MGI.symbol[match(pc_mouse$HGNC.symbol,pc$Gene)]
x <- pc_mouse$MGI.symbol[match(pc$Gene,pc_mouse$HGNC.symbol)]
x <- x[!is.na(x)]
x <- unique(x)

opc90_pc <- opc90_tpm[match(x,opc90_tpm$symbol),]
opc90_pc <- opc90_pc[complete.cases(opc90_pc),]

row.names(opc90_pc) <- opc90_pc$symbol
opc90_pc <- opc90_pc[,10:ncol(opc90_pc)]

png(filename = "./plots/opc90_pc_clustergram.png",width = 1200,height = 1200,res = 200)
pheatmap(mat = as.matrix(log2(opc90_pc[,opc90.info$type == "ten-cell"]+1)),scale = "column",cluster_rows = F,color = rev(brewer.pal(11,"RdBu")),show_rownames = F)
dev.off()

opc90_pc_subset <- opc90_pc[c(1:50,(nrow(opc90_pc)-49):nrow(opc90_pc)),]
png(filename = "./plots/opc90_pc_clustergram_subset.png",width = 1000,height = 1200,res = 150)
pheatmap(mat = as.matrix(log2(opc90_pc_subset[,opc90.info$type == "ten-cell"]+1)),scale = "column",cluster_rows = F,color = rev(brewer.pal(11,"RdBu")),show_rownames = F,gaps_row = rep(50,3),border_color = NA)
dev.off()

png(filename = "./plots/opc90_pc_clustergram_noScale.png",width = 1200,height = 1200,res = 200)
pheatmap(mat = as.matrix(log2(opc90_pc[,opc90.info$type == "ten-cell"]+1)),scale = "none",cluster_rows = F,color = brewer.pal(9,"Reds"),show_rownames = F)
dev.off()

opc90_pc_subset <- opc90_pc[c(1:50,(nrow(opc90_pc)-49):nrow(opc90_pc)),]
png(filename = "./plots/opc90_pc_clustergram_noScale_subset.png",width = 1000,height = 1200,res = 150)
pheatmap(mat = as.matrix(log2(opc90_pc_subset[,opc90.info$type == "ten-cell"]+1)),scale = "none",cluster_rows = F,color = brewer.pal(9,"Reds"),show_rownames = F,gaps_row = rep(50,3),border_color = NA)
dev.off()




pc1_opc90F_sample_load <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female"),function(x) pc$PC1 * opc90_tpm_gsc[,x])
pc2_opc90F_sample_load <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female"),function(x) pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])))

pc1_opc90M_sample_load <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male"),function(x) pc$PC1 * opc90_tpm_gsc[,x])
pc2_opc90M_sample_load <- sapply(9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male"),function(x) pc$PC2 * (opc90_tpm_gsc[,x]-(pc$PC1 * opc90_tpm_gsc[,x])))

row.names(pc1_opc90F_sample_load) <- pc$Gene
row.names(pc2_opc90F_sample_load) <- pc$Gene
row.names(pc1_opc90M_sample_load) <- pc$Gene
row.names(pc2_opc90M_sample_load) <- pc$Gene

colnames(pc1_opc90F_sample_load) <- names(opc90_tpm_gsc)[9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female")]
colnames(pc2_opc90F_sample_load) <- names(opc90_tpm_gsc)[9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "female")]
colnames(pc1_opc90M_sample_load) <- names(opc90_tpm_gsc)[9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")]
colnames(pc2_opc90M_sample_load) <- names(opc90_tpm_gsc)[9+which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")]

m <- cbind(pc1_opc90F_sample_load,pc1_opc90M_sample_load)
m <- m[complete.cases(m),]
ann <- data.frame(row.names = c(colnames(pc1_opc90F_sample_load),colnames(pc1_opc90M_sample_load)),
                  sex = c(rep("F",28),rep("M",28)),
                  PC1 = c(pc1_opc90F_sample,pc1_opc90M_sample))

png(filename = "./plots/opc90_pc_clustergram_load.png",width = 1400,height = 1000,res = 150)
pheatmap(mat = t(m[nrow(m):1,]),cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),breaks = seq(-4,4,length.out = 12),show_colnames = F,main = "90 dpi - 10cRNA-seq projections - PC1",
         annotation_row = ann,gaps_row = rep(28,2))
dev.off()
png(filename = "./plots/opc90_pc_clustergram_load_subset.png",width = 1400,height = 1000,res = 100)
pheatmap(mat = t(m[c(nrow(m):(nrow(m)-49),50:1),]),cluster_cols = F,gaps_col = rep(50,3),color = rev(brewer.pal(11,"RdBu")),main = "90 dpi - 10cRNA-seq projections - PC1",border_color = NA,breaks = seq(-4,4,length.out = 12),annotation_row = ann,cluster_rows = F,gaps_row = rep(28,2))
dev.off()

mt <- t(m)
mt <- scale(mt)
mt <- mt[,c(order(colSums(mt))[1:50],order(colSums(mt))[(ncol(mt)-49):ncol(mt)])]
order(c(pc1_opc90F_sample,pc1_opc90M_sample))
png(filename = "./plots/opc90_pc_clustergram_load_subset_max.png",width = 1400,height = 1000,res = 100)
pheatmap(mat = mt[order(c(pc1_opc90F_sample,pc1_opc90M_sample)),],cluster_cols = F,gaps_col = rep(50,3),color = rev(brewer.pal(11,"RdBu")),main = "90 dpi - 10cRNA-seq projections - PC1",border_color = NA,breaks = seq(-4,4,length.out = 12),annotation_row = ann,cluster_rows = T,scale = "none")
dev.off()