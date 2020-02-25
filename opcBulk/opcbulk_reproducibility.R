# DESeq2 on bulk RNA-sequencing data

library(pheatmap)
library(RColorBrewer)
library(ggbiplot)
library(DESeq2)

source("./functions/normalizeTPM.R")

# File paths
bulk.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/rsem_opcBulk.csv"
bulk.info.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/info_opcBulk.csv"

# Import
bulk <- read.csv(bulk.path,stringsAsFactors = F)
bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)

# Assign row names
row.names(bulk.info) <- bulk.info$name

# Keep only tdT+ samples
bulk <- bulk[,c(1:9,9+which(bulk.info$celltype == "tdTpositive"))]
bulk.info <- bulk.info[bulk.info$celltype == "tdTpositive",]

# Remove outliers
outliers <- c("opcBulk_90dpi_CKO_20647_tdTpositive_GAGCT","opcBulk_90dpi_WT_21167_tdTpositive")
bulk <- bulk[,!(names(bulk) %in% outliers)]
bulk.info <- bulk.info[!(bulk.info$name %in% outliers),]

# Rename samples if they're from the same mouse as an outlier
for (iOutlier in outliers) {
  if (length(strsplit(x = iOutlier,split = "_")[[1]]) == 6) {
    iRootName <- paste(strsplit(x = iOutlier,split = "_")[[1]][-6],collapse = "_")
    names(bulk)[which(grepl(iRootName,names(bulk)))] <- iRootName
    bulk.info$name[which(grepl(iRootName,bulk.info$name))] <- iRootName
  }
}

bulk <- bulk[,c(1:9,9+which(bulk.info$day == 150 & bulk.info$genotype == "CKO"))]
bulk.info <- bulk.info[which(bulk.info$day == 150 & bulk.info$genotype == "CKO"),]

bulk_tpm <- normalizeTPM(bulk,10:ncol(bulk))
bulk_tpm_log2 <- cbind(bulk_tpm[,1:9],log2(bulk_tpm[,10:ncol(bulk_tpm)] + 1))

bulk_tpm_log2_filter <- bulk_tpm_log2[rowSums(bulk_tpm_log2[,10:ncol(bulk_tpm_log2)] > 0) > 0,]

as.numeric(cor.test(bulk_tpm_log2_filter[,16],bulk_tpm_log2_filter[,17])$estimate)

pdf(file = "plots/linear_amplification_reproducible_example.pdf",width = 8,height = 8)
plot(bulk_tpm_log2_filter[,16],bulk_tpm_log2_filter[,17],pch = 20,lwd=0,col = "#00000033",axes=F,xlab = expression("Replicate 1 expression in Log"[2]*"(TPM + 1)"),
                                                                                           ylab = expression("Replicate 2 expression in Log"[2]*"(TPM + 1)"))
text(x = 10,y = 13.5,expression("r = 0.980"),adj=0)
axis(1)
axis(2,las = 2)
dev.off()

png(file = "plots/linear_amplification_reproducible_example.png",width = 2000,height = 2000,res = 400)
par(mar=c(4.1,4.1,0.5,0.5))
plot(bulk_tpm_log2_filter[,16],bulk_tpm_log2_filter[,17],pch = 20,lwd=0,col = "#00000022",axes=F,xlab = expression("Replicate 1 expression in Log"[2]*"(TPM + 1)"),
     ylab = expression("Replicate 2 expression in Log"[2]*"(TPM + 1)"),
     xlim = c(0,15),ylim = c(0,15))
text(x = 11,y = 14.5,expression(italic("r")*" = 0.980"),adj=0)
axis(1)
axis(2,las = 2)
dev.off()

rvals <- unlist(lapply(seq(from = 10,to = 22,by = 2),function(i) as.numeric(cor.test(bulk_tpm_log2_filter[,i],bulk_tpm_log2_filter[,i+1])$estimate)))

rvals_rand <- unlist(lapply(seq(from = 10,to = 22,by = 2),function(i) as.numeric(cor.test(sample(bulk_tpm_log2_filter[,i]),sample(bulk_tpm_log2_filter[,i+1]))$estimate)))

pdf(file = "plots/linear_amplification_reproducible_rvalues.pdf",width = 2,height = 4)
par(mar=c(1,4.5,1,0.5))
plot(rep(0,7),rvals,ylim = c(0,1),axes=F,xlab = "",ylab = expression("Pearson correlation"),frame = F,xaxs = "i")
axis(2,las=2,xpd=T)
axis(1,label=F,tcl=0)
dev.off()

pdf(file = "plots/linear_amplification_reproducible_rvalues_rand.pdf",width = 2,height = 4)
par(mar=c(1,4.5,1,0.5))
plot(-3:3,rvals,xlim = c(-15,15),ylim = c(-0.2,1.1),axes=F,xlab = "",ylab = expression("Pearson correlation"),frame = F,yaxs = "i",xaxs = "i")
points(-3:3,rvals_rand,col = "#888888")
axis(2,las=2,xpd=T)
axis(1,label=F,tcl=0)
dev.off()


all_combos <- combn(10:23,2)[,-c(1,26,47,64,77,86,91)]
rvals_other <- unlist(apply(all_combos,2,function(i) as.numeric(cor.test(bulk_tpm_log2_filter[,i[1]],bulk_tpm_log2_filter[,i[2]])$estimate)))
pdf(file = "plots/linear_amplification_reproducible_rvalues_other.pdf",width = 2,height = 4)
par(mar=c(1,4.5,1,0.5))
plot(seq(-3,3,length.out = 84),rvals_other,col = "#888888",xlim = c(-15,15),ylim = c(0,1.1),axes=F,xlab = "",ylab = expression("Pearson correlation"),frame = F,yaxs = "i",xaxs = "i")
points(-3:3,rvals)
axis(2,las=2,xpd=T)
axis(1,label=F,tcl=0)
dev.off()

