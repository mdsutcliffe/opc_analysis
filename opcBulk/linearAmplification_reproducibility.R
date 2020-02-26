# Bulk LCM reproducibility

source("./functions/normalizeTPM.R")

# File paths
bulk.path <- "./data/rsem_opcBulk.csv"
bulk.info.path <- "./data/info_opcBulk.csv"

# Import
bulk <- read.csv(bulk.path,stringsAsFactors = F)
bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)

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

bulk_tpm <- normalizeTPM(rsem = bulk,index_counts = 10:ncol(bulk))
bulk_log2 <- cbind(bulk_tpm[,1:9],log2(bulk_tpm[,10:ncol(bulk_tpm)] + 1))
bulk_log2 <- cbind(data.frame(row.names = bulk_log2$symbol),bulk_log2[,9 + which(bulk.info$genotype == "CKO" & bulk.info$day == 150)])
bulk_log2 <- bulk_log2[rowSums(bulk_log2) > 0,]

bulk_cor <- sapply(X = seq(1,14,2),FUN = function(i) cor(bulk_log2[,i],bulk_log2[,i + 1]))
random_cor <- sapply(X = seq(1,14,2),FUN = function(i) cor(sample(bulk_log2[,i]),sample(bulk_log2[,i + 1])))

iMedian <- match(median(bulk_cor),bulk_cor)

x_scatter <- seq(from = 0.9,to = 1.1,length.out = length(bulk_cor))

pdf(file = "plots/linear_amplification_reproducibility_correlation.pdf",width = 1.75,height = 3,pointsize = 6)
par(mar = c(3.6,3.6,0.8,0.8))
plot(x = x_scatter,y = random_cor,
     pch = 1,
     lwd = 0.5,
     xlim = c(0.7,2.7),
     ylim = c(-0.2,1),
     axes = F,
     xlab = "",
     ylab = "")
points(x = x_scatter,y = bulk_cor,pch = 16,lwd = 0.5)
axis(side = 2,las = 1)
legend(x = "topright",legend = c("Tumor replicates","Randomized"),pch = c(16,1),pt.lwd = c(0.5,0.5),box.lwd = 0.5)
title(ylab = expression("Pearson correlation, "*italic("R")),mgp = c(2.5,1,0))
dev.off()

pdf(file = "plots/linear_amplification_reproducibility_example.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(3.6,3.6,0.8,0.8))
plot(bulk_log2[,iMedian * 2 - 1],bulk_log2[,iMedian * 2],
     pch = 20,
     lwd = 0,
     col = "#00000033",
     frame = F,
     axes = F,
     las = 1,
     xlim = c(0,15),
     ylim = c(0,15),
     xlab = "",
     ylab = "")
text(x = 13,y = 10,expression(italic("R")*" = 0.980"))
axis(side = 1)
axis(side = 2,las = 2)
title(xlab = expression("Replicate 1 expression in Log"[2]*"(TPM + 1)"),
      ylab = expression("Replicate 2 expression in Log"[2]*"(TPM + 1)"),mgp = c(2.5,1,0))
dev.off()

