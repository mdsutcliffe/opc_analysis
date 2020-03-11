# Bulk LCM reproducibility

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

pdf(file = "plots/linear_amplification_reproducibility_correlation.pdf",width = 1.5625,height = 3.125,pointsize = 7)
par(mai = c(0.5,0.5,0,0),mgp = c(1.8,0.6,0))
plot(x = rep(x_scatter,2),y = c(random_cor,bulk_cor),
     pch = c(rep(1,length(x_scatter)),rep(16,length(x_scatter))),
     lwd = 0.5,
     xlim = c(0.7,2.7),
     ylim = c(-0.2,1),
     xaxs = "i",
     axes = F,
     xlab = "",
     ylab = "Pearson correlation, R")
axis(side = 2,las = 1,lwd = 0.5)
text(x = 1.2,y = median(bulk_cor),labels = "Tumor replicates",adj = c(0,0.5),cex = 6/7)
text(x = 1.2,y = median(random_cor),labels = "Randomized",adj = c(0,0.5),cex = 6/7)
dev.off()

pdf(file = "plots/linear_amplification_reproducibility_example.pdf",width = 3.125,height = 3.125,pointsize = 7)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(bulk_log2[,iMedian * 2 - 1],bulk_log2[,iMedian * 2],
     pch = 16,
     col = "#00000022",
     frame = F,
     axes = F,
     xlim = c(0,15),
     ylim = c(0,15),
     xlab = "Replicate 1 expression in Log2(TPM + 1)",
     ylab = "Replicate 2 expression in Log2(TPM + 1)")
text(x = 12,y = 6,"R = 0.980")
axis(side = 1,lwd = 0.5)
axis(side = 2,las = 2,lwd = 0.5)
dev.off()

