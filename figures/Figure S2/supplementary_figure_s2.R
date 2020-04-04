# Supplementary Figure S2

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

bulk_tpm <- normalizeTPM(rsem = bulk,index_counts = 10:ncol(bulk))
bulk_log2 <- cbind(bulk_tpm[,1:9],log2(bulk_tpm[,10:ncol(bulk_tpm)] + 1))
bulk_log2 <- cbind(data.frame(row.names = bulk_log2$symbol),bulk_log2[,9 + which(bulk.info$genotype == "CKO" & bulk.info$day == 150)])
bulk_log2 <- bulk_log2[rowSums(bulk_log2) > 0,]

bulk_cor <- sapply(X = seq(1,14,2),FUN = function(i) cor(bulk_log2[,i],bulk_log2[,i + 1]))
random_cor <- sapply(X = seq(1,14,2),FUN = function(i) cor(sample(bulk_log2[,i]),sample(bulk_log2[,i + 1])))

iMedian <- match(median(bulk_cor),bulk_cor)

x_scatter <- seq(from = 0.9,to = 1.1,length.out = length(bulk_cor))

pdf(file = "./Figure S2/figure_s2a.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(bulk_log2[,iMedian * 2 - 1],bulk_log2[,iMedian * 2],
     pch = 16,
     col = "#00000022",
     cex = 0.5,
     lwd = 0.5/0.75,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     xlim = c(0,16),
     ylim = c(0,16),
     xlab = "Replicate 1 expression in Log2(TPM + 1)",
     ylab = "Replicate 2 expression in Log2(TPM + 1)",
     xpd = T)
text(x = 12,y = 6,paste("R =",sprintf("%.3f",bulk_cor[iMedian])))
axis(side = 1,at = seq(0,16,4),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,16,4),lwd = 0.5/0.75,las = 2)
dev.off()

pdf(file = "./Figure S2/figure_s2b.pdf",width = 1.953125,height = 3.125,pointsize = 7,bg = "white")
par(mai = c(0.5,0.5,0,0),mgp = c(2.3,0.8,0))
plot(x = NULL,y = NULL,
     xlim = c(0.4,2.6),ylim = c(-0.2,1),
     xaxs = "i",yaxs = "i",
     axes = F,
     xlab = NA,ylab = "Pearrson correlation, R")
graphics::boxplot(x = cbind(random_cor,bulk_cor),
                  xaxs = "i",yaxs = "i",
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,at = 0:3,labels = c(NA,"Randomized","Tumor\nreplicates",NA),padj = 1,lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75)
dev.off()