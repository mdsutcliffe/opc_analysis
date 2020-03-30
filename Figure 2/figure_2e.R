# Figure 2E

# source("./functions/opcBulk_DESeq2.R")
load("./build/opcBulk_DESeq2_results.RData")

de150 <- bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj),c("log2FoldChange","padj")]

pdf(file = "./Figure 2/figure_2e.pdf",width = 3.125*1.25,height = 2,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.25,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = de150$log2FoldChange[de150$padj >= 0.05],
     y = -log10(x = de150$padj[de150$padj >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbd22",
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,100,20),labels = c(NA,seq(20,100,20)),pos = 0,las = 1,lwd = 0.5/0.75)
points(x = de150$log2FoldChange[de150$padj < 0.05 & de150$log2FoldChange < 0],
       y = -log10(x = de150$padj[de150$padj < 0.05 & de150$log2FoldChange < 0]),
       pch = 16,cex = 0.5,
       col = "#2166ac22")
points(x = de150$log2FoldChange[de150$padj < 0.05 & de150$log2FoldChange > 0],
       y = -log10(x = de150$padj[de150$padj < 0.05 & de150$log2FoldChange > 0]),
       pch = 16,cex = 0.5,
       col = "#b2182b22")
text(x = 0.5,y = 100,labels = "-Log10(p-value)",adj = c(0,0.5),xpd = T)
dev.off()