# Figure 3B

# source("./functions/opcBulk_DESeq2.R")
load("./build/opcBulk_DESeq2_results.RData")

de12 <- bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj),c("log2FoldChange","padj")]


pdf(file = "./Figure 3/figure_3b.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = de12$log2FoldChange[de12$padj >= 0.05],
     y = -log10(x = de12$padj[de12$padj >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbd22",
     xlim = c(-10,10),ylim = c(0,60),
     xlab = "Log2 fold change",ylab = NA,main = "Bulk 12 dpi samples",font.main = 1,cex.main = 8/7,
     frame = F,axes = F,
     xaxs = "i",yaxs = "i")
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,60,10),labels = c(NA,seq(10,60,10)),pos = 0,las = 1,lwd = 0.5/0.75)
points(x = de12$log2FoldChange[de12$padj < 0.05 & de12$log2FoldChange < 0],
       y = -log10(x = de12$padj[de12$padj < 0.05 & de12$log2FoldChange < 0]),
       pch = 16,cex = 0.5,
       col = "#2166ac22")
points(x = de12$log2FoldChange[de12$padj < 0.05 & de12$log2FoldChange > 0],
       y = -log10(x = de12$padj[de12$padj < 0.05 & de12$log2FoldChange > 0]),
       pch = 16,cex = 0.5,
       col = "#b2182b22")
text(x = 0.5,y = 60,labels = "Log10(pvalue)",adj = c(0,0.5),xpd = T)
dev.off()