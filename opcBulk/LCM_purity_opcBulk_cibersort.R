# Verify LCM purity by cibersort

library(RColorBrewer)

source("./opcBulk/import_opcBulk.R")
source("./functions/signatureMatrix.R")
source("./functions/normalizeTPM.R")

signature <- signatureMatrix(geneList = bulk$rsem$symbol,removeNFO = T)
commonGenes <- intersect(bulk$tpm$symbol,signature$tpm$symbol)

bulk$tpm <- normalizeTPM(rsem = bulk$rsem[match(commonGenes,bulk$rsem$symbol),],index_counts = 10:ncol(bulk$rsem))

wt12 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 12 & bulk$info$genotype == "WT"))]
wt150 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"))]

wt12 <- wt12[,c(3,10:ncol(wt12))]
wt150 <- wt150[,c(3,10:ncol(wt150))]

write.table(x = wt12,file = "./temp/mixture_bulk_12_wt.txt",quote = F,sep = "\t",row.names = F)
write.table(x = wt150,file = "./temp/mixture_bulk_150_wt.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort, get results, and place in external folder

f.cibersort_12_wt <- "./external/CIBERSORTx_bulk_12_wt.txt"
f.cibersort_150_wt <- "./external/CIBERSORTx_bulk_150_wt.txt"

cibersort_12_wt <- read.table(file = f.cibersort_12_wt,header = T,sep = "\t")
cibersort_150_wt <- read.table(file = f.cibersort_150_wt,header = T,sep = "\t")

row.names(cibersort_12_wt) <- cibersort_12_wt$Mixture
row.names(cibersort_150_wt) <- cibersort_150_wt$Mixture

cibersort_12_wt <- cibersort_12_wt[,2:(which(names(cibersort_12_wt) == "P.value") - 1)] / cibersort_12_wt$Absolute.score..sig.score.
cibersort_150_wt <- cibersort_150_wt[,2:(which(names(cibersort_150_wt) == "P.value") - 1)] / cibersort_150_wt$Absolute.score..sig.score.

pdf(file = "./plots/bulk_LCM_purity_12_WT_cibersort.pdf",width = 6,height = 6)
par(mar = c(1.8,4,1,8),mgp = c(2.9,1,0),xpd = T)
barplot(t(as.matrix(cibersort_12_wt)),ylab = "Fraction",names.arg = rep(x = "",nrow(cibersort_12_wt)),col = rev(brewer.pal(6,"Set1")),las = 1)
title(xlab = "12 dpi WT bulk samples",mgp = c(0.5,0,0))
legend(x = "topright",legend = names(cibersort_12_wt),pch = 15,col = rev(brewer.pal(6,"Set1")),inset = c(-0.4,0))
dev.off()

pdf(file = "./plots/bulk_LCM_purity_150_WT_cibersort.pdf",width = 6,height = 6)
par(mar = c(1.8,4,1,8),mgp = c(2.9,1,0),xpd = T)
barplot(t(as.matrix(cibersort_150_wt)),ylab = "Fraction",names.arg = rep(x = "",nrow(cibersort_150_wt)),col = rev(brewer.pal(6,"Set1")),las = 1)
title(xlab = "150 dpi WT bulk samples",mgp = c(0.5,0,0))
legend(x = "topright",legend = names(cibersort_150_wt),pch = 15,col = rev(brewer.pal(6,"Set1")),inset = c(-0.4,0))
dev.off()


# Boxplots
pdf(file = "./plots/cibersort_WT_150_boxplot.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
boxplot(x = cibersort_150_wt[,rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))],horizontal = T,las = 1,xlab = "Relative proportion",frame = F,xaxs = "i",yaxs = "i",axes = F,ylim = c(0,0.8),xlim = c(0.7,6.3),par(mar = c(3.1,4.6,0.5,1),mgp = c(2,1,0)),lwd = 0.5)
axis(side = 1,at = seq(0,0.8,0.2),lwd = 0.5)
axis(side = 2,at = seq(0,6)+0.5,labels = F,lwd = 0.5)
axis(side = 2,at = seq(1,6),labels = rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia")),las = 1,tick = F,lwd = 0.5,par(mgp = c(0,0.2,0)))
dev.off()

pdf(file = "./plots/cibersort_WT_12_boxplot.pdf",width = 2.5,height = 2,pointsize = 7,useDingbats = F)
boxplot(x = cibersort_150_wt[,rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))],horizontal = T,las = 1,xlab = "Relative proportion",frame = F,xaxs = "i",yaxs = "i",axes = F,ylim = c(0,0.8),xlim = c(0.7,5.3),par(mar = c(3.1,4.6,0.5,1),mgp = c(2,1,0)),lwd = 0.5)
axis(side = 1,at = seq(0,0.8,0.2),lwd = 0.5)
axis(side = 2,at = seq(0,5)+0.5,labels = F,lwd = 0.5)
axis(side = 2,at = seq(1,5),labels = rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia")),las = 1,tick = F,lwd = 0.5,par(mgp = c(0,0.2,0)))
dev.off()


pdf(file = "./plots/cibersort_WT_150_boxplot_full.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
boxplot(x = cibersort_150_wt[,rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))],
        horizontal = T,las = 1,xlab = "Relative proportion",frame = F,xaxs = "i",yaxs = "i",
        axes = F,ylim = c(0,0.8),xlim = c(0.7,6.2),par(mai = c(0.5,0.5,0,0)),lwd = 0.5)
axis(side = 1,at = seq(0,0.8,0.2),lwd = 0.5)
axis(side = 2,at = seq(0,6)+0.5,labels = F,lwd = 0.5)
axis(side = 2,at = seq(1,6),labels = rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia")),las = 1,tick = F,lwd = 0.5,par(mgp = c(0,0.2,0)))
dev.off()

pdf(file = "./plots/cibersort_bulk150_boxplot.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = c(),y = c(),
               xlim = c(0,0.8),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_150_wt[,rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))],
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

cairo_pdf(filename = "./plots/cibersort_bulk150_boxplot_cairo.pdf",width = 2.25,height = 2.25,pointsize = 7,family = "Arial",)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = c(0,1),y = c(0,1),
               xlim = c(0,0.8),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_150_wt[,rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))],
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lwd = 0.5,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5)
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia")),par(mgp = c(0,0.15,0)))
axis(side = 2,las = 1,lwd = 0.5,at = 0:6 + 0.5,labels = NA)
dev.off()