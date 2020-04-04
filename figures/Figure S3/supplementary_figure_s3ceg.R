# Supplementary Figure S3 C+E+G

source("./import/import_opcBulk.R")

f.cibersort <- "./Figure S3/CIBERSORTx_FigureS3_endothelial.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.
names(cibersort_signature)[names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"

cibersort_bulk150 <- cibersort_signature[c(bulk$info$day == 150 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]
cibersort_bulk12 <- cibersort_signature[c(bulk$info$day == 12 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]
cibersort_bulk90 <- cibersort_signature[c(bulk$info$day == 90 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]

pdf(file = "./Figure S3/figure_s3c.pdf",width = 3.125,height = 2.34375,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.5,0.75,0.25,0.25),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Bulk 150 dpi control samples",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk150[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

pdf(file = "./Figure S3/figure_s3e.pdf",width = 3.125,height = 2.34375,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.5,0.75,0.25,0.25),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Bulk 12 dpi control samples",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk12[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

pdf(file = "./Figure S3/figure_s3g.pdf",width = 3.125,height = 2.34375,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.5,0.75,0.25,0.25),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Bulk 90 dpi control samples",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk90[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()