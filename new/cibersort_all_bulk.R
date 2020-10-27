# Figure 2C

source("./import/import_opcBulk.R")
source("./functions/normalizeTPM.R")
# source("./functions/signatureMatrix.R")
# 
# signature <- signatureMatrix(geneList = bulk$tpm$symbol)
# 
# commonGenes <- intersect(bulk$tpm$symbol,signature$tpm$symbol)
# 
# bulk$tpm <- normalizeTPM(rsem = bulk$rsem[match(commonGenes,bulk$rsem$symbol),],
#                          index_counts = 10:ncol(bulk$rsem))
# 
# write.table(x = signature$matrix,file = "./Figure 2/signature.txt",
#             append = F,quote = F,sep = "\t",row.names = F,col.names = T)
# write.table(x = bulk$tpm[,c(3,10:ncol(bulk$tpm))],
#             file = "./Figure 2/mixture_opcBulk.txt",
#             append = F,quote = F,sep = "\t",row.names = F,col.names = T)

# Run Cibersort HERE

f.cibersort <- "./figures/Figure 2/CIBERSORTx_opcBulk.txt"

cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")

row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)] / cibersort$Absolute.score..sig.score.

pdf(file = "./plots/cibersort_bulk_12dpi_wt.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 12 dpi WT",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 12 & bulk$info$genotype == "WT",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()


pdf(file = "./plots/cibersort_bulk_12dpi_cko.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 12 dpi CKO",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 12 & bulk$info$genotype == "CKO",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()



pdf(file = "./plots/cibersort_bulk_90dpi_wt.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 90 dpi WT",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 90 & bulk$info$genotype == "WT",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()


pdf(file = "./plots/cibersort_bulk_90dpi_cko.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 90 dpi CKO",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 90 & bulk$info$genotype == "CKO",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()



pdf(file = "./plots/cibersort_bulk_150dpi_wt.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 150 dpi WT",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 150 & bulk$info$genotype == "WT",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()


pdf(file = "./plots/cibersort_bulk_150dpi_cko.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
graphics::plot(x = c(),y = c(),
               xlim = c(0,80),ylim = c(0.5,6.5),
               xlab = "Relative proportion",ylab = NA,main = "Sample unmixing - 150 dpi CKO",font.main = 1,cex.main = 8/7,
               xaxs = "i",yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort[bulk$info$day == 150 & bulk$info$genotype == "CKO",rev(c("OPC","Astrocyte","Neuron","MO","Endothelial","Microglia"))] * 100,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  lty = 1,
                  range = 2,
                  axes = F,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()














library(pheatmap)
library(RColorBrewer)

f.cibersort <- "./figures/Figure 2/CIBERSORTx_opcBulk.txt"

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 12 & bulk$info$genotype == "WT",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_12dpi_wt.pdf",
              width = 3.125,
              height = 3.125)

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 12 & bulk$info$genotype == "CKO",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_12dpi_cko.pdf",
              width = 3.125,
              height = 3.125)

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 90 & bulk$info$genotype == "WT",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_90dpi_wt.pdf",
              width = 3.125,
              height = 3.125)

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 90 & bulk$info$genotype == "CKO",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_90dpi_cko.pdf",
              width = 3.125,
              height = 3.125)

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 150 & bulk$info$genotype == "WT",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_150dpi_wt.pdf",
              width = 3.125,
              height = 3.125)

res <- read.table(file = f.cibersort,header = T,sep = "\t",row.names = 1)

res <- res[bulk$info$day == 150 & bulk$info$genotype == "CKO",] / max(res[,1:6])

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),breaks = seq(0,1,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_bulk_150dpi_cko.pdf",
              width = 3.125,
              height = 3.125)