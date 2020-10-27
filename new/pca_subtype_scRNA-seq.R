
source("./functions/normalizeTPM.R")
source("./functions/pca_subtype.R")

load("./external/gene_detection_rate_processed.RData")

# Plot scRNA-seq -----
sc <- list()
sc$rsem <- rsemcounts_opc_gse75330
sc$tpm <- normalizeTPM(rsem = sc$rsem,index_counts = 10:ncol(sc$rsem))
sc$pc <- pca_subtype(x = sc)

pdf(file = "./plots/single_cell_pca_subtype_75330.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = sc$pc$projection1,
     y = sc$pc$projection2,
     pch = 16,
     xlim = c(-60,60),ylim = c(0,140),
     frame = F,
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Proliferation",main = "scRNA-seq - GSE75330",font.main = 1,cex.main = 8/7,
     axes = F)
lines(x = c(0,0),y = c(0,130))
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2))
lines(x = c(5,52),y = rep(-10.5,2))
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()


sc <- list()
sc$rsem <- rsemcounts_opc_gse60361
sc$tpm <- normalizeTPM(rsem = sc$rsem,index_counts = 10:ncol(sc$rsem))
sc$pc <- pca_subtype(x = sc)

pdf(file = "./plots/single_cell_pca_subtype_60361.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = sc$pc$projection1,
     y = sc$pc$projection2,
     pch = 16,
     xlim = c(-60,60),ylim = c(0,140),
     frame = F,
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Proliferation",main = "scRNA-seq - GSE60361",font.main = 1,cex.main = 8/7,
     axes = F)
lines(x = c(0,0),y = c(0,130))
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2))
lines(x = c(5,52),y = rep(-10.5,2))
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()


# Plot virtual 10c samples from scRNA-seq -----
set.seed(0)
tensc <- list()
tensc$rsem <- rsemcounts_opc_gse75330[,1:9]
for (i in 1:56) {
  iCells <- sample(x = 10:ncol(rsemcounts_opc_gse75330),size = 10,replace = TRUE)
  tensc$rsem <- cbind(tensc$rsem,rowSums(rsemcounts_opc_gse75330[,iCells]))
}

tensc$tpm <- normalizeTPM(rsem = tensc$rsem,index_counts = 10:ncol(tensc$rsem))

tensc$pc <- pca_subtype(x = tensc)

pdf(file = "./plots/single_cell_pca_subtype_75330_virtual10c.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = tensc$pc$projection1,
     y = tensc$pc$projection2,
     pch = 16,
     xlim = c(-60,60),ylim = c(0,140),
     frame = F,
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Proliferation",main = "Virtual 10c samples from scRNA-seq",font.main = 1,cex.main = 8/7,
     axes = F)
lines(x = c(0,0),y = c(0,130))
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2))
lines(x = c(5,52),y = rep(-10.5,2))
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()


set.seed(0)
tensc <- list()
tensc$rsem <- rsemcounts_opc_gse60361[,1:9]
for (i in 1:56) {
  iCells <- sample(x = 10:ncol(rsemcounts_opc_gse60361),size = 10,replace = TRUE)
  tensc$rsem <- cbind(tensc$rsem,rowSums(rsemcounts_opc_gse60361[,iCells]))
}

tensc$tpm <- normalizeTPM(rsem = tensc$rsem,index_counts = 10:ncol(tensc$rsem))

tensc$pc <- pca_subtype(x = tensc)

pdf(file = "./plots/single_cell_pca_subtype_60361_virtual10c.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = tensc$pc$projection1,
     y = tensc$pc$projection2,
     pch = 16,
     xlim = c(-60,60),ylim = c(0,140),
     frame = F,
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Proliferation",main = "Virtual 10c samples from scRNA-seq",font.main = 1,cex.main = 8/7,
     axes = F)
lines(x = c(0,0),y = c(0,130))
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2))
lines(x = c(5,52),y = rep(-10.5,2))
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()