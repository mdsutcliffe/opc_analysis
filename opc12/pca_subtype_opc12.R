# PCA subtype analysis from PMID 28697342

source("./opc12/import_opc12.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)

pdf(file = "./plots/pca_subtype_opc12.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = opc12_projection1[opc12$info$type == "ten-cell"],
     y = opc12_projection2[opc12$info$type == "ten-cell"],
     pch = 16,
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5,axes = F)
lines(x = c(0,0),y = c(0,140),lty = 2,lwd = 0.5)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2),lwd = 0.5)
lines(x = c(5,52),y = rep(-10.5,2),lwd = 0.5)
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()

# Agglomerations
set.seed(0)
nCells <- sample(x = seq(700,1000,10),size = 100,replace = T)
agglomerate_pc <- t(sapply(X = nCells,FUN = function(nC) {
        cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC/10,replace = T)
        
        iOPC12 <- cbind(opc12$tpm[,1:9],log2(rowMeans(opc12$tpm[,cellGroups]) / 100 + 1))
        
        iOPC12$symbol <- orthologs$HGNC.symbol[match(x = iOPC12$symbol,table = orthologs$MGI.symbol)]
        
        iOPC12 <- iOPC12[match(x = pc$Gene,table = iOPC12$symbol),]
        
        iOPC12 <- iOPC12[,10:ncol(iOPC12)]
        
        iPC1 <- pc$PC1 * iOPC12
        iPC2 <- pc$PC2 * (iOPC12 - (pc$PC1 * iOPC12))
        
        iProjection1 <- sum(x = iPC1,na.rm = T)
        iProjection2 <- sum(x = iPC2,na.rm = T)
        
        return(c(iProjection1,iProjection2))
}))
plot(agglomerate_pc)
pdf(file = "./plots/pca_subtype_opc12_agglomerate.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = agglomerate_pc,
     pch = 16,
     col = "#00000022",
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5,axes = F)
lines(x = c(0,0),y = c(0,140),lty = 2,lwd = 0.5)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2),lwd = 0.5)
lines(x = c(5,52),y = rep(-10.5,2),lwd = 0.5)
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()



# Agglomerations
set.seed(0)
nCells <- sample(x = 700:1000,size = 100,replace = T)
agglomerate_pc <- t(sapply(X = nCells,FUN = function(nC) {
        cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC,replace = T)
        
        iOPC12 <- cbind(opc12$tpm[,1:9],log2(rowMeans(opc12$tpm[,cellGroups]/10) / 100 + 1))
        
        iOPC12$symbol <- orthologs$HGNC.symbol[match(x = iOPC12$symbol,table = orthologs$MGI.symbol)]
        
        iOPC12 <- iOPC12[match(x = pc$Gene,table = iOPC12$symbol),]
        
        iOPC12 <- iOPC12[,10:ncol(iOPC12)]
        
        iPC1 <- pc$PC1 * iOPC12
        iPC2 <- pc$PC2 * (iOPC12 - (pc$PC1 * iOPC12))
        
        iProjection1 <- sum(x = iPC1,na.rm = T)
        iProjection2 <- sum(x = iPC2,na.rm = T)
        
        return(c(iProjection1,iProjection2))
}))
plot(agglomerate_pc)
pdf(file = "./plots/pca_subtype_opc12_agglomerate_singlecell.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = agglomerate_pc,
     pch = 16,
     col = "#00000022",
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5,axes = F)
lines(x = c(0,0),y = c(0,140),lty = 2,lwd = 0.5)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2),lwd = 0.5)
lines(x = c(5,52),y = rep(-10.5,2),lwd = 0.5)
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()

# Agglomerations for heatmap
set.seed(0)
nCells <- sample(x = seq(700,1000,10),size = 100,replace = T)
agglomerate_pc <- t(sapply(X = nCells,FUN = function(nC) {
        cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC/10,replace = T)
        
        iOPC12 <- cbind(opc12$tpm[,1:9],log2(rowMeans(opc12$tpm[,cellGroups]) / 100 + 1))
        
        iOPC12$symbol <- orthologs$HGNC.symbol[match(x = iOPC12$symbol,table = orthologs$MGI.symbol)]
        
        iOPC12 <- iOPC12[match(x = pc$Gene,table = iOPC12$symbol),]
        
        iOPC12 <- iOPC12[,10:ncol(iOPC12)]
        
        iPC1 <- pc$PC1 * iOPC12
        iPC2 <- pc$PC2 * (iOPC12 - (pc$PC1 * iOPC12))
        
        iProjection1 <- sum(x = iPC1,na.rm = T)
        iProjection2 <- sum(x = iPC2,na.rm = T)
        
        return(iPC1)
}))

agglomerate_pc_hm <- agglomerate_pc[,complete.cases(bulk_pc1)]
pdf(file = "./plots/pc_opc12_agglomerations_heatmap.pdf",width = 8,height = 4,pointsize = 7,useDingbats = F)
pheatmap(mat = agglomerate_pc_hm[,ncol(agglomerate_pc_hm):1],
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         cluster_rows = F,cluster_cols = F,
         show_rownames = F,show_colnames = F,
         main = "12 dpi - 700-1k cell simulations")
dev.off()
