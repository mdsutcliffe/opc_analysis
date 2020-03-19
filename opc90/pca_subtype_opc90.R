# PCA subtype analysis from PMID 28697342

source("./opc90/import_opc90.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc90_pc <- cbind(opc90$tpm[,1:9],log2(opc90$tpm[,10:ncol(opc90$tpm)] / 100 + 1))

opc90_pc$symbol <- orthologs$HGNC.symbol[match(x = opc90_pc$symbol,table = orthologs$MGI.symbol)]

opc90_pc <- opc90_pc[match(x = pc$Gene,table = opc90_pc$symbol),]

opc90_pc1 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc90_pc2 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc90_projection1 <- colSums(x = opc90_pc1,na.rm = T)
opc90_projection2 <- colSums(x = opc90_pc2,na.rm = T)

pdf(file = "./plots/pca_subtype_opc90.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = opc90_projection1[opc90$info$type == "ten-cell"],
     y = opc90_projection2[opc90$info$type == "ten-cell"],
     pch = 16,
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5/0.75,axes = F)
lines(x = c(0,0),y = c(0,140),lty = 2,lwd = 0.5/0.75)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2),lwd = 0.5/0.75)
lines(x = c(5,52),y = rep(-10.5,2),lwd = 0.5/0.75)
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()