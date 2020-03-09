source("./opcBulk/import_opcBulk.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

bulk_pc <- cbind(bulk$tpm[,1:9],log2(bulk$tpm[,10:ncol(bulk$tpm)] / 100 + 1))

bulk_pc$symbol <- orthologs$HGNC.symbol[match(x = bulk_pc$symbol,table = orthologs$MGI.symbol)]

bulk_pc <- bulk_pc[match(x = pc$Gene,table = bulk_pc$symbol),]

bulk_pc1 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
bulk_pc2 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

bulk_projection1 <- colSums(x = bulk_pc1,na.rm = T)
bulk_projection2 <- colSums(x = bulk_pc2,na.rm = T)

pdf(file = "./plots/pca_subtype_bulk150.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = bulk_projection1[bulk$info$day == 150],
     y = bulk_projection2[bulk$info$day == 150],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 150] == "WT",1,16),
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
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pca_subtype_bulk12.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = bulk_projection1[bulk$info$day == 12],
     y = bulk_projection2[bulk$info$day == 12],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 12] == "WT",1,16),
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
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pca_subtype_bulk90.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = bulk_projection1[bulk$info$day == 90],
     y = bulk_projection2[bulk$info$day == 90],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 90] == "WT",1,16),
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
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()