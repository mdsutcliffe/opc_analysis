# Figure 2D

source("./import/import_opcBulk.R")
source("./functions/pca_subtype.R")

bulk$pc <- pca_subtype(bulk)

pdf(file = "./Figure 2/figure_2d.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = bulk$pc$projection1[bulk$info$day == 150],
     y = bulk$pc$projection2[bulk$info$day == 150],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 150] == "WT",1,16),
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     axes = F)
lines(x = c(0,0),y = c(0,140))
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2))
lines(x = c(5,52),y = rep(-10.5,2))
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
legend(x = "topright",legend = c("Control","Mutant"),pch = c(1,16),bty = "n")
dev.off()