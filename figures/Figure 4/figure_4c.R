# Figure 4C

source("./import/import_opc90.R")
source("./functions/pca_subtype.R")

opc90$pc <- pca_subtype(opc90)

pdf(file = "./Figure 4/figure_4c.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = opc90$pc$projection1[opc90$info$type == "ten-cell"],
     y = opc90$pc$projection2[opc90$info$type == "ten-cell"],
     pch = 16,
     xlim = c(-60,60),ylim = c(0,140),
     frame = F,
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Proliferation",main = "10cRNA-seq 90 dpi samples",font.main = 1,cex.main = 8/7,
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