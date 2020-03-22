# Figure S7

source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))
opc90_pc <- cbind(opc90$tpm[,1:9],log2(opc90$tpm[,10:ncol(opc90$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]
opc90_pc$symbol <- orthologs$HGNC.symbol[match(x = opc90_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]
opc90_pc <- opc90_pc[match(x = pc$Gene,table = opc90_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc90_pc1 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)

opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))
opc90_pc2 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
opc90_projection1 <- colSums(x = opc90_pc1,na.rm = T)

opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)
opc90_projection2 <- colSums(x = opc90_pc2,na.rm = T)

x = seq(-30,30,2.5)

y_female <- hist(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],breaks = x,plot = F)
y_male <- hist(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],breaks = x,plot = F)

y_both <- sapply(X = 1:(length(x)-1),FUN = function(i) min(c(y_female$counts[i],y_male$counts[i])))
y_female_only <- y_female$counts - y_both
y_male_only <- y_male$counts - y_both

pdf(file = "./plots/figure_s7a.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = NA,y = NA,
     xlim = range(x),ylim = c(0,8),
     frame = F,axes = F,
     xlab = NA,ylab = "Frequency",
     xaxs = "i",yaxs = "i")
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_female_only,col = "#4dac26",border = NA)
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_male_only,col = "#d01c8b",border = NA)
rect(xleft = x,ybottom = 0,xright = x + 2.5,ytop = y_both,col = "#8e6458",border = NA)
legend(x = "topright",legend = c("female","male"),fill = c("#4dac26","#d01c8b"),border = F,bty = "n")
axis(side = 1,at = seq(-30,30,10),labels = c(-30,NA,NA,0,NA,NA,30),las = 1,mgp = c(1.6,0.6,0),lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75)
lines(x = c(-5,-52)/2,y = rep(-10.5,2)*8/140,lwd = 0.5/0.75)
lines(x = c(5,52)/2,y = rep(-10.5,2)*8/140,lwd = 0.5/0.75)
text(x = -27.5/2,y = -15*8/140,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -15*8/140,labels = "Mesenchymal",cex = 6/7)
dev.off()

y_female <- hist(x = opc90_projection1[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],breaks = x,plot = F)
y_male <- hist(x = opc90_projection1[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],breaks = x,plot = F)

y_both <- sapply(X = 1:(length(x)-1),FUN = function(i) min(c(y_female$counts[i],y_male$counts[i])))
y_female_only <- y_female$counts - y_both
y_male_only <- y_male$counts - y_both

pdf(file = "./plots/figure_s7b.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = NA,y = NA,
     xlim = range(x),ylim = c(0,8),
     frame = F,axes = F,
     xlab = NA,ylab = "Frequency",
     xaxs = "i",yaxs = "i")
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_female_only,col = "#4dac26",border = NA)
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_male_only,col = "#d01c8b",border = NA)
rect(xleft = x,ybottom = 0,xright = x + 2.5,ytop = y_both,col = "#8e6458",border = NA)
legend(x = "topright",legend = c("female","male"),fill = c("#4dac26","#d01c8b"),border = F,bty = "n")
axis(side = 1,at = seq(-30,30,10),labels = c(-30,NA,NA,0,NA,NA,30),las = 1,mgp = c(1.6,0.6,0),lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75)
lines(x = c(-5,-52)/2,y = rep(-10.5,2)*8/140,lwd = 0.5/0.75)
lines(x = c(5,52)/2,y = rep(-10.5,2)*8/140,lwd = 0.5/0.75)
text(x = -27.5/2,y = -15*8/140,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -15*8/140,labels = "Mesenchymal",cex = 6/7)
dev.off()
