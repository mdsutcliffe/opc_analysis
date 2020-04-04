# Supplementary Figure S8

source("./import/import_opc12.R")
source("./import/import_opc90.R")
source("./functions/pca_subtype.R")

opc12$pc <- pca_subtype(opc12)
opc90$pc <- pca_subtype(opc90)

x = seq(-30,30,2.5)

y_female <- hist(x = opc12$pc$projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],breaks = x,plot = F)
y_male <- hist(x = opc12$pc$projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],breaks = x,plot = F)

y_both <- sapply(X = 1:(length(x)-1),FUN = function(i) min(c(y_female$counts[i],y_male$counts[i])))
y_female_only <- y_female$counts - y_both
y_male_only <- y_male$counts - y_both

pdf(file = "./Figure S8/figure_s8a.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = range(x),ylim = c(0,8),
     frame = F,axes = F,
     xlab = NA,ylab = "Frequency",main = "10cRNA-seq 12 dpi samples",font.main = 1,cex.main = 8/7,
     xaxs = "i",yaxs = "i")
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_female_only,col = "#4dac26",border = NA)
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_male_only,col = "#d01c8b",border = NA)
rect(xleft = x,ybottom = 0,xright = x + 2.5,ytop = y_both,col = "#8e6458",border = NA)
legend(x = "topright",legend = c("female","male"),fill = c("#4dac26","#d01c8b"),border = F,bty = "n")
axis(side = 1,at = seq(-30,30,10),labels = c(-30,NA,NA,0,NA,NA,30),las = 1,mgp = c(1.6,0.6,0),lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75)
lines(x = c(-5,-52)/2,y = rep(-8.5,2)*8/140,lwd = 0.5/0.75)
lines(x = c(5,52)/2,y = rep(-8.5,2)*8/140,lwd = 0.5/0.75)
text(x = -27.5/2,y = -12.75*8/140,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -12.75*8/140,labels = "Mesenchymal",cex = 6/7)
dev.off()

y_female <- hist(x = opc90$pc$projection1[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],breaks = x,plot = F)
y_male <- hist(x = opc90$pc$projection1[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],breaks = x,plot = F)

y_both <- sapply(X = 1:(length(x)-1),FUN = function(i) min(c(y_female$counts[i],y_male$counts[i])))
y_female_only <- y_female$counts - y_both
y_male_only <- y_male$counts - y_both

pdf(file = "./Figure S8/figure_s8b.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = range(x),ylim = c(0,8),
     frame = F,axes = F,
     xlab = NA,ylab = "Frequency",main = "10cRNA-seq 90 dpi samples",font.main = 1,cex.main = 8/7,
     xaxs = "i",yaxs = "i")
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_female_only,col = "#4dac26",border = NA)
rect(xleft = x,ybottom = y_both,xright = x + 2.5,ytop = y_both + y_male_only,col = "#d01c8b",border = NA)
rect(xleft = x,ybottom = 0,xright = x + 2.5,ytop = y_both,col = "#8e6458",border = NA)
legend(x = "topright",legend = c("female","male"),fill = c("#4dac26","#d01c8b"),border = F,bty = "n")
axis(side = 1,at = seq(-30,30,10),labels = c(-30,NA,NA,0,NA,NA,30),las = 1,mgp = c(1.6,0.6,0),lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75)
lines(x = c(-5,-52)/2,y = rep(-8.5,2)*8/140,lwd = 0.5/0.75)
lines(x = c(5,52)/2,y = rep(-8.5,2)*8/140,lwd = 0.5/0.75)
text(x = -27.5/2,y = -12.75*8/140,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -12.75*8/140,labels = "Mesenchymal",cex = 6/7)
dev.off()