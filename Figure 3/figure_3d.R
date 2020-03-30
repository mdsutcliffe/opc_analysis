# Figure 3D

library(uwot)

source('./import/import_opc12.R')

opc12F <- read.table(file = "./data/opc12_true/cleancounts_female_true.tsv")
opc12M <- read.table(file = "./data/opc12_true/cleancounts_male_true.tsv")

umap_opc12 <- cbind(data.frame(row.names = opc12$log2$symbol),opc12$log2[,9+which(opc12$info$type == "ten-cell")])

set.seed(12)
embedding_opc12 <- uwot::umap(X = t(umap_opc12))

pdf(file = "./Figure 3/figure_3d.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = embedding_opc12,
     pch = ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8520",1,
                  ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8519",16,
                         ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "M8516",1,16))),
     col = ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8520","#4dac26",
                  ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8519","#4dac26",
                         ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "M8516","#d01c8b","#d01c8b"))),
     lwd = 0.5/0.75,
     xlim = c(-6,6),ylim = c(-4,4),
     xlab = "UMAP-1",ylab = "UMAP-2",main = "10cRNA-seq 12 dpi samples",font.main = 1,cex.main = 8/7,
     frame = F,axes = F,
     xaxs = "i",yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
points(x = c(4,4,4.667,4.667),y = c(4,4.667,4,4.667)*4/6,pch = c(1,1,16,16),col = c("#4dac26","#d01c8b","#4dac26","#d01c8b"))
text(x = c(5,5),y = c(4,4.667)*4/6,labels = c("Female","Male"),adj = c(0,0.5),cex = 6/7)
text(x = c(4,4.667),y = c(5,5)*4/6,labels = c("#1","#2"),adj = c(0.5,0),cex = 6/7)
dev.off()
