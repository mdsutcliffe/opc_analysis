# Run UMAP on 90 dpi 10cRNA-seq

library(uwot)

source("./opc90/import_opc90.R")
source("./opc90/RHEGs_opc90.R")

opc90F <- read.table(file = "./data/opc90_true/cleancounts_female_true.tsv")
opc90M <- read.table(file = "./data/opc90_true/cleancounts_male_true.tsv")

umap_opc90 <- cbind(data.frame(row.names = opc90$log2$symbol),opc90$log2[,9+which(opc90$info$type == "ten-cell")])

umap_opc90_rheg <- umap_opc90[opc90_rheg,]

embedding_opc90_rheg <- uwot::umap(X = t(umap_opc90_rheg))

set.seed(90)
embedding_opc90 <- uwot::umap(X = t(umap_opc90))

pdf(file = "./plots/UMAP_opc90.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = embedding_opc90,
     pch = ifelse(opc90$info$mouse[opc90$info$type == "ten-cell"] == "F7460",1,
                  ifelse(opc90$info$mouse[opc90$info$type == "ten-cell"] == "F6340",16,
                         ifelse(opc90$info$mouse[opc90$info$type == "ten-cell"] == "M8170_1",0,15))),
     lwd = 0.5,
     xlim = c(-3,3),
     ylim = c(-2,2),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5,las = 1)
axis(side = 2,lwd = 0.5,las = 1)
dev.off()