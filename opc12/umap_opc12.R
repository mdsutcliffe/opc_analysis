
setwd("~")
library(umapr)

setwd("~/Github/opc_analysis/")
source("./opc12/import_opc12.R")
source("./opc12/RHEGs_opc12.R")

opc12F <- read.table(file = "./data/opc12_true/cleancounts_female_true.tsv")
opc12M <- read.table(file = "./data/opc12_true/cleancounts_male_true.tsv")

umap_opc12 <- cbind(data.frame(row.names = opc12$log2$symbol),opc12$log2[,9+which(opc12$info$type == "ten-cell")])

umap_opc12_rheg <- umap_opc12[opc12_rheg,]

embedding_opc12_rheg <- uwot::umap(X = t(umap_opc12_rheg))

set.seed(12)
embedding_opc12 <- uwot::umap(X = t(umap_opc12))

pdf(file = "./plots/umap_opc12_final.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0))
plot(x = embedding_opc12,
     xlab = "UMAP-1",ylab = "UMAP-2",
     pch = ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8520",1,
                  ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "F8519",16,
                         ifelse(opc12$info$mouse[opc12$info$type == "ten-cell"] == "M8516",0,15))),
     frame = F,lwd = 0.5,
     las = 1,
     xaxs = "i",
     yaxs = "i",
     xlim = c(-6,6),ylim = c(-4,4))
dev.off()