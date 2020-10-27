library(pheatmap)
library(RColorBrewer)

# Intersection -----

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

hallmark <- read.table("./external/hallmark_emt.txt",header = T,stringsAsFactors = F)[,1]

intersect(hallmark,pc$Gene)

# 12 dpi -----

source("./import/import_opc12.R")

x <- opc12

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

x_pc <- cbind(x$tpm[,1:9],log2(x$tpm[,10:ncol(x$tpm)] / 100 + 1))

x_pc$symbol <- orthologs$HGNC.symbol[match(x = x_pc$symbol,table = orthologs$MGI.symbol)]

x_pc <- x_pc[match(x = pc$Gene,table = x_pc$symbol),]

x_pc1 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
x_pc2 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

x_projection1 <- colSums(x = x_pc1,na.rm = T)
x_projection2 <- colSums(x = x_pc2,na.rm = T)

row.names(x_pc1) <- x_pc$symbol

png(filename = "./plots/pca_heatmap_expanded_12dpi.png",width = 8000,height = 1000,units = "px",res = 100)
pheatmap(mat = t(x_pc1[complete.cases(x_pc1),]),cluster_cols = F,color = rev(brewer.pal(n = 11,name = "RdBu")))
dev.off()

# 90 dpi -----

source("./import/import_opc90.R")

x <- opc90

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

x_pc <- cbind(x$tpm[,1:9],log2(x$tpm[,10:ncol(x$tpm)] / 100 + 1))

x_pc$symbol <- orthologs$HGNC.symbol[match(x = x_pc$symbol,table = orthologs$MGI.symbol)]

x_pc <- x_pc[match(x = pc$Gene,table = x_pc$symbol),]

x_pc1 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
x_pc2 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

x_projection1 <- colSums(x = x_pc1,na.rm = T)
x_projection2 <- colSums(x = x_pc2,na.rm = T)

row.names(x_pc1) <- x_pc$symbol

png(filename = "./plots/pca_heatmap_expanded_90dpi.png",width = 8000,height = 1000,units = "px",res = 100)
pheatmap(mat = t(x_pc1[complete.cases(x_pc1),]),cluster_cols = F,color = rev(brewer.pal(n = 11,name = "RdBu")))
dev.off()