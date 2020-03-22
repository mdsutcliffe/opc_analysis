# Figure S5

library(pheatmap)
library(RColorBrewer)

# PCA subtype analysis from PMID 28697342

source("./opc12/import_opc12.R")
source("./opcBulk/import_opcBulk.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))
bulk_pc <- cbind(bulk$tpm[,1:9],log2(bulk$tpm[,10:ncol(bulk$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]
bulk_pc$symbol <- orthologs$HGNC.symbol[match(x = bulk_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]
bulk_pc <- bulk_pc[match(x = pc$Gene,table = bulk_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
bulk_pc1 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)

opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))
bulk_pc2 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
bulk_projection1 <- colSums(x = bulk_pc1,na.rm = T)

opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)
bulk_projection2 <- colSums(x = bulk_pc2,na.rm = T)

opc12_pc1 <- opc12_pc1[complete.cases(opc12_pc1),]
bulk_pc1 <- bulk_pc1[complete.cases(bulk_pc1),]

opc12_pc1 <- opc12_pc1[,opc12$info$type == "ten-cell"]
bulk_pc1 <- bulk_pc1[,bulk$info$day == 12]

cluster_opc12 <- hclust(d = dist(x = t(opc12_pc1)),method = "ward.D2")
cluster_bulk <- hclust(d = dist(x = t(bulk_pc1)),method = "ward.D2")

pc1_hm <- rbind(t(opc12_pc1[,cluster_opc12$order]),
                t(bulk_pc1[,cluster_bulk$order]))

pdf(file = "./plots/figure_s5.pdf",width = 4,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
pheatmap(mat = pc1_hm[,ncol(pc1_hm):1],
         color = rev(brewer.pal(n = 11,name = "RdBu")),
         breaks = seq(-2,2,length.out = 12),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         gaps_row = rep(56,2))
dev.off()





