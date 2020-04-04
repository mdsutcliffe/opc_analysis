# Supplementary Figure S5

library(pheatmap)
library(RColorBrewer)

source("./import/import_opc12.R")
source("./import/import_opcBulk.R")
source("./functions/pca_subtype.R")

opc12$pc <- pca_subtype(opc12)
bulk$pc <- pca_subtype(bulk)

cluster_opc12 <- hclust(d = dist(x = t(opc12$pc$pc1[complete.cases(opc12$pc$pc1),opc12$info$type == "ten-cell"])),method = "ward.D2")
cluster_bulk <- hclust(d = dist(x = t(bulk$pc$pc1[complete.cases(bulk$pc$pc1),bulk$info$day == 12])),method = "ward.D2")

pc1_hm <- rbind(t(opc12$pc$pc1[complete.cases(opc12$pc$pc1),which(opc12$info$type == "ten-cell")[cluster_opc12$order]]),
                t(bulk$pc$pc1[complete.cases(bulk$pc$pc1),which(bulk$info$day == 12)[cluster_bulk$order]]))

pheatmap(mat = pc1_hm[,ncol(pc1_hm):1],
         color = rev(brewer.pal(n = 11,name = "RdBu")),
         breaks = seq(-2,2,length.out = 12),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         cellwidth = 0.5,
         gaps_row = rep(56,2),
         fontsize = 7,
         filename = "./Figure S5/figure_s5.pdf",
         width = 3.125*1.75,
         height = 3.125)





