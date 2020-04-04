# Figure S7

source("./functions/opc12_rheg.R")

opc12_rheg_mat <- opc12$log2[match(x = opc12$candidates$RHEG,table = opc12$log2$symbol),c(1:9,9+which(opc12$info$type == "ten-cell"))]
row.names(opc12_rheg_mat) <- opc12_rheg_mat$symbol
opc12_rheg_mat <- opc12_rheg_mat[,10:ncol(opc12_rheg_mat)]
opc12_rheg_mat <- as.matrix(opc12_rheg_mat)

pheatmap(mat = opc12_rheg_mat,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         cellwidth = 12,cellheight = 6,
         border_color = NA,
         clustering_method = "ward.D2",
         fontsize = 6,
         scale = "row",
         show_colnames = F,
         filename = "./Figure S7/figure_s7.pdf",
         width = 12,
         height = 6*138/72+1)