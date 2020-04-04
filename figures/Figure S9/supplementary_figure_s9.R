# Figure S9

source("./functions/opc90_rheg.R")

opc90_rheg_mat <- opc90$log2[match(x = opc90$candidates$RHEG,table = opc90$log2$symbol),c(1:9,9+which(opc90$info$type == "ten-cell"))]
row.names(opc90_rheg_mat) <- opc90_rheg_mat$symbol
opc90_rheg_mat <- opc90_rheg_mat[,10:ncol(opc90_rheg_mat)]
opc90_rheg_mat <- as.matrix(opc90_rheg_mat)

pheatmap(mat = opc90_rheg_mat,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         cellwidth = 12,cellheight = 6,
         border_color = NA,
         clustering_method = "ward.D2",
         fontsize = 6,
         scale = "row",
         show_colnames = F,
         filename = "./Figure S9/figure_s9.pdf",
         width = 12,
         height = 6*231/72+1)