# Get opc90 candidate genes

library(pheatmap)
library(RColorBrewer)
library(grid)

source("./functions/normalizeTPM.R")
source("./opc90/import_opc90.R")

resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearancesF <- table(unlist(x = resF))
nAppearancesM <- table(unlist(x = resM))

opc90_uniqueF <- setdiff(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
opc90_uniqueM <- setdiff(x = names(nAppearancesM)[nAppearancesM >= 75],y = names(nAppearancesF)[nAppearancesF >= 75])
opc90_rheg <- intersect(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])

opc90_rheg_mat <- opc90$log2[match(x = opc90_rheg,table = opc90$log2$symbol),c(1:9,9+which(opc90$info$type == "ten-cell"))]
row.names(opc90_rheg_mat) <- opc90_rheg_mat$symbol
opc90_rheg_mat <- opc90_rheg_mat[,10:ncol(opc90_rheg_mat)]
opc90_rheg_mat <- as.matrix(opc90_rheg_mat)

# opc90_phm_annotation <- data.frame(row.names = colnames(opc90_rheg_mat),
#                                    sex = opc90$info$sex[opc90$info$type == "ten-cell"])
# 
# pdf(file = "./plots/RHEG_heatmap_opc90.pdf",width = 8,height = 8)
# pheatmap(mat = opc90_rheg_mat,
#          color = rev(brewer.pal(11,"RdBu")),
#          breaks = seq(-4,4,length.out = 12),
#          border_color = NA,
#          clustering_method = "ward.D2",
#          show_rownames = F,
#          annotation_col = opc90_phm_annotation,
#          scale = "row",
#          labels_col = rep(x = "    ",ncol(opc90_rheg_mat)))
# grid.text(label = "10-cell samples",y = 0.02)
# grid.text(label = paste(nrow(opc90_rheg_mat),"RHEGs"),x = 0.875,rot = 90)
# dev.off()
# 
# pdf(file = "./plots/RHEG_heatmap_opc90_inspect.pdf",width = 20,height = 28)
# pheatmap(mat = opc90_rheg_mat,
#          color = rev(brewer.pal(11,"RdBu")),
#          breaks = seq(-4,4,length.out = 12),
#          border_color = NA,
#          clustering_method = "ward.D2",
#          show_rownames = T,
#          annotation_col = opc90_phm_annotation,
#          scale = "row",
#          labels_col = rep(x = "    ",ncol(opc90_rheg_mat)))
# grid.text(label = "10-cell samples",y = 0.005)
# dev.off()

