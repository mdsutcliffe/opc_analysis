# Get opc12 candidate genes

library(pheatmap)
library(RColorBrewer)
library(grid)

source("./functions/normalizeTPM.R")
source("./opc12/import_opc12.R")

resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearancesF <- table(unlist(x = resF))
nAppearancesM <- table(unlist(x = resM))

opc12_uniqueF <- setdiff(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
opc12_uniqueM <- setdiff(x = names(nAppearancesM)[nAppearancesM >= 75],y = names(nAppearancesF)[nAppearancesF >= 75])
opc12_rheg <- intersect(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])

opc12_rheg_mat <- opc12$log2[match(x = opc12_rheg,table = opc12$log2$symbol),c(1:9,9+which(opc12$info$type == "ten-cell"))]
row.names(opc12_rheg_mat) <- opc12_rheg_mat$symbol
opc12_rheg_mat <- opc12_rheg_mat[,10:ncol(opc12_rheg_mat)]
opc12_rheg_mat <- as.matrix(opc12_rheg_mat)

# opc12_phm_annotation <- data.frame(row.names = colnames(opc12_rheg_mat),
#                                    sex = opc12$info$sex[opc12$info$type == "ten-cell"])
# 
# pdf(file = "./plots/RHEG_heatmap_opc12.pdf",width = 8,height = 8)
# pheatmap(mat = opc12_rheg_mat,
#          color = rev(brewer.pal(11,"RdBu")),
#          breaks = seq(-4,4,length.out = 12),
#          border_color = NA,
#          clustering_method = "ward.D2",
#          show_rownames = F,
#          annotation_col = opc12_phm_annotation,
#          scale = "row",
#          labels_col = rep(x = "    ",ncol(opc12_rheg_mat)))
# grid.text(label = "10-cell samples",y = 0.02)
# grid.text(label = paste(nrow(opc12_rheg_mat),"RHEGs"),x = 0.875,rot = 90)
# dev.off()
# 
# pdf(file = "./plots/RHEG_heatmap_opc12_inspect.pdf",width = 20,height = 28)
# pheatmap(mat = opc12_rheg_mat,
#          color = rev(brewer.pal(11,"RdBu")),
#          breaks = seq(-4,4,length.out = 12),
#          border_color = NA,
#          clustering_method = "ward.D2",
#          show_rownames = T,
#          annotation_col = opc12_phm_annotation,
#          scale = "row",
#          labels_col = rep(x = "    ",ncol(opc12_rheg_mat)))
# grid.text(label = "10-cell samples",y = 0.005)
# dev.off()

# Correlations
cor(as.numeric(opc12$log2[opc12$log2$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),
    as.numeric(opc12$log2[opc12$log2$symbol == "Vip",9+which(opc12$info$type == "ten-cell")]))

cor(x = as.numeric(opc12$rsem[opc12$rsem$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),
    y = as.numeric(opc12$rsem[opc12$rsem$symbol == "Tgfbr1",9+which(opc12$info$type == "ten-cell")]),
    method = "pearson")

cor.test(x = as.numeric(opc12$tpm[opc12$tpm$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),
    y = as.numeric(opc12$tpm[opc12$tpm$symbol == "Tgfbr2",9+which(opc12$info$type == "ten-cell")]),
    method = "spearman")


cor(scale(as.numeric(opc12$log2[opc12$log2$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")])),
    scale(as.numeric(opc12$log2[opc12$log2$symbol == "Vip",9+which(opc12$info$type == "ten-cell")])))

cor(as.numeric(opc12$rsem[opc12$rsem$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),
    as.numeric(opc12$rsem[opc12$rsem$symbol == "Tgfbr1",9+which(opc12$info$type == "ten-cell")]))


cor(as.numeric(opc12$rsem[opc12$rsem$symbol == "Nbl1",10:ncol(opc12$rsem)]),
    as.numeric(opc12$rsem[opc12$rsem$symbol == "Tgfbr1",10:ncol(opc12$rsem)]))

spearmancors <- apply(X = opc12$tpm[,9+which(opc12$info$type == "ten-cell")],MARGIN = 1,function(x) cor(as.numeric(opc12$tpm[opc12$tpm$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),as.numeric(x)))
plot(x = as.numeric(opc12$log2[opc12$log2$symbol == "Nbl1",9+which(opc12$info$type == "ten-cell")]),
         y = as.numeric(opc12$log2[order(spearmancors,decreasing = T)[3],9+which(opc12$info$type == "ten-cell")]))
