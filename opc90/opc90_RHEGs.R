# Get opc90 candidate genes
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)
library(grid)

source("./functions/normalizeTPM.R")

f.opc90 <- "./data/rsem_opc90.csv"
f.opc90.info <- "./data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))

opc90_tpm_log2 <- cbind(opc90_tpm[,1:9],log2(opc90_tpm[,10:ncol(opc90_tpm)] + 1))

res_opc90F <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
res_opc90M <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearances_opc90F <- table(unlist(x = res_opc90F))
nAppearances_opc90M <- table(unlist(x = res_opc90M))

uniqueF_opc90 <- setdiff(x = names(nAppearances_opc90F)[nAppearances_opc90F >= 75],y = names(nAppearances_opc90M)[nAppearances_opc90M >= 75])
uniqueM_opc90 <- setdiff(x = names(nAppearances_opc90M)[nAppearances_opc90M >= 75],y = names(nAppearances_opc90F)[nAppearances_opc90F >= 75])
rheg_opc90 <- intersect(x = names(nAppearances_opc90F)[nAppearances_opc90F >= 75],y = names(nAppearances_opc90M)[nAppearances_opc90M >= 75])

opc90_tpm_log2_rheg <- opc90_tpm_log2[match(x = rheg_opc90,table = opc90_tpm_log2$symbol),]

opc90_tpm_log2_rheg_tencell <- as.matrix(opc90_tpm_log2_rheg[,9+which(opc90.info$type == "ten-cell")])
row.names(opc90_tpm_log2_rheg_tencell) <- opc90_tpm_log2_rheg$symbol

opc90_phm_annotation <- data.frame(row.names = colnames(opc90_tpm_log2_rheg_tencell),
                                   sex = opc90.info$sex[opc90.info$type == "ten-cell"])

pdf(file = "./plots/RHEG_heatmap_opc90.pdf",width = 8,height = 8)
pheatmap(mat = opc90_tpm_log2_rheg_tencell,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         border_color = NA,
         clustering_method = "ward.D2",
         show_rownames = F,
         annotation_col = opc90_phm_annotation,
         scale = "row",
         labels_col = rep(x = "    ",ncol(opc90_tpm_log2_rheg_tencell)))
grid.text(label = "10-cell samples",y = 0.02)
grid.text(label = paste(nrow(opc90_tpm_log2_rheg_tencell),"RHEGs"),x = 0.875,rot = 90)
dev.off()
