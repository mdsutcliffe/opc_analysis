# Get opc12 candidate genes
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)
library(grid)

source("./functions/normalizeTPM.R")

f.opc12 <- "./data/rsem_opc12.csv"
f.opc12.info <- "./data/info_opc12.csv"

opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)

opc12_tpm <- normalizeTPM(rsem = opc12,index_counts = 10:ncol(opc12))

opc12_tpm_log2 <- cbind(opc12_tpm[,1:9],log2(opc12_tpm[,10:ncol(opc12_tpm)] + 1))

res_opc12F <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
res_opc12M <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearances_opc12F <- table(unlist(x = res_opc12F))
nAppearances_opc12M <- table(unlist(x = res_opc12M))

uniqueF_opc12 <- setdiff(x = names(nAppearances_opc12F)[nAppearances_opc12F >= 75],y = names(nAppearances_opc12M)[nAppearances_opc12M >= 75])
uniqueM_opc12 <- setdiff(x = names(nAppearances_opc12M)[nAppearances_opc12M >= 75],y = names(nAppearances_opc12F)[nAppearances_opc12F >= 75])
rheg_opc12 <- intersect(x = names(nAppearances_opc12F)[nAppearances_opc12F >= 75],y = names(nAppearances_opc12M)[nAppearances_opc12M >= 75])

opc12_tpm_log2_rheg <- opc12_tpm_log2[match(x = rheg_opc12,table = opc12_tpm_log2$symbol),]

opc12_tpm_log2_rheg_tencell <- as.matrix(opc12_tpm_log2_rheg[,9+which(opc12.info$type == "ten-cell")])
row.names(opc12_tpm_log2_rheg_tencell) <- opc12_tpm_log2_rheg$symbol

opc12_phm_annotation <- data.frame(row.names = colnames(opc12_tpm_log2_rheg_tencell),
                                   sex = opc12.info$sex[opc12.info$type == "ten-cell"])

pdf(file = "./plots/RHEG_heatmap_opc12.pdf",width = 8,height = 8)
pheatmap(mat = opc12_tpm_log2_rheg_tencell,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         border_color = NA,
         clustering_method = "ward.D2",
         show_rownames = F,
         annotation_col = opc12_phm_annotation,
         scale = "row",
         labels_col = rep(x = "    ",ncol(opc12_tpm_log2_rheg_tencell)))
grid.text(label = "10-cell samples",y = 0.02)
grid.text(label = paste(nrow(opc12_tpm_log2_rheg_tencell),"RHEGs"),x = 0.875,rot = 90)
dev.off()
