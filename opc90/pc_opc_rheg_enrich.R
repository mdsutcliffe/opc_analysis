
library(biomaRt)
library(pheat)

source("./opc90/RHEGs_opc90.R")

load("./build/speciesConversion.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:3]

m2h <- convertMouseToHuman(geneList = opc90$tpm$symbol,human = human,mouse = mouse)

opc90_pc <- cbind(opc90$tpm[,1:9],log2(opc90$tpm[,10:ncol(opc90$tpm)]/100 + 1))

# Get orthologous genes
pc <- pc[!is.na(match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_pc$symbol)),]
opc90_pc <- opc90_pc[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_pc$symbol),]

opc90_pc1 <- sapply(X = 10:ncol(opc90_pc),FUN = function(x) pc$PC1 * opc90_pc[,x])
opc90_pc2 <- sapply(X = 10:ncol(opc90_pc),FUN = function(x) pc$PC2 * (opc90_pc[,x] - (pc$PC1 * opc90_pc[,x])))

F_unique <- convertMouseToHuman(geneList = opc90_uniqueF,human = human,mouse = mouse)
M_unique <- convertMouseToHuman(geneList = opc90_uniqueM,human = human,mouse = mouse)
rheg_human <- convertMouseToHuman(geneList = opc90_rheg,human = human,mouse = mouse)

sum(F_unique$HGNC.symbol %in% pc$Gene)
sum(M_unique$HGNC.symbol %in% pc$Gene)

fisher.test(matrix(c(68,638-68,1203-68,15000),nrow = 2,ncol = 2))
row.names(opc90_pc1) <- pc$Gene
colnames(opc90_pc1) <- names(opc90$rsem)[10:ncol(opc90$rsem)]
pc1_opc90_tencell_annotation <- data.frame(row.names = colnames(opc90_pc1)[opc90$info$type == "ten-cell"],
                                           sex = opc90$info$sex[opc90$info$type == "ten-cell"])
colannotation <- data.frame(row.names = pc$Gene,
                            femaleUnique = as.character(pc$Gene %in% F_unique$HGNC.symbol),
                            maleUnique = as.character(pc$Gene %in% M_unique$HGNC.symbol),
                            RHEG = as.character(pc$Gene %in% rheg_human$HGNC.symbol))

pdf(file = "./plots/PC_heatmap_opc90.pdf",width = 8,height = 6)
hm_pc1_opc90 <- pheatmap(mat = t(opc90_pc1[nrow(opc90_pc1):1,opc90$info$type == "ten-cell"]),
                         color = rev(brewer.pal(11,"RdBu")),
                         breaks = seq(-3,3,length.out = 12),
                         border_color = NA,
                         cluster_cols = F,
                         clustering_method = "ward.D2",
                         show_colnames = F,
                         show_rownames = F,
                         annotation_row = pc1_opc90_tencell_annotation,
                         annotation_col = colannotation,annotation_colors = list(femaleUnique = c("TRUE" = "#000000ff","FALSE" = "00000000"),maleUnique = c("TRUE" = "#000000ff","FALSE" = "#00000000"),RHEG = c("TRUE" = "#000000","FALSE" = "#00000000"),sex = c("male" = "#377eb8","female" = "#e41a1c" )))
grid.text(label = "PC1",y = 0.025)
grid.text(label = "10-cell samples",x = 0.875,rot = 90)
dev.off()
