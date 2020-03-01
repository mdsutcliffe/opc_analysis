# PCA subtype analysis from PMID 31554641

library(biomaRt)
library(pheatmap)
library(RColorBrewer)

source("./opc90/import_opc90.R")
source("./opc90/RHEGs_opc90.R")

load("./build/speciesConversion.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:3]

m2h <- convertMouseToHuman(geneList = opc90$tpm$symbol,human = human,mouse = mouse)

opc90_pc <- cbind(opc90$tpm[,1:9],log2(opc90$tpm[,10:ncol(opc90$tpm)]/100 + 1))

# Get orthologous genes
pc <- pc[!is.na(match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_pc$symbol)),]
opc90_pc <- opc90_pc[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_pc$symbol),]


iDuplicate <- 1
while(any(duplicated(opc90_pc$symbol))) {
        opc90_pc$symbol[duplicated(opc90_pc$symbol)] <- paste(opc90_pc$symbol[duplicated(opc90_pc$symbol)],iDuplicate,sep = "_")
}

row.names(opc90_pc) <- opc90_pc$symbol

opc90_pc1 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc90_pc2 <- apply(X = opc90_pc[,10:ncol(opc90_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

pc1_rheg <- opc90_rheg[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_rheg)]
pc1_uniqueF <- opc90_uniqueF[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_uniqueF)]
pc1_uniqueM <- opc90_uniqueF[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc90_uniqueM)]

pc1_annotation_col <- data.frame(row.names = row.names(opc90_pc1),
                             RHEG = as.character(!is.na(pc1_rheg)),
                             F.specific = as.character(!is.na(pc1_uniqueF)),
                             M.specific = as.character(!is.na(pc1_uniqueM)))
pc1_annotation_row <- data.frame(row.names = colnames(opc90_pc1)[opc90$info$type == "ten-cell"],
                                 sex = opc90$info$sex[opc90$info$type == "ten-cell"])
pc1_annotation_colors <- list(RHEG = c("TRUE" = "#000000","FALSE" = "00000000"),
                           F.specific = c("TRUE" = "#000000","FALSE" = "00000000"),
                           M.specific = c("TRUE" = "#000000","FALSE" = "00000000"),
                           sex = c(female = "#e41a1c",male = "#377eb8"))

pdf(file = "./plots/pc_opc90.pdf",width = 1.75,height = 1.75,pointsize = 6,useDingbats = F)
par(mar=c(3.5,3.5,1,1),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(opc90_pc1)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       y = colSums(opc90_pc2)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       col = "#e41a1c",pch = 16)
points(x = colSums(opc90_pc1)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       y = colSums(opc90_pc2)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       col = "#377eb8",pch = 16)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pc_opc90_descriptive_labels.pdf",width = 2.25,height = 2,pointsize = 6,useDingbats = F)
par(mar=c(3.5,5.5,1,3),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,at = c(-60,60),labels = c("Proneural","Mesenchymal"),lwd = 0.5)
axis(side = 2,at = c(0,140),labels = c("Quiescent","Proliferative"),lwd = 0.5,las = 1)
points(x = colSums(opc90_pc1)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       y = colSums(opc90_pc2)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       col = "#e41a1c",pch = 16)
points(x = colSums(opc90_pc1)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       y = colSums(opc90_pc2)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       col = "#377eb8",pch = 16)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pc_opc90_heatmap.pdf",width = 8,height = 6)
hm_pc1_opc90 <- pheatmap(mat = t(opc90_pc1[nrow(opc90_pc1):1,opc90$info$type == "ten-cell"]),
                         color = rev(brewer.pal(11,"RdBu")),
                         breaks = seq(-3,3,length.out = 12),
                         border_color = NA,
                         cluster_cols = F,
                         clustering_method = "ward.D2",
                         show_colnames = F,
                         show_rownames = F,
                         annotation_col = pc1_annotation_col,
                         annotation_row = pc1_annotation_row,
                         annotation_colors = pc1_annotation_colors)
grid.text(label = "PC1 genes",y = 0.025)
grid.text(label = "10-cell samples",x = 0.86,y = 0.3,rot = 90)
dev.off()

fisher.test(x = matrix(data = c(sum(!is.na(pc1_rheg)),nrow(pc)-sum(!is.na(pc1_rheg)),length(opc90_rheg) - sum(!is.na(pc1_rheg)),nrow(m2h) - (nrow(pc)-sum(!is.na(pc1_rheg))) - (length(opc90_rheg) - sum(!is.na(pc1_rheg)))),nrow = 2,ncol = 2))

fisher.test(x = matrix(data = c(sum(!is.na(pc1_uniqueF)),nrow(pc)-sum(!is.na(pc1_uniqueF)),length(opc90_uniqueF) - sum(!is.na(pc1_uniqueF)),nrow(m2h) - (nrow(pc)-sum(!is.na(pc1_uniqueF))) - (length(opc90_uniqueF) - sum(!is.na(pc1_uniqueF)))),nrow = 2,ncol = 2)) 

fisher.test(x = matrix(data = c(sum(!is.na(pc1_uniqueM)),nrow(pc)-sum(!is.na(pc1_uniqueM)),length(opc90_uniqueM) - sum(!is.na(pc1_uniqueM)),nrow(m2h) - (nrow(pc)-sum(!is.na(pc1_uniqueM))) - (length(opc90_uniqueM) - sum(!is.na(pc1_uniqueM)))),nrow = 2,ncol = 2)) 
