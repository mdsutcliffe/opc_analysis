# PCA subtype analysis from PMID 31554641

library(biomaRt)

source("./opc12/import_opc12.R")

load("./build/speciesConversion.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:3]

m2h <- convertMouseToHuman(geneList = opc12$tpm$symbol,human = human,mouse = mouse)

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)]/100 + 1))

# Get orthologous genes
pc <- pc[!is.na(match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc12_pc$symbol)),]
opc12_pc <- opc12_pc[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],opc12_pc$symbol),]

opc12_pc1 <- sapply(X = 10:ncol(opc12_pc),FUN = function(x) pc$PC1 * opc12_pc[,x])
opc12_pc2 <- sapply(X = 10:ncol(opc12_pc),FUN = function(x) pc$PC2 * (opc12_pc[,x] - (pc$PC1 * opc12_pc[,x])))

pdf(file = "./plots/pc_opc12.pdf",width = 1.75,height = 1.75,pointsize = 6,useDingbats = F)
par(mar=c(3.5,3.5,1,1),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(opc12_pc1)[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       y = colSums(opc12_pc2)[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       col = "#e41a1c",pch = 16)
points(x = colSums(opc12_pc1)[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       y = colSums(opc12_pc2)[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       col = "#377eb8",pch = 16)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pc_opc12_descriptive_labels.pdf",width = 2.25,height = 2,pointsize = 6,useDingbats = F)
par(mar=c(3.5,5.5,1,3),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,at = c(-60,60),labels = c("Proneural","Mesenchymal"),lwd = 0.5)
axis(side = 2,at = c(0,140),labels = c("Quiescent","Proliferative"),lwd = 0.5,las = 1)
points(x = colSums(opc12_pc1)[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       y = colSums(opc12_pc2)[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       col = "#e41a1c",pch = 16)
points(x = colSums(opc12_pc1)[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       y = colSums(opc12_pc2)[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       col = "#377eb8",pch = 16)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
dev.off()