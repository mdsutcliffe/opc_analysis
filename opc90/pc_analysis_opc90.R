library(biomaRt)

source("./opc90/import_opc90.R")

load("./build/speciesConversion.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:3]

h2m <- convertHumanToMouse(geneList = pc$Gene,human = human,mouse = mouse)
m2h <- convertMouseToHuman(geneList = opc90$tpm$symbol,human = human,mouse = mouse)

pc_opc90 <- cbind(opc90$tpm[,1:9],log2(opc90$tpm[,10:ncol(opc90$tpm)]/100 + 1))

pc <- pc[!is.na(match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],pc_opc90$symbol)),]

pc_opc90 <- pc_opc90[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],pc_opc90$symbol),]

pc1_opc90 <- sapply(X = 10:ncol(pc_opc90),FUN = function(x) pc$PC1 * pc_opc90[,x])
pc2_opc90 <- sapply(X = 10:ncol(pc_opc90),FUN = function(x) pc$PC2 * (pc_opc90[,x] - (pc$PC1 * pc_opc90[,x])))

colSums(pc1_opc90)[opc90$info$type == "ten-cell"]

pdf(file = "./plots/pc_opc90.pdf",width = 2,height = 2,pointsize = 6)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(pc1_opc90)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       y = colSums(pc2_opc90)[opc90$info$type == "ten-cell" & opc90$info$sex == "female"],
       col = "#e41a1c",pch = 16)
points(x = colSums(pc1_opc90)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       y = colSums(pc2_opc90)[opc90$info$type == "ten-cell" & opc90$info$sex == "male"],
       col = "#377eb8",pch = 16)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
dev.off()