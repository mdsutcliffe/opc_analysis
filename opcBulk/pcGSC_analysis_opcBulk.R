# PC subtype analysis

library(biomaRt)

source("./opcBulk/import_opcBulk.R")
source("./functions/collapseIsoforms.R")

load("./build/speciesConversion.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:3]

h2m <- convertHumanToMouse(geneList = pc$Gene,human = human,mouse = mouse)
m2h <- convertMouseToHuman(geneList = bulk$tpm$symbol,human = human,mouse = mouse)

pc_bulk <- cbind(bulk$tpm[,1:9],log2(bulk$tpm[,10:ncol(bulk$tpm)]/100 + 1))

pc <- pc[!is.na(match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],pc_bulk$symbol)),]

pc_bulk <- pc_bulk[match(m2h$MGI.symbol[match(pc$Gene,m2h$HGNC.symbol)],pc_bulk$symbol),]

pc1_bulk <- sapply(X = 10:ncol(pc_bulk),FUN = function(x) pc$PC1 * pc_bulk[,x])
pc2_bulk <- sapply(X = 10:ncol(pc_bulk),FUN = function(x) pc$PC2 * (pc_bulk[,x] - (pc$PC1 * pc_bulk[,x])))


pdf(file = "./plots/pc_12.pdf",width = 2,height = 2,pointsize = 6)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 12],
       y = colSums(pc2_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 12],
       pch = 1,lwd = 0.5)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       y = colSums(pc2_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       pch = 16,lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pc_90.pdf",width = 2,height = 2,pointsize = 6)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 90],
       y = colSums(pc2_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 90],
       pch = 1,lwd = 0.5)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 90],
       y = colSums(pc2_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 90],
       pch = 16,lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()

pdf(file = "./plots/pc_150.pdf",width = 2,height = 2,pointsize = 6)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 150],
       y = colSums(pc2_bulk)[bulk$info$genotype == "WT" & bulk$info$day == 150],
       pch = 1,lwd = 0.5)
points(x = colSums(pc1_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 150],
       y = colSums(pc2_bulk)[bulk$info$genotype == "CKO" & bulk$info$day == 150],
       pch = 16,lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()
