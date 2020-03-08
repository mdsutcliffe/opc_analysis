source("./opcBulk/import_opcBulk.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

bulk_pc <- cbind(bulk$tpm[,1:9],log2(bulk$tpm[,10:ncol(bulk$tpm)] / 100 + 1))

bulk_pc$symbol <- orthologs$HGNC.symbol[match(x = bulk_pc$symbol,table = orthologs$MGI.symbol)]

bulk_pc <- bulk_pc[match(x = pc$Gene,table = bulk_pc$symbol),]

bulk_pc1 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
bulk_pc2 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

bulk_projection1 <- colSums(x = bulk_pc1,na.rm = T)
bulk_projection2 <- colSums(x = bulk_pc2,na.rm = T)


pdf(file = "./plots/pc_150_final.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0))
plot(x = bulk_projection1[bulk$info$day == 150],
     y = bulk_projection2[bulk$info$day == 150],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 150] == "WT",1,16),
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = "PC1",
     ylab = "PC2",
     las = 1,
     lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()



pdf(file = "./plots/pc_12_final.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0))
plot(x = bulk_projection1[bulk$info$day == 12],
     y = bulk_projection2[bulk$info$day == 12],
     pch = ifelse(bulk$info$genotype[bulk$info$day == 12] == "WT",1,16),
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = "PC1",
     ylab = "PC2",
     las = 1,
     lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()




pdf(file = "./plots/pc_12.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = bulk_projection1[bulk$info$genotype == "WT" & bulk$info$day == 12],
       y = bulk_projection2[bulk$info$genotype == "WT" & bulk$info$day == 12],
       pch = 1,lwd = 0.5)
points(x = bulk_projection1[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       y = bulk_projection2[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       pch = 16,lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()



pdf(file = "./plots/pc_12_descriptive_labels.pdf",width = 2,height = 1.75,pointsize = 6,useDingbats = F)
par(mar=c(3.5,5.5,1,3),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,at = c(-60,0,60),labels = c("Proneural","0","Mesenchymal"),lwd = 0.5)
axis(side = 2,at = c(0,140),labels = c("Quiescent","Proliferative"),lwd = 0.5,las = 1)
points(x = bulk_projection1[bulk$info$genotype == "WT" & bulk$info$day == 12],
       y = bulk_projection2[bulk$info$genotype == "WT" & bulk$info$day == 12],
       pch = 1,lwd = 0.5)
points(x = bulk_projection1[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       y = bulk_projection2[bulk$info$genotype == "CKO" & bulk$info$day == 12],
       pch = 16,lwd = 0.5)
legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()
