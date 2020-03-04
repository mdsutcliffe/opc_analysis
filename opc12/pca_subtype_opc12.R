source("./opc12/import_opc12.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)

pdf(file = "./plots/pc_opc12.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar=c(3.25,3.25,1,1),mgp = c(2.25,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,lwd = 0.5)
axis(side = 2,lwd = 0.5,las = 1)
points(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       y = opc12_projection2[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       pch = 1,col = "#e41a1c",lwd = 0.5)
points(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       y = opc12_projection2[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       pch = 16,col = "#377eb8",lwd = 0.5)
legend(x = "topright",legend = c("Female","Male"),col = c("#e41a1c","#377eb8"),pch = c(16,16),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()



pdf(file = "./plots/pc_opc12_descriptive_labels.pdf",width = 2,height = 1.75,pointsize = 6,useDingbats = F)
par(mar=c(3.5,5.5,1,3),mgp = c(2.5,1,0))
plot(x = c(),y = c(),
     xlim = c(-60,60),
     ylim = c(0,140),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "PC2")
axis(side = 1,at = c(-60,0,60),labels = c("Proneural","0","Mesenchymal"),lwd = 0.5)
axis(side = 2,at = c(0,140),labels = c("Quiescent","Proliferative"),lwd = 0.5,las = 1)
points(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       y = opc12_projection2[opc12$info$type == "ten-cell" & opc12$info$sex == "female"],
       pch = 16,col = "#e41a1c",lwd = 0.5)
points(x = opc12_projection1[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       y = opc12_projection2[opc12$info$type == "ten-cell" & opc12$info$sex == "male"],
       pch = 16,col = "#377eb8",lwd = 0.5)
legend(x = "topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),pt.lwd = 0.5,box.lwd = 0.5)
dev.off()
