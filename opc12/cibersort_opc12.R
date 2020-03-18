source("./opc12/import_opc12.R")


commonGenes <- intersect(opc12$tpm$symbol,sig_avg_tpm$symbol)

opc12$tpm <- normalizeTPM(rsem = opc12$rsem[match(commonGenes,opc12$rsem$symbol),],index_counts = 10:ncol(opc12$rsem))

write.table(x = opc12$tpm[,c(3,10:ncol(opc12$tpm))],file = "./temp/mixture_opc12.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)



f.cibersort <- "./external/CIBERSORTx_opc12.txt"
cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)]
cibersort <- cibersort * 100

pheatmap(cibersort[opc12$info$type == "ten-cell",],scale = "column")





library(RColorBrewer)

source("./opc12/import_opc12.R")
source("./functions/normalizeTPM.R")

f.cibersort <- "./external/CIBERSORTx_opc12.txt"
cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)]
# cibersort <- cibersort * 100

pheatmap(cibersort[opc12$info$type == "ten-cell",],col = brewer.pal(9,"Reds"))

f.cibersort <- "./external/CIBERSORTx_opc12_pericyte.txt"
cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)]
# cibersort <- cibersort * 100

pheatmap(cibersort[opc12$info$type == "ten-cell",],col = brewer.pal(9,"Reds"))

f.cibersort <- "./external/CIBERSORTx_opc12_pericyte_noEndothelial.txt"
cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)]
# cibersort <- cibersort * 100

pheatmap(cibersort[opc12$info$type == "ten-cell",],col = brewer.pal(9,"Reds"))

pdf(file = "./plots/cibersort_opc12_pericyte_noEndo.pdf",width = 3,height = 3,pointsize = 6,useDingbats = F)
par(mar = c(1.8,4,1,8),mgp = c(2.9,1,0),xpd = T)
barplot(t(as.matrix(cibersort[order(cibersort$OPC,decreasing = T),c("OPC","Astrocyte","Neuron","MO","Microglia","Pericyte")])),ylab = "Fraction",names.arg = rep(x = "",nrow(cibersort)),col = rev(brewer.pal(6,"Set1")),las = 1)
title(xlab = "90 dpi ten-cell samples",mgp = c(0.5,0,0))
legend(x = "topright",legend = c("OPC","Astrocyte","Neuron","MO","Microglia","Pericyte"),pch = 15,col = rev(brewer.pal(6,"Set1")),inset = c(-0.4,0))
dev.off()