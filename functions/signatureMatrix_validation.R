# Signature matrix validation by cibersort

source("./functions/signatureMatrix.R")
source("./opcBulk/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)

write.table(x = signature$tpm,file = "./temp/mixture_signature.txt",quote = F,sep = "\t",row.names = F)
write.table(x = signature$matrix,file = "~/signature_old.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort, get results, and place in external folder

f.cibersort <- "./external/CIBERSORTx_signature_validation.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.

pdf(file = "./plots/cibersort_signature_validation.pdf",width = 6,height = 6)
par(mar = c(1.8,4,1,8),mgp = c(2.9,1,0),xpd = T)
barplot(t(as.matrix(cibersort_signature)),ylab = "Fraction",names.arg = rep(x = "",ncol(cibersort_150_wt)),col = rev(brewer.pal(6,"Set1")),las = 1)
title(xlab = "Signatures as mixtures",mgp = c(0.5,0,0))
legend(x = "topright",legend = names(cibersort_signature),pch = 15,col = rev(brewer.pal(6,"Set1")),inset = c(-0.4,0))
dev.off()

cibersort_signature_noMO <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature_noMO) <- cibersort_signature_noMO$Mixture
cibersort_signature_noMO <- cibersort_signature_noMO[c(1:3,5:6),c(2:4,6:7)] / rowSums(cibersort_signature_noMO[c(1:3,5:6),c(2:4,6:7)])
row.names(cibersort_signature_noMO) <- gsub(pattern = "_mixture",replacement = "",x = row.names(cibersort_signature_noMO))

pdf(file = "./plots/purification_bar.pdf",width = 4,height = 3,pointsize = 7)
par(mar = c(3.1,5.2,1,1),mgp = c(2,1,0))
b <- barplot(t(as.matrix(cibersort_signature_noMO[,ncol(cibersort_signature_noMO):1])),
        horiz = T,las = 1,xlab = "Relative proportion",
        axes = F,ylim = c(0,6),lwd = 0.5)
axis(side = 1)
dev.off()

axis(side = 2)


par(mar = c(2,6,1,1))
barplot(t(as.matrix(cibersort_signature[c("OPC_mixture","Astrocyte_mixture","Neuron_mixture","Endothelial_mixture","Microglia_mixture"),c("OPC","Astrocyte","Neuron","Endothelial","Microglia"),])),col = rev(brewer.pal(6,"Set1")),las = 1,horiz = T,ylim = c(0,7))
axis(side = 2,labels = F)

bar
