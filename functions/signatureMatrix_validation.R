# Signature matrix validation by cibersort

source("./functions/signatureMatrix.R")
source("./opcBulk/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)

write.table(x = signature$tpm,file = "./temp/mixture_signature.txt",quote = F,sep = "\t",row.names = F)
write.table(x = signature$matrix,file = "./temp/signatureMatrix.txt",quote = F,sep = "\t",row.names = F)

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