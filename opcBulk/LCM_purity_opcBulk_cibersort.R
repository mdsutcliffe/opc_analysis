# Verify LCM purity by cibersort

source("./opcBulk/import_opcBulk.R")
source("./functions/signatureMatrix.R")
source("./functions/normalizeTPM.R")

signature <- signatureMatrix(geneList = bulk$rsem$symbol,removeNFO = T)
commonGenes <- intersect(bulk$tpm$symbol,signature$tpm$symbol)

bulk$tpm <- normalizeTPM(rsem = bulk$rsem[match(commonGenes,bulk$rsem$symbol),],index_counts = 10:ncol(bulk$rsem))

wt12 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 12 & bulk$info$genotype == "WT"))]
wt150 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"))]

wt12 <- wt12[,c(3,10:ncol(wt12))]
wt150 <- wt150[,c(3,10:ncol(wt150))]

write.table(x = wt12,file = "./temp/mixture_bulk_12_wt.txt",quote = F,sep = "\t",row.names = F)
write.table(x = wt150,file = "./temp/mixture_bulk_150_wt.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort, get results, and place in external folder

