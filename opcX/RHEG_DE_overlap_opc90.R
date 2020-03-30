source("./opcBulk/import_opcBulk.R")
source("./opc90/import_opc90.R")
source("./opc90/RHEGs_opc90.R")

load("./build/bulk_DESeq2_results.RData")

fList <- list.files(path = "./data/opc90_rdatafiles",full.names = T)
arvListF <- list()
arvListM <- list()
for (i in 1:length(fList)) {
  load(fList[i])
  arvListF[[i]] <- resFemale$arv
  arvListM[[i]] <- resMale$arv
}
gene_set_F <- table(unlist(sapply(1:100,function(x) arvListF[[x]]$symbol)))
gene_set_M <- table(unlist(sapply(1:100,function(x) arvListM[[x]]$symbol)))
gene_set_F <- names(gene_set_F)[gene_set_F == 100]
gene_set_M <- names(gene_set_M)[gene_set_M == 100]

gene_set_opc90 <- intersect(gene_set_F,gene_set_M)
gene_set_opcBulk <- row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),])

common_genes <- intersect(gene_set_opc90,gene_set_opcBulk)

a_b <- length(intersect(opc90_rheg,bulk$deseq2$de90$genesDE))
a_only <- length(setdiff(opc90_rheg,bulk$deseq2$de90$genesDE))
b_only <- length(setdiff(bulk$deseq2$de90$genesDE,opc90_rheg))
neither <- length(setdiff(common_genes,union(opc90_rheg,bulk$deseq2$de90$genesDE)))

fisher.test(x = matrix(data = c(a_b,a_only,b_only,neither),nrow = 2,ncol = 2))

binom.test(x = length(opc90_rheg),
           n = nrow(bulk$rsem),
           p = length(bulk$deseq2$de90$genesDE)/nrow(bulk$rsem),
           alternative = "greater")
