
tot_pool <- length(intersect(intersect(row.names(opc90F),row.names(opc90M)),row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),])))

fisher.test(x = matrix(data = c(4,134,465,10000-4-134-465),nrow = 2,ncol = 2))


fList <- list.files(path = "./data/opc90_rdatafiles",full.names = T)

arvListF <- list()
arvListM <- list()
for (i in 1:length(fList)) {
  load(fList[i])
  arvListF[[i]] <- resFemale$arv
  arvListM[[i]] <- resMale$arv
}

gene_set_F <- table(unlist(sapply(1:100,function(x) arvListF[[x]]$symbol)))
gene_set_F <- names(gene_set_F)[gene_set_F == 100]

gene_set_M <- table(unlist(sapply(1:100,function(x) arvListM[[x]]$symbol)))
gene_set_M <- names(gene_set_M)[gene_set_M == 100]

gene_set <- intersect(gene_set_F,gene_set_M)
tot_pool <- length(intersect(gene_set,row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),])))
fisher.test(x = matrix(data = c(4,134,465,tot_pool-4-134-465),nrow = 2,ncol = 2))
