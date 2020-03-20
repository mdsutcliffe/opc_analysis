source("./opc90/RHEGs_opc90.R")

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

f_m <- length(opc90_rheg)
f_only <- length(opc90_uniqueF)
m_only <- length(opc90_uniqueM)
neither <- length(setdiff(gene_set_opc90,Reduce(f = union,x = list(opc90_rheg,opc90_uniqueF,opc90_uniqueM))))

fisher.test(x = matrix(data = c(f_m,f_only,m_only,neither),nrow = 2,ncol = 2))
