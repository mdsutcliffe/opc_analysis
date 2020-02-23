library(biomaRt)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertHumanToMouse <- function(geneList,human,mouse) {
  
  conversion = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = geneList, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  
  return(conversion)
  
}

convertMouseToHuman <- function(geneList,human,mouse) {
  
  conversion = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = geneList, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)
  
  return(conversion)
  
}

save.image(file = "./build/speciesConversion.RData")