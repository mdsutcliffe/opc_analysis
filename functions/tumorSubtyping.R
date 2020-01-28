library(biomaRt)

f.tumorSubtype <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/tumor_subtype.csv"

tumorSubtype <- read.csv(file = f.tumorSubtype,header = T)

convertHumanToMouse <- function(geneList) {
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  conversion = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = geneList, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  return(conversion)
}
  
convertHumanToMouse(geneList = tumorSubtype$GeneSymbol)