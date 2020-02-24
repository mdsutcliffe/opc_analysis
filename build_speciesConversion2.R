opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

head(listAttributes(human))
getBM(attributes = "hgnc_symbol",filters = )


conversion = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = geneList, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)


cv <- getBM(attributes = "hsapiens_homolog_associated_gene_name",filters = "mgi_symbol",values = opc90$symbol,mart = mouse)

a <- sort(listAttributes(mouse)[,1])
