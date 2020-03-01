library(annotables)
library(tidyverse)
library(biomaRt)

annotables_mouse <- annotables::grcm38 %>% 
  group_by(ensgene) %>% 
  summarize_all(funs(. %>% unique %>% paste(collapse="; ")))

annotables_human <- annotables::grch38 %>%
  group_by(ensgene) %>%
  summarize_all(funs(. %>% unique %>% paste(collapse="; ")))

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

m2h <- getLDS(attributes = "mgi_symbol",filters = "mgi_symbol",values = annotables_mouse$symbol,mart = mouse,attributesL = "hgnc_symbol",martL = human,uniqueRows = T)
h2m <- getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol",values = annotables_human$symbol,mart = human,attributesL = "mgi_symbol",martL = mouse,uniqueRows = T)

orthologs <- rbind(m2h[,c("MGI.symbol","HGNC.symbol")],
                   h2m[,c("HGNC.symbol","MGI.symbol")])

orthologs <- orthologs[!duplicated(orthologs),]

save(orthologs,file = "./build/orthologs.RData")