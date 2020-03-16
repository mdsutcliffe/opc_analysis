source("./opc12/import_opc12.R")


commonGenes <- intersect(opc12$tpm$symbol,sig_avg_tpm$symbol)

opc12$tpm <- normalizeTPM(rsem = opc12$rsem[match(commonGenes,opc12$rsem$symbol),],index_counts = 10:ncol(opc12$rsem))

write.table(x = opc12$tpm[,c(3,10:ncol(opc12$tpm))],file = "./temp/mixture_opc12.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)



f.cibersort <- "./external/CIBERSORTx_opc12.txt"
cibersort <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,2:(which(names(cibersort) == "P.value") - 1)]
cibersort <- cibersort * 100

pheatmap(cibersort[opc12$info$type == "ten-cell",],scale = "column")
