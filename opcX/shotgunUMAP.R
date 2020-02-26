# Throw everything into UMAP

source("./opcBulk/import_opcBulk.R")
source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")

opc <- cbind(bulk$log2$symbol,
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",10:ncol(opc12$log2)],
             opc90$log2[opc90$log2$chr != "ERCC",10:ncol(opc90$log2)])

