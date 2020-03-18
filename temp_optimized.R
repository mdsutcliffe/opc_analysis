# Cibersort deconvolution

source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")
source("./opcBulk/import_opcBulk.R")

cibersort_mixture <- cbind(data.frame(symbol = bulk$tpm$symbol),
                           bulk$tpm[,10:ncol(bulk$tpm)],
                           opc12$tpm[!(opc12$tpm$chr == "ERCC"),10:ncol(opc12$tpm)],
                           opc90$tpm[!(opc90$tpm$chr == "ERCC"),10:ncol(opc90$tpm)])

write.table(x = cibersort_mixture,file = "~/mixture.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
