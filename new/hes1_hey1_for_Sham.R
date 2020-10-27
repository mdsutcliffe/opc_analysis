source("./import/import_opc12.R")
source("./import/import_opc90.R")

info12 <- opc12$info

info12$Hes1 <- as.numeric(opc12$log2[opc12$log2$symbol == "Hes1",10:ncol(opc12$log2)])
info12$Hey1 <- as.numeric(opc12$log2[opc12$log2$symbol == "Hey1",10:ncol(opc12$log2)])
info12$Hes1_plus_Hey1 <- info12$Hes1 + info12$Hey1

info90 <- opc90$info

info90$Hes1 <- as.numeric(opc90$log2[opc90$log2$symbol == "Hes1",10:ncol(opc90$log2)])
info90$Hey1 <- as.numeric(opc90$log2[opc90$log2$symbol == "Hey1",10:ncol(opc90$log2)])
info90$Hes1_plus_Hey1 <- info90$Hes1 + info90$Hey1

write.table(x = info12,file = "./new/opc12_info.tsv",quote = F,row.names = F,col.names = T,sep = "\t")
write.table(x = info90,file = "./new/opc90_info.tsv",quote = F,row.names = F,col.names = T,sep = "\t")
