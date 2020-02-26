
source("./opc12/RHEGs_opc12.R")
source("./opc90/RHEGs_opc90.R")


write.table(x = opc12_uniqueF,file = "./data/opc12_uniqueF.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc12_uniqueM,file = "./data/opc12_uniqueM.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc12_rheg,file = "./data/opc12_rheg.txt",quote = F,sep = "\n",row.names = F,col.names = F)

write.table(x = opc90_uniqueF,file = "./data/opc90_uniqueF.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc90_uniqueM,file = "./data/opc90_uniqueM.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc90_rheg,file = "./data/opc90_rheg.txt",quote = F,sep = "\n",row.names = F,col.names = F)
