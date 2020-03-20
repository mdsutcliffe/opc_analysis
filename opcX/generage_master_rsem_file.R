# Construct master csv file with all samples

f.opc12 <- "./data/rsem_opc12.csv"
opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)

f.opc90 <- "./data/rsem_opc90.csv"
opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)

# Bulk
{
  # File paths
  bulk.path <- "./data/rsem_opcBulk.csv"
  bulk.info.path <- "./data/info_opcBulk.csv"
  
  # Import
  bulk <- read.csv(bulk.path,stringsAsFactors = F)
  bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)
  
  # Keep only tdT+ samples
  bulk <- bulk[,c(1:9,9+which(bulk.info$celltype == "tdTpositive"))]
  bulk.info <- bulk.info[bulk.info$celltype == "tdTpositive",]
  
  # Bulk had no ERCCs
  bulk_ercc <- cbind(opc12[opc12$chr == "ERCC",1:9],matrix(data = 0,nrow = nrow(opc12[opc12$chr == "ERCC",1:9]),ncol = nrow(bulk.info)))
  names(bulk_ercc)[10:ncol(bulk_ercc)] <- names(bulk)[10:ncol(bulk)]
  
  bulk <- rbind(bulk,bulk_ercc)
}

rsem <- cbind(opc12,opc90[,10:ncol(opc90)],bulk[,10:ncol(bulk)])
names(rsem)[10:ncol(rsem)] <- c(paste0("opc12-",sprintf("%02d",1:96)),
                                paste0("opc90-",sprintf("%02d",1:96)),
                                names(bulk[,10:ncol(bulk)]))

library(tidyverse)
write_csv(x = rsem,path = "./data/rsem_opc.csv")
