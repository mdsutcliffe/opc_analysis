# Import 10cRNA-seq 90 dpi samples

import_opc90 <- function() {
  
  source("./functions/normalizeTPM.R")
  
  f.opc90 <- "./data/rsem_opc90.csv"
  f.opc90.info <- "./data/info_opc90.csv"
  
  opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
  opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)
  
  opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))
  
  opc90_tpm_log2 <- cbind(opc90_tpm[,1:9],log2(opc90_tpm[,10:ncol(opc90_tpm)] + 1))
  
  return(list(rsem = opc90,
              info = opc90.info,
              tpm = opc90_tpm,
              log2 = opc90_tpm_log2))
}

opc90 <- import_opc90()