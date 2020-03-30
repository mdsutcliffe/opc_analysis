# Import 10cRNA-seq 12 dpi samples

import_opc12 <- function() {
  
  source("./functions/normalizeTPM.R")
  
  f.opc12 <- "./data/rsem_opc12.csv"
  f.opc12.info <- "./data/info_opc12.csv"
  
  opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
  opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)
  
  opc12_tpm <- normalizeTPM(rsem = opc12,index_counts = 10:ncol(opc12))
  
  opc12_tpm_log2 <- cbind(opc12_tpm[,1:9],log2(opc12_tpm[,10:ncol(opc12_tpm)] + 1))
  
  return(list(rsem = opc12,
              info = opc12.info,
              tpm = opc12_tpm,
              log2 = opc12_tpm_log2))
}

opc12 <- import_opc12()