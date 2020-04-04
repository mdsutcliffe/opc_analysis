library(tidyverse)
f.opc12 <- "./data/rsem_opc12.csv"
f.opc12.info <- "./data/info_opc12.csv"

opc12 <- read.csv(file = f.opc12,stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info,stringsAsFactors = F)

opc12_rearrange <- opc12[,c(1:9,9+c(which(opc12.info$type == "pooled" & opc12.info$sex == "female"),
  which(opc12.info$type == "pooled" & opc12.info$sex == "male"),
  which(opc12.info$type == "ten-cell" & opc12.info$sex == "female"),
  which(opc12.info$type == "ten-cell" & opc12.info$sex == "male")))]

names(opc12_rearrange)[10:105] <- paste("opc12",c(rep("pooled",40),rep("10cell",56)),c(rep("female",20),rep("male",20),rep("female",28),rep("male",28)),sep = "_")

write_csv(x = opc12_rearrange,path = "/Volumes/GoogleDrive/My Drive/Profiling manuscripts/OPC datasets/rsem_opc12_rearrange.csv")

opc12_rearrange_tpm <- normalizeTPM(rsem = opc12_rearrange,index_counts = 10:ncol(opc12_rearrange))

write_csv(x = opc12_rearrange_tpm,path = "/Volumes/GoogleDrive/My Drive/Profiling manuscripts/OPC datasets/tpm_opc12_rearrange.csv")

f.opc90 <- "./data/rsem_opc90.csv"
f.opc90.info <- "./data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

opc90_rearrange <- opc90[,c(1:9,9+c(which(opc90.info$type == "pooled" & opc90.info$sex == "female"),
                                    which(opc90.info$type == "pooled" & opc90.info$sex == "male"),
                                    which(opc90.info$type == "ten-cell" & opc90.info$sex == "female"),
                                    which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")))]

names(opc90_rearrange)[10:105] <- paste("opc90",c(rep("pooled",40),rep("10cell",56)),c(rep("female",20),rep("male",20),rep("female",28),rep("male",28)),sep = "_")

write_csv(x = opc90_rearrange,path = "/Volumes/GoogleDrive/My Drive/Profiling manuscripts/OPC datasets/rsem_opc90_rearrange.csv")

opc90_rearrange_tpm <- normalizeTPM(rsem = opc90_rearrange,index_counts = 10:ncol(opc90_rearrange))

write_csv(x = opc90_rearrange_tpm,path = "/Volumes/GoogleDrive/My Drive/Profiling manuscripts/OPC datasets/tpm_opc90_rearrange.csv")
