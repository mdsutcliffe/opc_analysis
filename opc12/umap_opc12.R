
setwd("~")
library(umapr)

setwd("~/Github/opc_analysis/")
source("./opc12/import_opc12.R")

opc12F <- read.table(file = "./data/opc12_true/cleancounts_female_true.tsv")
opc12M <- read.table(file = "./data/opc12_true/cleancounts_male_true.tsv")

umap_opc12 <- cbind(data.frame(row.names = opc12$log2$symbol),opc12$log2[,9+which(opc12$info$type == "ten-cell")])