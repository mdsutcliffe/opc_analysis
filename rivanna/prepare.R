library(scde)

f.opc12 <- "./rsem_opc12.csv"
f.opc12.info <- "./info_opc12.csv"

nSimulations <- 100

opc12 <- read.csv(file = f.opc12, stringsAsFactors = F)
opc12.info <- read.csv(file = f.opc12.info, stringsAsFactors = F)

nPooled <- sum(opc12.info$type == "pooled")/2
nTenCell <- sum(opc12.info$type == "ten-cell")/2

opc12 <- opc12[opc12$chr != "ERCC" & opc12$chr != "MT",]
row.names(opc12) <- opc12$symbol
opc12 <- opc12[,10:ncol(opc12)]

opc12_clean <- clean.counts(opc12,min.lib.size = 1000, min.reads = 1, min.detected = 1)
opc12_clean <- apply(opc12_clean,2,function(x) {storage.mode(x) <- 'integer'; x})

set.seed(0)

iFemale_pooled <- replicate(nSimulations, {
  i <- sample(which(opc12.info$sex == "female" & opc12.info$type == "pooled"),nPooled)
  i <- i[c(1:(nPooled-1),nPooled-1)]
})

iFemale_tencell <- replicate(nSimulations, {
  i <- sample(which(opc12.info$sex == "female" & opc12.info$type == "ten-cell"),nTenCell)
  i <- i[c(1:(nTenCell-1),nTenCell-1)]
})

iMale_pooled <- replicate(nSimulations, {
  i <- sample(which(opc12.info$sex == "male" & opc12.info$type == "pooled"),nPooled)
  i <- i[c(1:(nPooled-1),nPooled-1)]
})

iMale_tencell <- replicate(nSimulations, {
  i <- sample(which(opc12.info$sex == "male" & opc12.info$type == "ten-cell"),nTenCell)
  i <- i[c(1:(nTenCell-1),nTenCell-1)]
})

save.image(file = "opc12_prepared.RData")