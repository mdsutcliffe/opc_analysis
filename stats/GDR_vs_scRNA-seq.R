
source("./functions/normalizeTPM.R")

load("./external/gene_detection_rate_processed.RData")

head(rsemcounts_opc_gse60361)

mean(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5))
sd(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5))

mean(colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5))
sd(colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5))

source("./opc12/import_opc12.R")

mean(colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5))
sd(colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5))

ks.test(x = colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5),
       y = colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5),alternative = "g")

ks.test(x = colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5),
       y = colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5))

mean(c(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5),
       colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5)))
sd(c(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5),
     colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5)))
ks.test(x = c(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5),
             colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5)),
       y = colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5),alternative = "greater")

plot(density(colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5),adjust = 0.75))
lines(density(colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5),adjust = 0.75),col = "red")
lines(density(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5),adjust = 0.75),col = "blue")

plot(ecdf(colSums(opc12$tpm[,9+which(opc12$info$type == "ten-cell")] > 5)))
lines(ecdf(colSums(rsemcounts_opc_gse75330[,10:ncol(rsemcounts_opc_gse75330)] > 5)),col = "red")
lines(ecdf(colSums(rsemcounts_opc_gse60361[,10:ncol(rsemcounts_opc_gse60361)] > 5)),col = "blue")
