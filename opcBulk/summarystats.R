#summary stats, bulk

source("./opcBulk/import_opcBulk.R")
load("./build/bulk_DESeq2_results.RData")

ratio12wt <- apply(X = bulk$log2[,9+which(bulk$info$day == 12 & bulk$info$genotype == "WT")],MARGIN = 1,FUN = function(x) IQR(x) / median(x))
ratio90wt <- apply(X = bulk$log2[,9+which(bulk$info$day == 90 & bulk$info$genotype == "WT")],MARGIN = 1,FUN = function(x) IQR(x) / median(x))

ratio12cko <- apply(X = bulk$log2[,9+which(bulk$info$day == 12 & bulk$info$genotype == "CKO")],MARGIN = 1,FUN = function(x) IQR(x) / median(x))
ratio90cko <- apply(X = bulk$log2[,9+which(bulk$info$day == 90 & bulk$info$genotype == "CKO")],MARGIN = 1,FUN = function(x) IQR(x) / median(x))


genes_12 <- row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj),])
genes_90 <- row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),])

ratio12wt_filter <- ratio12wt[bulk$log2$symbol %in% genes_12]
ratio90wt_filter <- ratio90wt[bulk$log2$symbol %in% genes_90]
ratio12cko_filter <- ratio12cko[bulk$log2$symbol %in% genes_12]
ratio90cko_filter <- ratio90cko[bulk$log2$symbol %in% genes_90]


boxplot(x = list("12-wt" = ratio12wt_filter[!is.infinite(ratio12wt_filter)],
                 "12-cko" = ratio12cko_filter[!is.infinite(ratio12cko_filter)],
                 "90-wt" = ratio90wt_filter[!is.infinite(ratio90wt_filter)],
                 "90-cko" = ratio90cko_filter[!is.infinite(ratio90cko_filter)]),
        ylim = c(0,10),
        ylab = "IQR / median",las = 1)
boxplot(x = list("12-wt" = ratio12wt_filter[!is.infinite(ratio12wt_filter)],
                 "12-cko" = ratio12cko_filter[!is.infinite(ratio12cko_filter)],
                 "90-wt" = ratio90wt_filter[!is.infinite(ratio90wt_filter)],
                 "90-cko" = ratio90cko_filter[!is.infinite(ratio90cko_filter)]),
        ylim = c(0,1),
        ylab = "IQR / median",las = 1)

b12cko <- bulk$log2[,9+which(bulk$info$day == 12 & bulk$info$genotype == "CKO")]
cor12 <- c()
for (i in 1:(ncol(b12cko)-1)) {
  for (j in (i+1):ncol(b12cko)) {
    cor12 <- c(cor12,cor(b12cko[,i],b12cko[,j]))
  }
}

b90cko <- bulk$log2[,9+which(bulk$info$day == 90 & bulk$info$genotype == "CKO")]
cor90 <- c()
for (i in 1:(ncol(b90cko)-1)) {
  for (j in (i+1):ncol(b90cko)) {
    cor90 <- c(cor90,cor(b90cko[,i],b90cko[,j]))
  }
}
