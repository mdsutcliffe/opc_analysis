pca_subtype <- function(x) {
  load("./build/orthologs.RData")
  
  pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]
  
  x_pc <- cbind(x$tpm[,1:9],log2(x$tpm[,10:ncol(x$tpm)] / 100 + 1))
  
  x_pc$symbol <- orthologs$HGNC.symbol[match(x = x_pc$symbol,table = orthologs$MGI.symbol)]
  
  x_pc <- x_pc[match(x = pc$Gene,table = x_pc$symbol),]
  
  x_pc1 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
  x_pc2 <- apply(X = x_pc[,10:ncol(x_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))
  
  x_projection1 <- colSums(x = x_pc1,na.rm = T)
  x_projection2 <- colSums(x = x_pc2,na.rm = T)
  
  return(list(pc1 = x_pc1,
              pc2 = x_pc2,
              projection1 = x_projection1,
              projection2 = x_projection2))
}
