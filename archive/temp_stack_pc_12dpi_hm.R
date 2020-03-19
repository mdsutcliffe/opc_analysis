source("./opcBulk/import_opcBulk.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

bulk_pc <- cbind(bulk$tpm[,1:9],log2(bulk$tpm[,10:ncol(bulk$tpm)] / 100 + 1))

bulk_pc$symbol <- orthologs$HGNC.symbol[match(x = bulk_pc$symbol,table = orthologs$MGI.symbol)]

bulk_pc <- bulk_pc[match(x = pc$Gene,table = bulk_pc$symbol),]

bulk_pc1 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
bulk_pc2 <- apply(X = bulk_pc[,10:ncol(bulk_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

bulk_projection1 <- colSums(x = bulk_pc1,na.rm = T)
bulk_projection2 <- colSums(x = bulk_pc2,na.rm = T)

bulk12_pc1_hm <- t(bulk_pc1[complete.cases(bulk_pc1),bulk$info$day == 12])






# PCA subtype analysis from PMID 28697342

source("./opc12/import_opc12.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)

set.seed(0)
nCells <- sample(x = seq(700,1000,10),size = 100,replace = T)
agglomerate_pc <- t(sapply(X = nCells,FUN = function(nC) {
  cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC/10,replace = T)
  
  iOPC12 <- cbind(opc12$tpm[,1:9],log2(rowMeans(opc12$tpm[,cellGroups]) / 100 + 1))
  
  iOPC12$symbol <- orthologs$HGNC.symbol[match(x = iOPC12$symbol,table = orthologs$MGI.symbol)]
  
  iOPC12 <- iOPC12[match(x = pc$Gene,table = iOPC12$symbol),]
  
  iOPC12 <- iOPC12[,10:ncol(iOPC12)]
  
  iPC1 <- pc$PC1 * iOPC12
  iPC2 <- pc$PC2 * (iOPC12 - (pc$PC1 * iOPC12))
  
  iProjection1 <- sum(x = iPC1,na.rm = T)
  iProjection2 <- sum(x = iPC2,na.rm = T)
  
  return(iPC1)
}))
agglomerate_pc_hm <- agglomerate_pc[,!is.na(agglomerate_pc[1,])]

x <- rbind(bulk12_pc1_hm,agglomerate_pc_hm[1:10,])
row.names(x) <- 1:nrow(x)
pdf(file = "./plots/temp_stack_pc12.pdf",width = 10,height = 4)
pheatmap(x[,ncol(x):1],cluster_rows = F,cluster_cols = F,show_rownames = F,color = rev(brewer.pal(11,"RdBu")),breaks = seq(-4,4,length.out = 12),gaps_row = c(7,12),
         annotation_row = data.frame(row.names = row.names(x),group = c(rep("bulkCKO",7),rep("bulkWT",5),rep("sim700-1k",10))),show_colnames = F)
dev.off()
pdf(file = "./plots/temp_stack_pc12.pdf",width = 10,height = 4)
pheatmap(x[,ncol(x):1],cluster_rows = F,cluster_cols = F,show_rownames = F,color = rev(brewer.pal(11,"RdBu")),breaks = seq(-1,1,length.out = 12),gaps_row = c(7,12),
         annotation_row = data.frame(row.names = row.names(x),group = c(rep("bulkCKO",7),rep("bulkWT",5),rep("sim700-1k",10))),show_colnames = F)
dev.off()

