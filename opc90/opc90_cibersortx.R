# 90 dpi deconvolution
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(readxl)

source("./functions/normalizeTPM.R")

f.opc90 <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/rsem_opc90.csv"
f.opc90.info <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))

f <- list.files(path = "./temp/GSE52564_RAW",full.names = T)
f_base <- basename(f)

sig <- read_xls(path = f[1])[,1]
names(sig)[1] <- "symbol"
sig_samples <- do.call(cbind,sapply(X = 1:length(f),FUN = function(i) read_xls(path = f[i])[,2]))

f_base_cut <- sapply(X = 1:length(f_base),FUN = function(i) strsplit(x = f_base[i],split = "[_.]")[[1]][2])
sig <- cbind(sig,sig_samples)
names(sig)[2:ncol(sig)] <- f_base_cut

# remove whole cortex
sig <- sig[,1:15]
f_base_cut <- f_base_cut[1:14]

cell_types <- sapply(X = f_base_cut,FUN = function(x) substr(x = x,start = 1,stop = nchar(x)-1),USE.NAMES = F)
cell_types <- factor(cell_types)
sig_avg <- cbind(sig[,1,drop = F],do.call(cbind,lapply(X = 1:7,FUN = function(i) rowMeans(sig[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)

sig_avg_tpm <- sig_avg
sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6)

common_genes <- intersect(opc90_tpm$symbol,sig_avg_tpm$symbol)

sig_avg_tpm <- sig_avg_tpm[match(common_genes,sig_avg_tpm$symbol),]
opc90_tpm <- opc90_tpm[match(common_genes,opc90_tpm$symbol),]

sig_Genes <- lapply(X = levels(cell_types),FUN = function(i) {
  iSig <- sig_avg_tpm[,c(1,1+which(levels(cell_types) == i))]
  
  if (i %in% c("OPC","NFO","MO")) {
    iOther <- sig_avg_tpm[,c(1,1+which(!(levels(cell_types) %in% c("OPC","NFO","MO"))))]
  } else {
    iOther <- sig_avg_tpm[,c(1,1+which(levels(cell_types) != i))]
  }
  
  iSig$fc <- iSig[,2] / rowMeans(iOther[,2:ncol(iOther)])
  
  iSig <- iSig[order(iSig$fc,decreasing = T),]
  iSig <- iSig[iSig[,2] > 50,]
  
  return(iSig$symbol[1:1000])
})

sig_mat <- do.call(cbind,sig_Genes)
colnames(sig_mat) <- levels(cell_types)

for (i in 2:nrow(sig_mat)) {
  for (j in 1:ncol(sig_mat)) {
    if (sig_mat[i,j] %in% sig_mat[1:(i-1),]) {
      sig_mat[i,j] = NA
    }
    if (j > 1 & sig_mat[i,j] %in% sig_mat[i,1:(j-1)]) {
      sig_mat[i,j] = NA
    }
  }
}

sig_Genes_unique <- apply(sig_mat,2,function(x) {
  y <- na.omit(x)
  return(y[1:40])
})

allSigGenes <- c()
for (i in 1:ncol(sig_Genes_unique)) {
  allSigGenes <- c(allSigGenes,sig_Genes_unique[1:40,i])
}

sig_avg_tpm <- sig_avg_tpm[match(allSigGenes,sig_avg_tpm$symbol),]
opc90_tpm <- opc90_tpm[match(allSigGenes,opc90_tpm$symbol),]

opc90_tpm <- opc90_tpm[,c(3,9+which(opc90.info$type == "ten-cell"))]
names(opc90_tpm)[2:ncol(opc90_tpm)] <- sub("e.","e_",names(opc90_tpm)[2:ncol(opc90_tpm)])

write.table(x = sig_avg_tpm,file = "./temp/signatureMatrix.txt",quote = F,row.names = F,sep = "\t")
write.table(x = opc90_tpm,file = "./temp/mixture_opc90.txt",quote = F,row.names = F,sep = "\t")
# 

temp <- sig_avg_tpm
row.names(temp) <- temp[,1]
temp <- temp[,2:ncol(temp)]
temp <- log2(temp + 1)

png(filename = "./plots/cibersort_input.png",width = 1200,height = 1200,res = 200)
pheatmap(temp[allSigGenes,],cluster_rows = F,cluster_cols = F,scale = "column",color = rev(brewer.pal(11,"RdBu")),breaks = seq(-3,3,length.out = 12))
dev.off()

# results
x <- read.table(file = "~/Downloads/CIBERSORTx_Job2_Results.txt",header = T,sep = "\t")
row.names(x) <- x$Mixture
x <- x[,2:8]
png(filename = "./plots/cibersort_opc90.png",width = 1200,height = 1200,res = 200)
pheatmap(x)
dev.off()


