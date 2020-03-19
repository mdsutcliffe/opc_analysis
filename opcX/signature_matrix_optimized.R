library(readxl)
library(pheatmap)

load("./build/annotables_mouse.RData")

removeType <- c("WC","Microvascular","NFO","Pericyte")

f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
f.barres.base <- basename(f.barres)

sig <- data.frame(symbol = read_xls(path = f.barres[1])$gene.symbol)
sig <- cbind(sig,do.call(cbind,sapply(X = 1:length(f.barres),FUN = function(i) {
  x <- data.frame(read_xls(path = f.barres[i])[,2])
  name_x <- strsplit(x = f.barres.base[i],split = "[_.]")[[1]][2]
  name_x <- substr(x = name_x,start = 1,stop = nchar(name_x)-1)
  names(x) <- name_x
  return(x)
})))
cell_types <- names(sig)[2:ncol(sig)]
sig <- sig[complete.cases(sig) & !duplicated(sig),]

common_signature_genes <- intersect(sig$symbol,bulk$rsem$symbol)

sig <- sig[match(common_signature_genes,sig$symbol),]
sig <- sig[complete.cases(sig),]

cell_types <- factor(cell_types)

sig_avg <- cbind(sig[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
names(sig_avg)[names(sig_avg) == "Mural"] <- "Pericyte"

sig_avg_tpm <- sig_avg[match(common_signature_genes,sig_avg$symbol),]
sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6)
sig_avg_tpm <- sig_avg_tpm[,setdiff(names(sig_avg_tpm),removeType)]
# sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) ifelse(x > 1,x,1))

sig_genes <- sapply(X = 2:ncol(sig_avg_tpm),FUN = function(i) {
  iGroup <- sig_avg_tpm[,i]
  iOther <- as.numeric(rowMeans(sig_avg_tpm[,setdiff(2:ncol(sig_avg_tpm),i)]))
  
  ratio <- iGroup / iOther
  ratioOrder <- order(ratio,decreasing = T)
  sig_avg_tpm$symbol[ratioOrder[iGroup[ratioOrder] > 100][1:40]]
})

sig_genes_all <- unique(as.character(sig_genes))

sig_matrix <- sig_avg_tpm[match(sig_genes_all,sig_avg_tpm$symbol),]

pheatmap(sig_matrix[,2:ncol(sig_matrix)],cluster_rows = F,cluster_cols = F,scale = "row")

write.table(x = sig_matrix,file = "./temp/signature_matrix_optimized.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)