library(readxl)

load("./build/annotables_mouse.RData")

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
sig_tpm <- cbind(sig[,1,drop = F],apply(X = sig[,2:ncol(sig)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6))

cell_types <- factor(cell_types)
sig_avg <- cbind(sig_tpm[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_tpm[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
row.names(sig_avg) <- sig_avg$symbol
sig_avg <- sig_avg[,2:ncol(sig_avg)]



barplot(log2(as.matrix(sig_avg["Clu",])+1),las = 2,ylim = c(0,12),ylab = "Log2(TPM + 1)",main = "Clusterin expression")

barplot(log2(as.matrix(sig_avg["S100a10",])+1),las = 2,ylim = c(0,8),ylab = "Log2(TPM + 1)",main = "S100a10 expression")
