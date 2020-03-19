library(readxl)

load("./build/annotables_mouse.RData")

f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
f.pericyte <- list.files(path = "./external/GSE75668_RAW",full.names = T)

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

sig_pericyte <- data.frame(ensgene = read.table(file = f.pericyte[1],header = T,sep = "\t",stringsAsFactors = F)[,1])
sig_pericyte$symbol <- annotables_mouse$symbol[match(sig_pericyte$ensgene,annotables_mouse$ensgene)]
sig_pericyte <- cbind(sig_pericyte,do.call(cbind,sapply(X = 1:length(f.pericyte),FUN = function(i) {
  x <- read.table(file = f.pericyte[i],header = T,sep = "\t")[,2,drop = F]
  x[x <= 0.1] <- 0.1
  name_x <- strsplit(x = names(x),split = ".",fixed = T)[[1]][1]
  name_x <- substr(x = name_x,start = 1,stop = nchar(name_x) - 1)
  names(x) <- name_x
  return(x)
})))
cell_types <- c(cell_types,names(sig_pericyte)[3:ncol(sig_pericyte)])
sig_pericyte <- sig_pericyte[complete.cases(sig_pericyte) & !duplicated(sig_pericyte),!colnames(sig_pericyte) %in% "ensgene"]
sig_pericyte_tpm <- cbind(sig_pericyte[,1,drop = F],apply(X = sig_pericyte[,2:ncol(sig_pericyte)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6))

common_signature_genes <- intersect(sig_tpm$symbol,sig_pericyte_tpm$symbol)

sig_join <- cbind(sig_tpm[match(common_signature_genes,sig_tpm$symbol),],sig_pericyte_tpm[match(common_signature_genes,sig_pericyte_tpm$symbol),2:ncol(sig_pericyte_tpm)])

cell_types <- factor(cell_types)

sig_avg <- cbind(sig_join[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_join[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
names(sig_avg)[names(sig_avg) == "Mural"] <- "Pericyte"


removeType <- c("WC","NFO")

sig_genes <- read.table("./external/barres_cell_type.csv",header = T,sep = ",",stringsAsFactors = F)
sig_genes <- sig_genes[,setdiff(names(sig_genes),removeType)]

sig_mat <- sig_avg[match(as.character(unlist(sig_genes)),sig_avg$symbol),] 
sig_mat <- sig_mat[complete.cases(sig_mat),setdiff(names(sig_mat),removeType)]
sig_mat <- sig_mat[!duplicated(sig_mat$symbol),]

write.table(x = sig_mat,file = "./temp/signature_matrix_exact.txt",quote= F,sep = "\t",row.names = F,col.names = T)