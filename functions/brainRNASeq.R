
# Source
f.barres <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Analysis/data/GSE52564_RAW"

s.barres <- list.files(path = f.barres)
s.barres <- unlist(lapply(X = s.barres,FUN = function(x) strsplit(x = x,split = "_",T)[[1]][2]))
s.barres <- unlist(lapply(X = s.barres,FUN = function(x) strsplit(x = x,split = ".",T)[[1]][1]))

f.s.barres <- list.files(path = f.barres,full.names = T)

lapply(f.s.barres,function(x) read.csv(x))

importBrainRNASeq <- function(f,label) {
  x <-  as.data.frame(readxl::read_xls(path = f,sheet = 1))
  row.names(x) <- x$gene.symbol
  colnames(x)[2] <- label
  return(x)
}

x <- importBrainRNASeq(f = f.s.barres[1],label = s.barres[1])
