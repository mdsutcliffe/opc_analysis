library(annotables)
library(tidyverse)
library(readxl)

# Cibersort signature matrix
signatureMatrix <- function(geneList = NULL,nGenesEach = 40,removeType = c("WC","NFO"),oligoIndependence = TRUE) {
  load("./build/annotables_mouse.RData")
  
  f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
  f.barres_base <- basename(f.barres)
  
  f.pericyte <- list.files(path = "./external/GSE75668_RAW",full.names = T)
  
  sig <- read_xls(path = f.barres[1])[,1]
  names(sig)[1] <- "symbol"
  sig_samples <- do.call(cbind,sapply(X = 1:length(f.barres),FUN = function(i) read_xls(path = f.barres[i])[,2]))
  
  sig_pericyte <- read.table(file = f.pericyte[1],header = T,sep = "\t",stringsAsFactors = F)[,1]
  sig_pericyte_samples <- do.call(cbind,sapply(X = 1:length(f.pericyte),FUN = function(i) read.table(file = f.pericyte[i],header = T,sep = "\t")[,2,drop = F]))
  sig_pericyte <- cbind(data.frame(ensgene = sig_pericyte),sig_pericyte_samples)
  sig_pericyte <- sig_pericyte[!duplicated(sig_pericyte),]
  
  # Keep only annotated genes
  sig_pericyte$symbol <- annotables_mouse$symbol[match(sig_pericyte$ensgene,annotables_mouse$ensgene)]
  sig_pericyte <- sig_pericyte[complete.cases(sig_pericyte),]
  
  # Reorganize/rename columns
  names(sig_pericyte)[grepl("Mural",names(sig_pericyte))] <- paste0("Mural",1:sum(grepl("Mural",names(sig_pericyte))))
  sig_pericyte <- cbind(sig_pericyte[,"symbol",drop = F],sig_pericyte[,which(grepl("Mural",names(sig_pericyte)))])
  
  # Get signature names
  f.barres_base_cut <- sapply(X = 1:length(f.barres_base),FUN = function(i) strsplit(x = f.barres_base[i],split = "[_.]")[[1]][2])
  sig <- cbind(sig,sig_samples)
  names(sig)[2:ncol(sig)] <- f.barres_base_cut
  
  # Merge pericytes
  sig <- merge(sig,sig_pericyte,by = "symbol")
  
  # Average replicates
  cell_types <- sapply(X = names(sig)[2:ncol(sig)],FUN = function(x) substr(x = x,start = 1,stop = nchar(x)-1),USE.NAMES = F)
  cell_types <- factor(cell_types)
  sig_avg <- cbind(sig[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig[,1 + which(cell_types == levels(cell_types)[i])]))))
  names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
  
  # Normalize against gene intersection
  if (!is.null(geneList)) {
    commonGenes <- intersect(sig_avg$symbol,geneList)
  } else {
    commonGenes <- sig_avg$symbol
  }
  
  # Convert FPKM to TPM for later
  sig_avg_tpm <- sig_avg[match(commonGenes,sig_avg$symbol),]
  sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6)
  sig_avg <- sig_avg[match(commonGenes,sig_avg$symbol),]
  
  sig_Genes <- lapply(X = setdiff(levels(cell_types),removeType),FUN = function(i) {
    iSig <- sig_avg[,c(1,1+which(setdiff(levels(cell_types),removeType) == i))]
    
    iOther <- sig_avg[,c(1,1+which(setdiff(levels(cell_types),removeType) != i))]
    
    if (!oligoIndependence & i %in% c("OPC","NFO","MO")) {
      iOther <- sig_avg[,c(1,1+which(!(levels(cell_types) %in% c("OPC","NFO","MO"))))]
    }
    
    iSig$fc <- iSig[,2] / rowMeans(iOther[,2:ncol(iOther)])
    
    iSig <- iSig[order(iSig$fc,decreasing = T),]
    iSig <- iSig[iSig[,2] > 20,]
    
    return(iSig$symbol[1:1000])
  })
  
  sig_mat <- do.call(cbind,sig_Genes)
  colnames(sig_mat) <- setdiff(levels(cell_types),removeType)
  
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
    allSigGenes <- c(allSigGenes,sig_Genes_unique[1:nGenesEach,i])
  }
  
  sig_avg_tpm_sigGenes <- sig_avg_tpm[match(allSigGenes,sig_avg_tpm$symbol),]
  
  return(list(matrix = sig_avg_tpm_sigGenes,
              tpm = sig_avg_tpm,
              nGenes = nGenesEach,
              genes = allSigGenes))
}
