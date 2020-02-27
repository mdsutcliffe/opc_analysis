source("./opcBulk/import_opcBulk.R")

bulk_cleancounts <- function() {
  # Get counts matrix and remove gene annotation
  row.names(bulk$rsem) <- bulk$rsem$symbol
  bulk$rsem <- bulk$rsem[,10:ncol(bulk$rsem)]
  
  # Convert all columns to factors
  bulk$info[,names(bulk$info)] <- lapply(bulk$info[,names(bulk$info)],factor)
  
  bulk$info$genotype <- relevel(x = bulk$info$genotype,ref = "WT")
  
  # For DESeq2 compatability
  row.names(bulk$info) <- bulk$info$name
  
  # DESeq2 requires integer counts
  bulk$rsem <- round(bulk$rsem)
  
  bulk12 <- bulk$rsem[,which(bulk$info$day == 12)]
  bulk90 <- bulk$rsem[,which(bulk$info$day == 90)]
  bulk150 <- bulk$rsem[,which(bulk$info$day == 150)]
  
  bulk12.info <- bulk$info[which(bulk$info$day == 12),]
  bulk90.info <- bulk$info[which(bulk$info$day == 90),]
  bulk150.info <- bulk$info[which(bulk$info$day == 150),]
  
  bulk12 <- as.matrix(bulk12)
  bulk90 <- as.matrix(bulk90)
  bulk150 <- as.matrix(bulk150)
  
  mode(bulk12) <- "integer"
  mode(bulk90) <- "integer"
  mode(bulk150) <- "integer"
  
  bulk12 <- bulk12[rowSums(bulk12[,bulk12.info$genotype == "WT"]) > 0 & rowSums(bulk12[,bulk12.info$genotype == "CKO"]) > 0,]
  bulk90 <- bulk90[rowSums(bulk90[,bulk90.info$genotype == "WT"]) > 0 & rowSums(bulk90[,bulk90.info$genotype == "CKO"]) > 0,]
  bulk150 <- bulk150[rowSums(bulk150[,bulk150.info$genotype == "WT"]) > 0 & rowSums(bulk150[,bulk150.info$genotype == "CKO"]) > 0,]
  
  return(list(bulk12 = bulk12,
              bulk90 = bulk90,
              bulk150 = bulk150))
}

bulk_clean <- bulk_cleancounts()