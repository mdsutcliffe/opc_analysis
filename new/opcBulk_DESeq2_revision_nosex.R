library(DESeq2)

source("./import/import_opcBulk.R")

opcBulk_DESeq2 <- function() {
  
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
  
  # DESeq2
  .runDESeq2 <- function(rsem,info,day,design) {
    rsem <- rsem[,info$day == day]
    info <- info[info$day == day,]
    
    rsem <- as.matrix(rsem)
    mode(rsem) <- "integer"
    
    modelMatrix <- model.matrix(object = design,data = info)
    
    rsem_filter <- rsem[rowSums(rsem) > 0,]
    
    dds <- DESeqDataSetFromMatrix(countData = rsem_filter,colData = info,design = modelMatrix)
    dds <- DESeq(object = dds)
    
    res <- results(object = dds,contrast = list("genotypeCKO"))
    resOrdered <- res[order(res$padj),]
    
    genesDE <- row.names(resOrdered)[!is.na(resOrdered$padj) & resOrdered$padj < 0.05]
    
    message(paste("design =",
                  format(design),"\n  ",
                  sprintf(fmt = "%6s",day),"dpi :",sprintf(fmt = "%4s",length(genesDE)),"genes"))
    
    return(list("design" = design,
                "modelMatrix" = modelMatrix,
                "genesDE" = genesDE,
                "results" = resOrdered,
                "dds" = dds))
  }
  
  design <- formula(~ genotype)
  
  res12 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 12,design = design)
  res90 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 90,design = design)
  res150 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 150,design = design)
  
  return(list(de12F = res12,
              de90F = res90,
              de150F = res150))
}

bulk$deseq2 <- opcBulk_DESeq2()

save.image("./new/opcBulk_DESeq2_results_revision_nosex.RData")