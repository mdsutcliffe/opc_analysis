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
  .runDESeq2 <- function(rsem,info,day,sex,design) {
    rsem <- rsem[,info$day == day & info$sex == sex]
    info <- info[info$day == day & info$sex == sex,]
    
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
  
  res12F <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 12,sex = "female",design = design)
  res12M <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 12,sex = "male",design = design)
  res90F <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 90,sex = "female",design = design)
  res90M <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 90,sex = "male",design = design)
  res150F <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 150,sex = "female",design = design)
  res150M <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 150,sex = "male",design = design)
  
  return(list(de12F = res12F,
              de12M = res12M,
              de90F = res90F,
              de90M = res90M,
              de150F = res150F,
              de150M = res150M))
}

bulk$deseq2 <- opcBulk_DESeq2()

save.image("./new/opcBulk_DESeq2_results_revision.RData")