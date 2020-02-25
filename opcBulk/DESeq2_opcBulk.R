
source("./functions/normalizeTPM.R")
source("./opcBulk/import_opcBulk.R")

# Get counts matrix and remove gene annotation
row.names(bulk$rsem) <- bulk$rsem$symbol
bulk$rsem <- bulk$rsem[,10:ncol(bulk$rsem)]

# Convert all columns to factors
bulk$info[,names(bulk$info)] <- lapply(bulk$info[,names(bulk$info)],factor)

# For DESeq2 compatability
row.names(bulk$info) <- bulk$info$name

# DESeq2 requires integer counts
bulk$rsem <- round(bulk$rsem)

# DESeq2
runDESeq2 <- function(rsem,info,day,design) {
  rsem <- rsem[,info$day == day]
  info <- info[info$day == day,]
  
  rsem <- as.matrix(rsem)
  mode(rsem) <- "integer"
  
  modelMatrix <- model.matrix(object = design,data = info)
  
  rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
  
  dds <- DESeqDataSetFromMatrix(countData = rsem_filter,colData = info,design = modelMatrix)
  dds <- DESeq(object = dds,quiet = T,test = "LRT",reduced = model.matrix(object = formula(~ sex),data = info))
  
  res <- results(object = dds,contrast = list("genotypeCKO"))
  resOrdered <- res[order(res$padj),]
  
  genesDE <- row.names(resOrdered)[!is.na(resOrdered$padj) & resOrdered$padj < 0.05]
  
  message(paste("design =", format(design),"\n  ",sprintf(fmt = "%6s",day),"dpi :",sprintf(fmt = "%4s",length(genesDE)),"genes"))
  
  return(list("design" = design,
              "modelMatrix" = modelMatrix,
              "genesDE" = genesDE,
              "results" = resOrdered,
              "dds" = dds))
}

design <- formula(~ genotype + sex)
res_plus_12 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 12,design = design)
res_plus_90 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 90,design = design)
res_plus_150 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 150,design = design)