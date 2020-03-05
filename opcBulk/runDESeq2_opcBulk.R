library(DESeq2)

source("./opcBulk/import_opcBulk.R")

runDESeq2 <- function() {
  
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
    
    rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
    
    dds <- DESeqDataSetFromMatrix(countData = rsem_filter,colData = info,design = modelMatrix)
    dds <- DESeq(object = dds,quiet = T,test = "LRT",reduced = model.matrix(object = formula(~ sex),data = info))
    
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
  
  design <- formula(~ genotype + sex)
  res12 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 12,design = design)
  res90 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 90,design = design)
  res150 <- .runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 150,design = design)
  
  return(list(de12 = res12,
              de90 = res90,
              de150 = res150))
}

bulk$deseq2 <- runDESeq2()
# 
# # Plots -----
# 
# fc_12 <- bulk$deseq2$de12$results$log2FoldChange[!is.na(bulk$deseq2$de12$results$pvalue)]
# fc_90 <- bulk$deseq2$de90$results$log2FoldChange[!is.na(bulk$deseq2$de90$results$pvalue)]
# fc_150 <- bulk$deseq2$de150$results$log2FoldChange[!is.na(bulk$deseq2$de150$results$pvalue)]
# 
# p_12 <- bulk$deseq2$de12$results$padj[!is.na(bulk$deseq2$de12$results$pvalue)]
# p_90 <- bulk$deseq2$de90$results$padj[!is.na(bulk$deseq2$de90$results$pvalue)]
# p_150 <- bulk$deseq2$de150$results$padj[!is.na(bulk$deseq2$de150$results$pvalue)]
# 
# # Volcano plots
# {
#   pdf(file = "./plots/volcano_bulk12.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
#   par(mar=c(3.6,3.6,0.5,0.5))
#   plot(x = c(),y = c(),
#        xlim = c(-10,10),ylim = c(0,30),
#        xlab = "",ylab = "",
#        axes = F,xaxs = "i",yaxs = "i")
#   axis(side = 1,lwd = 0.5)
#   axis(side = 2,lwd = 0.5,las = 1)
#   points(x = fc_12[p_12 >= 0.05],
#          y = -log(x = p_12[p_12 >= 0.05],base = 10),
#          col = "#bdbdbd",
#          pch = 16)
#   points(x = fc_12[p_12 < 0.05 & fc_12 > 0],
#          y = -log(x = p_12[p_12 < 0.05 & fc_12 > 0],base = 10),
#          col = "#de2d26",
#          pch = 16)
#   points(x = fc_12[p_12 < 0.05 & fc_12 < 0],
#          y = -log(x = p_12[p_12 < 0.05 & fc_12 < 0],base = 10),
#          col = "#3182bd",
#          pch = 16)
#   legend(x = "topright",
#          legend = c("Increased","Decreased"),
#          col = c("#de2d26","#3182bd"),pch=c(16,16),box.lwd = 0.5)
#   title(xlab = expression("Log"[2]*" fold change"),
#         ylab = expression("-Log"[10]*"("*italic("p")*"-value)"),mgp = c(2.5,1,0))
#   dev.off()
#   
#   pdf(file = "./plots/volcano_bulk90.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
#   par(mar=c(3.6,3.6,0.5,0.5))
#   plot(x = c(),y = c(),
#        xlim = c(-10,10),ylim = c(0,60),
#        xlab = "",ylab = "",
#        axes = F,xaxs = "i",yaxs = "i")
#   axis(side = 1,lwd = 0.5)
#   axis(side = 2,lwd = 0.5,las = 1)
#   points(x = fc_90[p_90 >= 0.05],
#          y = -log(x = p_90[p_90 >= 0.05],base = 10),
#          col = "#bdbdbd",
#          pch = 16)
#   points(x = fc_90[p_90 < 0.05 & fc_90 > 0],
#          y = -log(x = p_90[p_90 < 0.05 & fc_90 > 0],base = 10),
#          col = "#de2d26",
#          pch = 16)
#   points(x = fc_90[p_90 < 0.05 & fc_90 < 0],
#          y = -log(x = p_90[p_90 < 0.05 & fc_90 < 0],base = 10),
#          col = "#3182bd",
#          pch = 16)
#   legend(x = "topright",
#          legend = c("Increased","Decreased"),
#          col = c("#de2d26","#3182bd"),pch=c(16,16),box.lwd = 0.5)
#   title(xlab = expression("Log"[2]*" fold change"),
#         ylab = expression("-Log"[10]*"("*italic("p")*"-value)"),
#         mgp = c(2.5,1,0))
#   dev.off()
#   
#   pdf(file = "./plots/volcano_bulk150.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
#   par(mar=c(3.6,3.6,0.5,0.5))
#   plot(x = c(),y = c(),
#        xlim = c(-10,10),ylim = c(0,100),
#        xlab = "",ylab = "",
#        axes = F,xaxs = "i",yaxs = "i")
#   axis(side = 1,lwd = 0.5)
#   axis(side = 2,lwd = 0.5,las = 1)
#   points(x = fc_150[p_150 >= 0.05],
#          y = -log(x = p_150[p_150 >= 0.05],base = 10),
#          col = "#bdbdbd",
#          pch = 16)
#   points(x = fc_150[p_150 < 0.05 & fc_150 > 0],
#          y = -log(x = p_150[p_150 < 0.05 & fc_150 > 0],base = 10),
#          col = "#de2d26",
#          pch = 16)
#   points(x = fc_150[p_150 < 0.05 & fc_150 < 0],
#          y = -log(x = p_150[p_150 < 0.05 & fc_150 < 0],base = 10),
#          col = "#3182bd",
#          pch = 16)
#   legend(x = "topright",
#          legend = c("Increased","Decreased"),
#          col = c("#de2d26","#3182bd"),pch=c(16,16),box.lwd = 0.5)
#   title(xlab = expression("Log"[2]*" fold change"),
#         ylab = expression("-Log"[10]*"("*italic("p")*"-value)"),
#         mgp = c(2.5,1,0))
#   dev.off()
# }
# 
# # MSigDB enrichments
# 
# library(hypeR)
# 
# msigdb_path <- msigdb_download_all(species = "Mus musculus",output_dir = "./external")
# hallmark <- msigdb_fetch(msigdb_path = msigdb_path,symbol = "H")
# 
# enrich_up_12 <- hypeR(signature = (row.names(bulk$deseq2$de12$results)[!is.na(bulk$deseq2$de12$results$pvalue)])[fc_12 > 0 & p_12 < 0.05],
#                       gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# enrich_down_12 <- hypeR(signature = (row.names(bulk$deseq2$de12$results)[!is.na(bulk$deseq2$de12$results$pvalue)])[fc_12 < 0 & p_12 < 0.05],
#                         gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# # No enrichments
# 
# enrich_up_90 <- hypeR(signature = (row.names(bulk$deseq2$de90$results)[!is.na(bulk$deseq2$de90$results$pvalue)])[fc_90 > 0 & p_90 < 0.05],
#                       gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# enrich_down_90 <- hypeR(signature = (row.names(bulk$deseq2$de90$results)[!is.na(bulk$deseq2$de90$results$pvalue)])[fc_90 < 0 & p_90 < 0.05],
#                         gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# 
# # Only one enrichment for down (angiogenesis)
# 
# enrich_up_150 <- hypeR(signature = (row.names(bulk$deseq2$de150$results)[!is.na(bulk$deseq2$de150$results$pvalue)])[fc_150 > 0 & p_150 < 0.05],
#                       gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# enrich_down_150 <- hypeR(signature = (row.names(bulk$deseq2$de150$results)[!is.na(bulk$deseq2$de150$results$pvalue)])[fc_150 < 0 & p_150 < 0.05],
#                         gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
# 
# hyp_dots(enrich_up_150,top = 100)
# hyp_dots(enrich_down_150,top = 100)
# 
