
f <- "./external/CIBERSORTx_Job34_Results.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[,1:(which(names(res) == "P.value") - 1)]
pheatmap(res_opc90,color = brewer.pal(9,"Reds"),border_color = NA,clustering_method = "ward.D2")
names(which(rowSums(res_opc90 < 0.25*max(res_opc90)) == 7))
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.25*max(res_opc90)) == 6,"undefined","defined"))
names(ann) <- "class"
pdf(file = "./plots/cibersort_optimized_pericyte.pdf",width = 10,height = 9)
pheatmap(res_opc90,annotation_row = ann,color = brewer.pal(9,"Reds"),clustering_method = "ward.D2")
dev.off()

{
  library(DESeq2)
  
  deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
  deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
  deseq_info <- ann
  deseq_info$celltype <- relevel(deseq_info$class,ref = "defined")
  
  deseq_opc90 <- round(deseq_opc90)
  deseq_opc90 <- as.matrix(deseq_opc90)
  mode(deseq_opc90) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~class)
  dds <- DESeq(object = dds)
  
  res <- results(object = dds)
  resOrdered <- res[order(res$padj),]
  
  degenes <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,])
  dehm <- cbind(data.frame(row.names = degenes),opc90$log2[match(degenes,opc90$log2$symbol),9+which(opc90$info$type == "ten-cell")])
  
  pdf(file = "./plots/cibersortDE_optimized_pericyte.pdf",width = 10,height = 9)
  pheatmap(dehm,scale = "row",annotation_col = ann,clustering_method = "ward.D2",border_color = NA,color = rev(brewer.pal(11,"RdBu")))
  dev.off()
}