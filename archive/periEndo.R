# Hybrid pericyte/endothelial

library(readxl)
library(pheatmap)

load("./build/annotables_mouse.RData")

removeType <- c("WC","Microvascular","NFO")

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

common_signature_genes <- intersect(sig$symbol,sig_pericyte$symbol)

sig_join <- cbind(sig[match(common_signature_genes,sig$symbol),],sig_pericyte[match(common_signature_genes,sig_pericyte$symbol),2:ncol(sig_pericyte)])

cell_types[cell_types == "Mural" | cell_types == "Endothelial"] <- "PeriEndo"
cell_types <- factor(cell_types)

sig_avg <- cbind(sig_join[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_join[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
# names(sig_avg)[names(sig_avg) == "Mural"] <- "Pericyte"

sig_avg_tpm <- sig_avg[match(common_signature_genes,sig_avg$symbol),]
sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6)
sig_avg_tpm <- sig_avg_tpm[,setdiff(names(sig_avg_tpm),removeType)]
sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) ifelse(x > 1,x,1))

sig_genes <- sapply(X = 2:ncol(sig_avg_tpm),FUN = function(i) {
  iGroup <- sig_avg_tpm[,i]
  iOther <- as.numeric(rowMeans(sig_avg_tpm[,setdiff(2:ncol(sig_avg_tpm),i)]))
  
  ratio <- iGroup / iOther
  ratioOrder <- order(ratio,decreasing = T)
  sig_avg_tpm$symbol[ratioOrder[iGroup[ratioOrder] > 100][1:40]]
})

sig_genes_all <- unique(as.character(sig_genes))

sig_matrix <- sig_avg_tpm[match(sig_genes_all,sig_avg_tpm$symbol),]

pheatmap(cbind(data.frame(row.names = sig_matrix$symbol),sig_matrix[,2:ncol(sig_matrix)]),cluster_rows = F,cluster_cols = F,scale = "row",color = rev(brewer.pal(11,"RdBu")),show_rownames = F)

write.table(x = sig_matrix,file = "./temp/signature_matrix_periendo.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)











f <- "./external/CIBERSORTx_optimized_periEndo.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_bulk12 <- resratio[8:12,]
res_bulk90 <- resratio[33:36,]
res_bulk150 <- resratio[20:25,]
res_bulk150CKO <- resratio[13:19,]

boxplot(res_bulk12[,1:(which(names(res) == "P.value") - 1)])
boxplot(res_bulk90[,1:(which(names(res) == "P.value") - 1)])
boxplot(res_bulk150[,1:(which(names(res) == "P.value") - 1)])
boxplot(res_bulk150CKO[,(1:(which(names(res) == "P.value") - 1))[c(3,2,5,4,1,6)]],horizontal = T,las = 1)

res_opc90 <- res[173:228,1:(which(names(res) == "P.value") - 1)]
pheatmap(res_opc90,color = brewer.pal(9,"Reds"),border_color = NA,clustering_method = "ward.D2")
names(which(rowSums(res_opc90 < 0.25*max(res_opc90)) == 7))
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.25*max(res_opc90)) == 7,"undefined","defined"))
names(ann) <- "class"
pdf(file = "./plots/cibersort_optimized_pericyte.pdf",width = 10,height = 9)
pheatmap(res_opc90,annotation_row = ann,color = brewer.pal(9,"Reds"))
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
