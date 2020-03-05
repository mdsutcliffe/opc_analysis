
source("opc90/import_opc90.R")
source("./functions/signatureMatrix.R")
source("./functions/normalizeTPM.R")

signature <- signatureMatrix(geneList = opc90$rsem$symbol,removeNFO = T)
commonGenes <- intersect(opc90$tpm$symbol,signature$tpm$symbol)

opc90$tpm <- normalizeTPM(rsem = opc90$rsem[match(commonGenes,opc90$rsem$symbol),],index_counts = 10:ncol(opc90$rsem))

opc90_mixture <- opc90$tpm[,c(1:9,9+which(opc90$info$type == "ten-cell"))]

opc90_mixture <- opc90_mixture[,c(3,10:ncol(opc90_mixture))]

write.table(x = opc90_mixture,file = "./temp/mixture_opc90.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort

f.cibersort_opc90 <- "./external/CIBERSORTx_opc90.txt"

cibersort_opc90 <- read.table(file = f.cibersort_opc90,header = T,sep = "\t")

row.names(cibersort_opc90) <- cibersort_opc90$Mixture

cibersort_opc90 <- cibersort_opc90[,2:(which(names(cibersort_opc90) == "P.value") - 1)]
# cibersort_opc90 <- cibersort_opc90[,2:(which(names(cibersort_opc90) == "P.value") - 1)] / cibersort_opc90$Absolute.score..sig.score.
# cibersort_opc90 <- cibersort_opc90[,2:(which(names(cibersort_opc90) == "P.value") - 1)] / colSums(opc90_mixture[opc90_mixture$symbol %in% signature$genes,2:ncol(opc90_mixture)])

pdf(file = "./plots/cibersort_opc90.pdf",width = 3,height = 3,pointsize = 6,useDingbats = F)
par(mar = c(1.8,4,1,8),mgp = c(2.9,1,0),xpd = T)
barplot(t(as.matrix(cibersort_opc90[order(cibersort_opc90$OPC,decreasing = T),])),ylab = "Fraction",names.arg = rep(x = "",nrow(cibersort_opc90)),col = rev(brewer.pal(6,"Set1")),las = 1)
title(xlab = "90 dpi ten-cell samples",mgp = c(0.5,0,0))
legend(x = "topright",legend = names(cibersort_opc90),pch = 15,col = rev(brewer.pal(6,"Set1")),inset = c(-0.4,0))
dev.off()



library(pheatmap)
library(RColorBrewer)

ann <- data.frame(row.names = row.names(cibersort_opc90),
                  undefined = as.character(rowSums(cibersort_opc90 < max(cibersort_opc90)/5) == 6))

pdf(file = "./plots/cibersort_opc90_results.pdf",width = 8,height = 8,pointsize = 16,useDingbats = F)
pheatmap(mat = cibersort_opc90,
         color = brewer.pal(9,"Reds"),
         show_rownames = F,
         annotation_row = ann,
         annotation_colors = list(undefined = c("TRUE" = "#000000FF","FALSE" = "#00000000")))
grid.text(label = "90 dpi ten-cell samples",x = 0.85,y = 0.4,rot = 90)
dev.off()


all(cibersort_opc90)

undefined_samples <- rowSums(cibersort_opc90 < max(cibersort_opc90)/5) == 6
undefined_samples <- names(undefined_samples)[undefined_samples]

# DE analysis

library(DESeq2)
deseq_opc90 <- cbind(data.frame(row.names = opc90_mixture$symbol),opc90_mixture[,2:ncol(opc90_mixture)])
deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
deseq_info <- data.frame(row.names = names(opc90_mixture[,2:ncol(opc90_mixture)]),
                         celltype = ifelse(rowSums(cibersort_opc90 < max(cibersort_opc90)/5) == 6,"undefined","defined"))
deseq_info$celltype <- relevel(deseq_info$celltype,ref = "defined")

deseq_opc90 <- round(deseq_opc90)
deseq_opc90 <- as.matrix(deseq_opc90)
mode(deseq_opc90) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~celltype)
dds <- DESeq(object = dds)

res <- results(object = dds)
resOrdered <- res[order(res$padj),]

as.data.frame(resOrdered[1:20,c("log2FoldChange","pvalue","padj")])
plotCounts(dds = dds,gene = "Rad51c",intgroup = "celltype")

res2 <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,]
res2

res3 <- data.frame(res2[,c(2,5,6)])

write.csv(x = res3,file = "./temp/opc90_celltypedefined_DEgenes.csv",quote = F)

row.names(res3)

ann <- data.frame(row.names = row.names(cibersort_opc90),
                  undefined = as.character(rowSums(cibersort_opc90 < max(cibersort_opc90)/5) == 6),
                  "Upf3a" = as.numeric(opc90$log2[opc90$log2$symbol == "Upf3a",9+which(opc90$info$type == "ten-cell")]))

pdf(file = "./plots/cibersort_opc90_deseq_heatmap_upf3a.pdf",width = 8,height = 8,pointsize = 16,useDingbats = F)
pheatmap(mat = cbind(data.frame(row.names = opc90$log2$symbol[opc90$log2$symbol %in% row.names(res3)]),opc90$log2[opc90$log2$symbol %in% row.names(res3),9+which(opc90$info$type == "ten-cell")]),#clustering_method = "ward.D2",
         annotation_col = ann,color = rev(brewer.pal(11,"RdBu")),scale = "row",annotation_colors = list(undefined = c("TRUE" = "#000000FF","FALSE" = "#00000000")),show_colnames = F)
dev.off()

plot(as.numeric(opc90$log2[opc90$log2$symbol == "Upf3a",9+which(opc90$info$type == "ten-cell")]),
     as.numeric(opc90$log2[opc90$log2$symbol == "Upf3b",9+which(opc90$info$type == "ten-cell")]))

cor.test(as.numeric(opc90$log2[opc90$log2$symbol == "Upf3a",9+which(opc90$info$type == "ten-cell")]),
     as.numeric(opc90$log2[opc90$log2$symbol == "Upf3b",9+which(opc90$info$type == "ten-cell")]))
