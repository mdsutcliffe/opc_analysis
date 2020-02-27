library(DESeq2)

source("./functions/normalizeTPM.R")
source("./opcBulk/import_opcBulk.R")

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
res_plus_12 <- runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 12,design = design)
res_plus_90 <- runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 90,design = design)
res_plus_150 <- runDESeq2(rsem = bulk$rsem,info = bulk$info,day = 150,design = design)



# Volcano plot
fc <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj >= 0.05]
p <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj >= 0.05]

fc_up <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange > 0]
p_up <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange > 0]

fc_down <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange < 0]
p_down <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange < 0]





# png("./plots/volcano_150.png",width = 2000,height = 2000,res = 300)
pdf("./plots/volcano_bulk150.pdf",width = 2,height = 2,pointsize = 6)
par(mar=c(4.1,4.1,0.5,0.5))
plot(fc,-log(p,base = 10),
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = expression("Log"[2]*" fold change"),
     ylab = expression("-Log"[10]*"("*italic("p")*"-value)"),
     las = 1,frame = F,lwd = 0.5)
points(fc_up,-log(p_up,base = 10),col = "#de2d26",lwd = 0.5)
points(fc_down,-log(p_down,base = 10),col = "#3182bd",lwd = 0.6)
legend(x = "topright",legend = c("Increased","Decreased"),col = c("#de2d26","#3182bd"),pch=c(1,1),pt.lwd = c(0.5,0.5),box.lwd = 0.5)
dev.off()


msigdb_path <- msigdb_download_all(species = "Mus musculus",output_dir = "./external")
hallmark <- msigdb_fetch(msigdb_path = msigdb_path,symbol = "H")

bulk_150_up <- hypeR(signature = row.names(res_plus_150$results[res_plus_150$results$log2FoldChange > 0 & !is.na(res_plus_150$results$padj) & res_plus_150$results$padj < 0.05,]),gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
pdf("./plots/enrichment_bulk_150_up.pdf",width = 5.5,height = 6,pointsize = 6)
hyp_dots(bulk_150_up,top = 100)
dev.off()

bulk_150_down <- hypeR(signature = row.names(res_plus_150$results[res_plus_150$results$log2FoldChange < 0 & !is.na(res_plus_150$results$padj) & res_plus_150$results$padj < 0.05,]),gsets = hallmark,fdr_cutoff = 0.05,do_plots = T)
pdf("./plots/enrichment_bulk_150_down.pdf",width = 5.5,height = 6,pointsize = 6)
hyp_dots(bulk_150_down,top = 100,)
dev.off()
