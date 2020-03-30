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

# Extract fold changes and p-values

fc_12 <- bulk$deseq2$de12$results$log2FoldChange[!is.na(bulk$deseq2$de12$results$padj)]
fc_90 <- bulk$deseq2$de90$results$log2FoldChange[!is.na(bulk$deseq2$de90$results$padj)]
fc_150 <- bulk$deseq2$de150$results$log2FoldChange[!is.na(bulk$deseq2$de150$results$padj)]

p_12 <- bulk$deseq2$de12$results$padj[!is.na(bulk$deseq2$de12$results$padj)]
p_90 <- bulk$deseq2$de90$results$padj[!is.na(bulk$deseq2$de90$results$padj)]
p_150 <- bulk$deseq2$de150$results$padj[!is.na(bulk$deseq2$de150$results$padj)]



# save.image(file = "./build/bulk_DESeq2_results.RData")

# bulk150 volcano plot - horizontally stretched
pdf(file = "./plots/volcano_bulk150_stretch.pdf",width = 3,height = 1.5,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0,0,0),mgp = c(1.6,0.6,0))
plot(x = fc_150[p_150 >= 0.05],
     y = -log10(x = p_150[p_150 >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbd22",
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     lwd = 0.5/0.75)
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,100,20),labels = c(NA,seq(20,100,20)),pos = 0,las = 1,lwd = 0.5/0.75)
points(x = fc_150[p_150 < 0.05 & fc_150 < 0],
       y = -log10(x = p_150[p_150 < 0.05 & fc_150 < 0]),
       pch = 16,cex = 0.5,
       col = "#2166ac22",
       lwd = 0.5/0.75)
points(x = fc_150[p_150 < 0.05 & fc_150 > 0],
       y = -log10(x = p_150[p_150 < 0.05 & fc_150 > 0]),
       pch = 16,cex = 0.5,
       col = "#b2182b22",
       lwd = 0.5/0.75)
text(x = 0.5,y = 100,labels = "Log10(pvalue)",adj = c(0,0.5),xpd = T)
dev.off()

# bulk150 volcano plot
pdf(file = "./plots/volcano_bulk150.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = fc_150[p_150 >= 0.05],
     y = -log10(x = p_150[p_150 >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbda0",
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     lwd = 0.5/0.75)
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,100,20),labels = c(NA,seq(20,100,20)),pos = 0,las = 1,lwd = 0.5/0.75)
points(x = fc_150[p_150 < 0.05 & fc_150 < 0],
       y = -log10(x = p_150[p_150 < 0.05 & fc_150 < 0]),
       pch = 16,cex = 0.5,
       col = "#2166aca0",
       lwd = 0.5/0.75)
points(x = fc_150[p_150 < 0.05 & fc_150 > 0],
       y = -log10(x = p_150[p_150 < 0.05 & fc_150 > 0]),
       pch = 16,cex = 0.5,
       col = "#b2182ba0",
       lwd = 0.5/0.75)
text(x = 0.5,y = 100,labels = "Log10(p-value)",adj = c(0,0.5),xpd = T)
dev.off()

# bulk12 volcano plot
pdf(file = "./plots/volcano_bulk12.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = fc_12[p_12 >= 0.05],
     y = -log10(x = p_12[p_12 >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbda0",
     xlim = c(-10,10),
     ylim = c(0,60),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     lwd = 0.5)
axis(side = 1,lwd = 0.5)
axis(side = 2,at = seq(0,60,20),labels = c(NA,seq(20,60,20)),pos = 0,las = 1,lwd = 0.5)
points(x = fc_12[p_12 < 0.05 & fc_12 < 0],
       y = -log10(x = p_12[p_12 < 0.05 & fc_12 < 0]),
       pch = 16,cex = 0.5,
       col = "#2166aca0",
       lwd = 0.5)
points(x = fc_12[p_12 < 0.05 & fc_12 > 0],
       y = -log10(x = p_12[p_12 < 0.05 & fc_12 > 0]),
       pch = 16,cex = 0.5,
       col = "#b2182ba0",
       lwd = 0.5)
text(x = 0.5,y = 100,labels = "Log10(p-value)",adj = c(0,0.5),xpd = T)
dev.off()

de90 <- bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj) & bulk$deseq2$de90$results$padj < 0.05,c("log2FoldChange","padj")]

# bulk90 volcano plot
pdf(file = "./plots/volcano_bulk90.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = fc_90[p_90 >= 0.05],
     y = -log10(x = p_90[p_90 >= 0.05]),
     col = "#bdbdbda0",
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     lwd = 0.5)
goi90 <- c("Top2a","Ccna2","Hjurp")
col_bulk90 <- ifelse(bulk$deseq2$de90$results$padj[!is.na(bulk$deseq2$de90$results$padj)] >= 0.05,"#bdbdbdaa",
                     ifelse(row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),]) %in% goi90,"#000000ff",
                            ifelse(bulk$deseq2$de90$results$log2FoldChange[!is.na(bulk$deseq2$de90$results$padj)] < 0,"#4393c3a0","#d6604da0")))
plot(x = bulk$deseq2$de90$results$log2FoldChange[!is.na(bulk$deseq2$de90$results$padj)],
     y = -log10(x = bulk$deseq2$de90$results$padj[!is.na(bulk$deseq2$de90$results$padj)]),
     pch = 16,
     col = col_bulk90,
     cex = 0.5
)
axis(side = 1,lwd = 0.5)
axis(side = 2,at = seq(0,100,20),labels = c(NA,seq(20,100,20)),pos = 0,las = 1,lwd = 0.5)
points(x = fc_90[p_90 < 0.05 & fc_90 < 0],
       y = -log10(x = p_90[p_90 < 0.05 & fc_90 < 0]),
       col = "#2166aca0",
       lwd = 0.5)
points(x = fc_90[p_90 < 0.05 & fc_90 > 0],
       y = -log10(x = p_90[p_90 < 0.05 & fc_90 > 0]),
       col = "#b2182ba0",
       lwd = 0.5)
goi <- c("Top2a","Ccna2","Hjurp")
points(x = de90[goi,"log2FoldChange"],
       y = -log10(x = de90[goi,"padj"]),
       pch = 16,cex = 0.5,
       col = "#000000",
       lwd = 0.5)
text(x = 0.5,y = 100,labels = "Log10(pvalue)",adj = c(0,0.5),xpd = T)
dev.off()

de90 <- bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),c("log2FoldChange","padj")]
goi90 <- c("Top2a","Ccna2","Hjurp")
pdf(file = "./plots/volcano_bulk90.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = NULL,y = NULL,
     xlim = c(-10,10),
     ylim = c(0,60),
     xlab = "Log2foldchange",
     ylab = NA,
     xaxs = "i",
     yaxs = "i",
     axes = F)
points(x = de90$log2FoldChange[de90$padj >= 0.05],
       y = -log10(de90$padj[de90$padj >= 0.05]),
       pch = 16,
       cex = 0.5,
       col = "#bdbdbd88")
points(x = de90$log2FoldChange[de90$padj < 0.05 & !(row.names(de90) %in% goi90)],
       y = -log10(de90$padj[de90$padj < 0.05 & !(row.names(de90) %in% goi90)]),
       pch = 16,
       cex = 0.5,
       col = ifelse(de90$log2FoldChange[de90$padj < 0.05 & !(row.names(de90) %in% goi90)] < 0,"#4393c388","#d6604d88"))
points(x = de90$log2FoldChange[de90$padj < 0.05 & (row.names(de90) %in% goi90)],
       y = -log10(de90$padj[de90$padj < 0.05 & (row.names(de90) %in% goi90)]),
       pch = 16,
       cex = 0.5,
       col = "#000000ff")
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,60,20),labels = c(NA,seq(20,60,20)),pos = 0,las = 1,lwd = 0.5/0.75)
text(x = 0.5,y = 60,labels = "Log10pvalue",adj = c(0,0.5),xpd = T)
dev.off()


de12 <- bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05,c("log2FoldChange","padj")]

# bulk12 volcano plot
pdf(file = "./plots/volcano_bulk12_annotated.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = fc_12[p_12 >= 0.05],
     y = -log10(x = p_12[p_12 >= 0.05]),
     pch = 16,cex = 0.5,
     col = "#bdbdbda0",
     xlim = c(-10,10),
     ylim = c(0,60),
     xlab = "Log2 fold change",
     ylab = NA,
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i",
     lwd = 0.5/0.75)
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,60,20),labels = c(NA,seq(20,60,20)),pos = 0,las = 1,lwd = 0.5/0.75)
points(x = de12$log2FoldChange[de12$log2FoldChange < 0],
       y = -log10(x = de12$padj[de12$log2FoldChange < 0]),
       pch = 16,cex = 0.5,
       col = "#4393c3a0",
       lwd = 0.5/0.75)
points(x = de12$log2FoldChange[de12$log2FoldChange > 0],
       y = -log10(x = de12$padj[de12$log2FoldChange > 0]),
       pch = 16,cex = 0.5,
       col = "#d6604da0",
       lwd = 0.5/0.75)
goi <- c("Robo3","Clspn","Mycn","Rpl30","Rpl34","Rplp1","Rps13","Rpl36a")
points(x = de12[goi,"log2FoldChange"],
       y = -log10(x = de12[goi,"padj"]),
       pch = 16,cex = 0.5,
       col = "#000000",
       lwd = 0.5/0.75)
text(x = 0.5,y = 60,labels = "Log10(pvalue)",adj = c(0,0.5),xpd = T)
dev.off()



# MSigDB enrichments

library(hypeR)

msigdb_path <- msigdb_download_all(species = "Mus musculus",output_dir = "./external/MSigDB")
msigdb_path <- list(output_dir = "./external/MSigDB",vs = "V7.0.1")
hallmark <- msigdb_fetch(msigdb_path = msigdb_path,symbol = "H")

enrich_up_12 <- hypeR(signature = (row.names(bulk$deseq2$de12$results)[!is.na(bulk$deseq2$de12$results$pvalue)])[fc_12 > 0 & p_12 < 0.05],
                      gsets = hallmark,fdr_cutoff = 0.05,do_plots = T)
enrich_down_12 <- hypeR(signature = (row.names(bulk$deseq2$de12$results)[!is.na(bulk$deseq2$de12$results$pvalue)])[fc_12 < 0 & p_12 < 0.05],
                        gsets = hallmark,fdr_cutoff = 0.05,do_plots = T)
# No enrichments

enrich_up_90 <- hypeR(signature = (row.names(bulk$deseq2$de90$results)[!is.na(bulk$deseq2$de90$results$pvalue)])[fc_90 > 0 & p_90 < 0.05],
                      gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
enrich_down_90 <- hypeR(signature = (row.names(bulk$deseq2$de90$results)[!is.na(bulk$deseq2$de90$results$pvalue)])[fc_90 < 0 & p_90 < 0.05],
                        gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)

# Only one enrichment for down (angiogenesis)

enrich_up_150 <- hypeR(signature = (row.names(bulk$deseq2$de150$results)[!is.na(bulk$deseq2$de150$results$pvalue)])[fc_150 > 0 & p_150 < 0.05],
                       gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
enrich_down_150 <- hypeR(signature = (row.names(bulk$deseq2$de150$results)[!is.na(bulk$deseq2$de150$results$pvalue)])[fc_150 < 0 & p_150 < 0.05],
                         gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)

hyp_dots(enrich_up_150,top = 100)
hyp_dots(enrich_down_150,top = 100)

hyp_dots(enrich_12_up,top = 100)
hyp_dots(enrich_12_down,top = 100)
hyp_dots(enrich_12_up,top = 100)
hyp_dots(enrich_90_down,top = 100)

# THIS
enrich_12_up <- hypeR(signature = row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange > 0,]),
                      gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 12)] > 0) > 0),
                      fdr_cutoff = 0.01)
enrich_12_down <- hypeR(signature = row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange < 0,]),
                        gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 12)] > 0) > 0),
                        fdr_cutoff = 0.01)
enrich_90_up <- hypeR(signature = row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj) & bulk$deseq2$de90$results$padj < 0.05 & bulk$deseq2$de90$results$log2FoldChange > 0,]),
                      gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 90)] > 0) > 0),
                      fdr_cutoff = 0.01)
enrich_90_down <- hypeR(signature = row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj) & bulk$deseq2$de90$results$padj < 0.05 & bulk$deseq2$de90$results$log2FoldChange < 0,]),
                        gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 90)] > 0) > 0),
                        fdr_cutoff = 0.01)

enrich_150_up <- hypeR(signature = row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange > 0,]),
                       gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0),
                       fdr_cutoff = 0.01)
enrich_150_down <- hypeR(signature = row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange < 0,]),
                         gsets = hallmark,bg = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0),
                         fdr_cutoff = 0.01)

write.table(x = enrich_12_up$as.data.frame(),file = "./temp/msigdb_12_up.tsv",quote = F,sep = "\t",row.names = F)
write.table(x = enrich_12_down$as.data.frame(),file = "./temp/msigdb_12_down.tsv",quote = F,sep = "\t",row.names = F)
write.table(x = enrich_90_up$as.data.frame(),file = "./temp/msigdb_90_up.tsv",quote = F,sep = "\t",row.names = F)
write.table(x = enrich_90_down$as.data.frame(),file = "./temp/msigdb_90_down.tsv",quote = F,sep = "\t",row.names = F)
write.table(x = enrich_150_up$as.data.frame(),file = "./temp/msigdb_150_up.tsv",quote = F,sep = "\t",row.names = F)
write.table(x = enrich_150_down$as.data.frame(),file = "./temp/msigdb_150_down.tsv",quote = F,sep = "\t",row.names = F)

# exact numbers

n150up <- nrow(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange > 0,])
n150down <- nrow(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange < 0,])

binom.test(x = n150up,n = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0),p = n150down/sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0))
binom.test(x = n150down,n = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0),p = n150up/sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0))

binom.test(x = n150up,n = length(bulk$deseq2$de150$genesDE),p = 0.5)
binom.test(x = n150down,n = length(bulk$deseq2$de150$genesDE),p = 0.5)

n12up <- nrow(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange > 0,])
n12down <- nrow(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange < 0,])

binom.test(x = n12up,n = length(bulk$deseq2$de12$genesDE),p = 0.5)
binom.test(x = n12down,n = length(bulk$deseq2$de12$genesDE),p = 0.5)

n90up <- nrow(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj) & bulk$deseq2$de90$results$padj < 0.05 & bulk$deseq2$de90$results$log2FoldChange > 0,])
n90down <- nrow(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj) & bulk$deseq2$de90$results$padj < 0.05 & bulk$deseq2$de90$results$log2FoldChange < 0,])

binom.test(x = n90up,n = length(bulk$deseq2$de90$genesDE),p = 0.5)
binom.test(x = n90down,n = length(bulk$deseq2$de90$genesDE),p = 0.5)

source("./functions/overlap.test.R")
overlap.test(df = data.frame(x = c(rep("bulk12",nrow(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05,])),
                                   rep("bulk150",nrow(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05,]))),
                             y = c(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05,]),
                                   row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05,]))),
             nGenes = sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0 & rowSums(bulk$rsem[,9+which(bulk$info$day == 12)] > 0) > 0),
             nSim = 10000)

fisher.test(x = matrix(data = c(204,
                                nrow(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05,])-204,
                                nrow(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05,])-204,
                                sum(rowSums(bulk$rsem[,9+which(bulk$info$day == 150)] > 0) > 0 & rowSums(bulk$rsem[,9+which(bulk$info$day == 12)] > 0) > 0) - 
                                  nrow(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05,]) - 
                                  nrow(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05,]) + 204),nrow = 2,ncol = 2))

up_12_150 <- intersect(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange > 0,]),
                       row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange > 0,]))

down_12_150 <- intersect(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange < 0,]),
                         row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange < 0,]))

up_12_down_150 <- intersect(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange > 0,]),
                            row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange < 0,]))

down_12_up_150 <- intersect(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj) & bulk$deseq2$de12$results$padj < 0.05 & bulk$deseq2$de12$results$log2FoldChange < 0,]),
                            row.names(bulk$deseq2$de150$results[!is.na(bulk$deseq2$de150$results$padj) & bulk$deseq2$de150$results$padj < 0.05 & bulk$deseq2$de150$results$log2FoldChange > 0,]))

denom <- length(intersect(row.names(bulk$deseq2$de12$results[!is.na(bulk$deseq2$de12$results$padj),]),row.names(bulk$deseq2$de90$results[!is.na(bulk$deseq2$de90$results$padj),])))
binom.test(x = 469,n = denom,p = 143/denom)



denom <- merge(bulk$tpm[,c(3,9 + which(bulk$info$day == 12))],opc12$tpm[,c(3,10:ncol(opc12$tpm))],by = "symbol")
sum(rowSums(denom[,2:ncol(denom)] > 0) > 0)
fisher.test(x = matrix(data = c(4,
                                length(bulk$deseq2$de12$genesDE) - 4,
                                length(opc12_rheg) - 4,
                                sum(rowSums(denom[,2:ncol(denom)] > 0) > 0) - length(union(bulk$deseq2$de12$genesDE,opc12_rheg))),nrow = 2,ncol = 2))


opc12F <- read.table(file = "./data/opc12_true/cleancounts_female_true.tsv")
opc12M <- read.table(file = "./data/opc12_true/cleancounts_male_true.tsv")
