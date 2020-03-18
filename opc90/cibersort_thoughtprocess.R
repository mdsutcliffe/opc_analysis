# Mixture file
{
  source("./opc12/import_opc12.R")
  source("./opc90/import_opc90.R")
  source("./opcBulk/import_opcBulk.R")
  
  cibersort_mixture <- cbind(data.frame(symbol = bulk$tpm$symbol),
                             bulk$tpm[,10:ncol(bulk$tpm)],
                             opc12$tpm[!(opc12$tpm$chr == "ERCC"),10:ncol(opc12$tpm)],
                             opc90$tpm[!(opc90$tpm$chr == "ERCC"),10:ncol(opc90$tpm)])
  
  write.table(x = cibersort_mixture,file = "~/mixture.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
}
# Signature matrix file
{
  library(readxl)
  
  load("./build/annotables_mouse.RData")
  
  f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
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
  sig_tpm <- cbind(sig[,1,drop = F],apply(X = sig[,2:ncol(sig)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6))
  
  cell_types <- factor(cell_types)
  
  sig_avg <- cbind(sig_tpm[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_tpm[,1 + which(cell_types == levels(cell_types)[i])]))))
  names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
  
  removeType <- c("WC","NFO","Pericyte")
  
  sig_genes <- read.table("./external/barres_cell_type.csv",header = T,sep = ",",stringsAsFactors = F)
  sig_genes <- sig_genes[,setdiff(names(sig_genes),removeType)]
  
  sig_mat <- sig_avg[match(as.character(unlist(sig_genes)),sig_avg$symbol),] 
  sig_mat <- sig_mat[complete.cases(sig_mat),setdiff(names(sig_mat),removeType)]
  sig_mat <- sig_mat[!duplicated(sig_mat$symbol),]
  
  write.table(x = sig_mat,file = "~/signature.txt",quote= F,sep = "\t",row.names = F,col.names = T)
}


library(pheatmap)
library(RColorBrewer)

f <- "~/CIBERSORTx_batch_corrected.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[173:228,1:(which(names(res) == "P.value") - 1)]
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.275*max(res_opc90)) == 6,"undefined","defined"))
names(ann) <- "class"
pheatmap(res_opc90,annotation_row = ann,color = brewer.pal(9,"Reds"),clustering_method = "ward.D2")
# batch correcting makes most signatures undefined

f <- "~/CIBERSORTx.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[173:228,1:(which(names(res) == "P.value") - 1)]
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.25*max(res_opc90)) == 6,"undefined","defined"))
names(ann) <- "class"
pheatmap(res_opc90,annotation_row = ann,color = brewer.pal(9,"Reds"),clustering_method = "ward.D2")

# run DESeq
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
  
  pheatmap(dehm,scale = "row",annotation_col = ann,clustering_method = "ward.D2",border_color = NA,color = rev(brewer.pal(11,"RdBu")))
}
# DESeq results are totally different


# Average FPKM first, then TPM-normalize
{
  library(readxl)
  
  load("./build/annotables_mouse.RData")
  
  f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
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
  
  cell_types <- factor(cell_types)
  
  sig_avg <- cbind(sig[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig[,1 + which(cell_types == levels(cell_types)[i])]))))
  names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
  
  sig_tpm <- cbind(sig_avg[,1,drop = F],apply(X = sig_avg[,2:ncol(sig_avg)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6))
  
  removeType <- c("WC","NFO","Pericyte")
  
  sig_genes <- read.table("./external/barres_cell_type.csv",header = T,sep = ",",stringsAsFactors = F)
  sig_genes <- sig_genes[,setdiff(names(sig_genes),removeType)]
  
  sig_mat <- sig_tpm[match(as.character(unlist(sig_genes)),sig_tpm$symbol),] 
  sig_mat <- sig_mat[complete.cases(sig_mat),setdiff(names(sig_mat),removeType)]
  sig_mat <- sig_mat[!duplicated(sig_mat$symbol),]
  
  write.table(x = sig_mat,file = "~/signature_avg_before_tpm.txt",quote= F,sep = "\t",row.names = F,col.names = T)
}

f <- "~/CIBERSORTx_avg_before_tpm.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[173:228,1:(which(names(res) == "P.value") - 1)]
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.275*max(res_opc90)) == 6,"undefined","defined"))
names(ann) <- "class"
pheatmap(res_opc90[,c(6,1,2,4,5,3)],annotation_row = ann,color = brewer.pal(9,"Reds"),clustering_method = "ward.D2")
# very similar to previous

# Normalize using intersection genes
{
  source("./opc12/import_opc12.R")
  source("./opc90/import_opc90.R")
  source("./opcBulk/import_opcBulk.R")
  
  library(readxl)
  
  load("./build/annotables_mouse.RData")
  
  f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
  f.barres.base <- basename(f.barres)
  
  sig <- data.frame(symbol = read_xls(path = f.barres[1])$gene.symbol)
  
  cibersort_mixture <- cbind(bulk$rsem[,1:9],
                             bulk$rsem[,10:ncol(bulk$rsem)],
                             opc12$rsem[opc12$rsem$chr != "ERCC",10:ncol(opc12$rsem)],
                             opc90$rsem[opc90$rsem$chr != "ERCC",10:ncol(opc90$rsem)])
  cibersort_mixture <- cibersort_mixture[cibersort_mixture$symbol %in% sig$symbol,]
  cibersort_mixture <- normalizeTPM(rsem = cibersort_mixture,index_counts = 10:ncol(cibersort_mixture))
  cibersort_mixture <- cibersort_mixture[,c(3,10:ncol(cibersort_mixture))]
  write.table(x = cibersort_mixture,file = "~/mixture_renormalize.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  
  sig <- cbind(sig,do.call(cbind,sapply(X = 1:length(f.barres),FUN = function(i) {
    x <- data.frame(read_xls(path = f.barres[i])[,2])
    name_x <- strsplit(x = f.barres.base[i],split = "[_.]")[[1]][2]
    name_x <- substr(x = name_x,start = 1,stop = nchar(name_x)-1)
    names(x) <- name_x
    return(x)
  })))
  cell_types <- names(sig)[2:ncol(sig)]
  sig <- sig[complete.cases(sig) & !duplicated(sig),]
  sig <- sig[sig$symbol %in% cibersort_mixture$symbol,]
  sig_tpm <- cbind(sig[,1,drop = F],apply(X = sig[,2:ncol(sig)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6))
  
  cell_types <- factor(cell_types)
  
  sig_avg <- cbind(sig_tpm[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_tpm[,1 + which(cell_types == levels(cell_types)[i])]))))
  names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
  
  removeType <- c("WC","NFO","Pericyte")
  
  sig_genes <- read.table("./external/barres_cell_type.csv",header = T,sep = ",",stringsAsFactors = F)
  sig_genes <- sig_genes[,setdiff(names(sig_genes),removeType)]
  
  sig_mat <- sig_avg[match(as.character(unlist(sig_genes)),sig_avg$symbol),] 
  sig_mat <- sig_mat[complete.cases(sig_mat),setdiff(names(sig_mat),removeType)]
  sig_mat <- sig_mat[!duplicated(sig_mat$symbol),]
  
  write.table(x = sig_mat,file = "~/signature_renormalize.txt",quote= F,sep = "\t",row.names = F,col.names = T)
}

f <- "~/CIBERSORTx_renormalize.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[173:228,1:(which(names(res) == "P.value") - 1)]
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.3*max(res_opc90)) == 6,"undefined","defined"))
names(ann) <- "class"
pheatmap(res_opc90[,c(6,1,2,4,5,3)],annotation_row = ann,color = brewer.pal(9,"Reds"),clustering_method = "ward.D2")


# Reverting to old out of desperation
f <- "~/CIBERSORTx_opc90_old.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[,1:(which(names(res) == "P.value") - 1)]
ann <- data.frame(ifelse(rowSums(res_opc90 < 0.203*max(res_opc90)) == 6,"undefined","defined"),
                  ifelse(rowSums(res_opc90$OPC > res_opc90[,!(names(res_opc90) %in% "OPC")]) == 5,"OPC-dom",""),
                  ifelse(rowSums(res_opc90$Endothelial > res_opc90[,!names(res_opc90) %in% "Endothelial"]) == 5,"EndoDom","nonEndo"))
names(ann) <- c("class","OPC","Endo")
pdf(file = "./plots/cibersort_opc90_heatmap_4E.pdf",width = 3,height = 2.25,pointsize = 7,family = "ArialMT")
pheatmap(mat = res_opc90[,c(6,1,2,4,5,3)],
         annotation_row = ann,
         color = brewer.pal(9,"Reds"),
         clustering_method = "ward.D2",
         border_color = NA,
         # clustering_distance_rows = "correlation",
         # clustering_distance_cols = "correlation",
         # legend_breaks = NULL,
         # legend_labels = c(""),
         # legend_breaks = c(0),
         lwd = 0.5/0.75,
         annotation_colors = list(class = c(defined = "#00000000",undefined = "#000000FF")),
         treeheight_row = 10,
         treeheight_col = 10,
         show_rownames = F)

dev.off()

ann <- data.frame(E = ifelse(rowSums(res_opc90$Endothelial > res_opc90[,!names(res_opc90) %in% "Endothelial"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6,"Endo",""),
                  O = ifelse(rowSums(res_opc90$OPC > res_opc90[,!names(res_opc90) %in% "OPC"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6,"OPC",""),
                  N = ifelse(rowSums(res_opc90$Neuron > res_opc90[,!names(res_opc90) %in% "Neuron"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6,"Neuron",""),
                  MO = ifelse(rowSums(res_opc90$MO > res_opc90[,!names(res_opc90) %in% "MO"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6,"MO",""),
                  U = ifelse(rowSums(res_opc90 < 0.203*max(res_opc90)) == 6,"U",""))



pheatmap(mat = res_opc90,
         annotation_row = ann,
         color = brewer.pal(9,"Reds"),
         clustering_method = "ward.D2",
         border_color = NA,
         # clustering_distance_rows = "correlation",
         # clustering_distance_cols = "correlation",
         # legend_breaks = NULL,
         # legend_labels = c(""),
         # legend_breaks = c(0),
         lwd = 0.5/0.75,
         
         treeheight_row = 10,
         treeheight_col = 10,
         show_rownames = F)






# Endothelial samples
{
  ann <- data.frame(E = ifelse(rowSums(res_opc90$Endothelial > res_opc90[,!names(res_opc90) %in% "Endothelial"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6,"Endo","Other"))
  pheatmap(mat = res_opc90[,c(6,1,2,4,5,3)],
           annotation_row = ann,
           color = brewer.pal(9,"Reds"),
           clustering_method = "ward.D2",
           border_color = NA,
           # clustering_distance_rows = "correlation",
           # clustering_distance_cols = "correlation",
           # legend_breaks = NULL,
           # legend_labels = c(""),
           # legend_breaks = c(0),
           lwd = 0.5/0.75,
           
           treeheight_row = 10,
           treeheight_col = 10,
           show_rownames = F)
  
  deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
  deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
  deseq_info <- ann
  deseq_info$E <- relevel(deseq_info$E,ref = "Other")
  
  deseq_opc90 <- round(deseq_opc90)
  deseq_opc90 <- as.matrix(deseq_opc90)
  mode(deseq_opc90) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~E)
  dds <- DESeq(object = dds)
  
  res <- results(object = dds)
  resOrdered <- res[order(res$padj),]
  
  degenes <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,])
  dehm <- cbind(data.frame(row.names = degenes),opc90$log2[match(degenes,opc90$log2$symbol),9+which(opc90$info$type == "ten-cell")])
  
  pdf("./plots/heatmap_E_vs_Other.pdf",width = 12,height = 40,pointsize = 6)
  pheatmap(dehm,scale = "row",annotation_col = ann,clustering_method = "ward.D2",border_color = NA,color = rev(brewer.pal(11,"RdBu")),
           show_rownames = T,show_colnames =F)
  dev.off()
  
  write.table(data.frame(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,c("log2FoldChange","padj")]),file = "./temp/E-group_DE_genes.csv",quote = F,row.names = T,col.names = T,sep = ",")

  upreg_endo <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 0,])
  
  sig_genes <- read.table("./external/barres_cell_type.csv",header = T,sep = ",",stringsAsFactors = F)
  upreg_endo %in% sig_genes$Pericyte
  }

sum(ann$OPC == "OPC-dom")/56*100
sum(ann$Endo == "Endo-dom" & ann$class != "undefined")/56*100

sum(rowSums(res_opc90$Neuron > res_opc90[,!names(res_opc90) %in% "Neuron"]) == 5 & rowSums(res_opc90 < 0.203*max(res_opc90)) != 6)/56*100
sum(rowSums(res_opc90 < 0.203*max(res_opc90)) == 6)/56*100

sum(ann$Endo == "EndoDom" & ann$class != "undefined")

annDeEndo <- ann[,"Endo",drop = F] & !annDeEndo[,"class"]


View(opc90$log2[,c(1:9,9+40 + which(ann$Endo == "EndoDom" & ann$class != "undefined"))])

endoHighest <- opc90$log2[,c(1:9,9+40 + which(ann$Endo == "EndoDom" & ann$class != "undefined"))]
endoHighest <- cbind(endoHighest[,1:9],rowMeans(endoHighest[,10:ncol(endoHighest)]))


# run DESeq on Endo
{
  library(DESeq2)
  
  deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
  deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
  deseq_info <- annDeEndo
  deseq_info$Endo <- relevel(deseq_info$Endo,ref = "nonEndo")
  
  deseq_opc90 <- round(deseq_opc90)
  deseq_opc90 <- as.matrix(deseq_opc90)
  mode(deseq_opc90) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~Endo)
  dds <- DESeq(object = dds)
  
  res <- results(object = dds)
  resOrdered <- res[order(res$padj),]
  
  degenes <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,])
  dehm <- cbind(data.frame(row.names = degenes),opc90$log2[match(degenes,opc90$log2$symbol),9+which(opc90$info$type == "ten-cell")])
  
  pheatmap(dehm,scale = "row",annotation_col = annDeEndo,clustering_method = "ward.D2",border_color = NA,color = rev(brewer.pal(11,"RdBu")))
}

# run DESeq
{
  library(DESeq2)
  
  deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
  deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
  deseq_info <- ann[,"class",drop = F]
  deseq_info$class <- relevel(deseq_info$class,ref = "defined")
  
  deseq_opc90 <- round(deseq_opc90)
  deseq_opc90 <- as.matrix(deseq_opc90)
  mode(deseq_opc90) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~class)
  dds <- DESeq(object = dds)
  
  res <- results(object = dds)
  resOrdered <- res[order(res$padj),]
  
  degenes <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,])
  dehm <- cbind(data.frame(row.names = degenes),opc90$log2[match(degenes,opc90$log2$symbol),9+which(opc90$info$type == "ten-cell")])
  
  pdf(file = "./plots/U-group_DE_heatmap.pdf",width = 6,height = 24,pointsize = 7)
  pheatmap(dehm,scale = "row",annotation_col = ann[,"class",drop = F],clustering_method = "ward.D2",border_color = NA,color = rev(brewer.pal(11,"RdBu")))
  dev.off()
}

ann$class

library(vioplot)
ann$class
x <- as.numeric(opc90$log2[opc90$log2$symbol == "Rad51c",9+which(opc90$info$type == "ten-cell")])

pdf(file = "./plots/violin_Rad51c.pdf",width = 2.25/3,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.1,0,0),mgp = c(1.6,0.6,0),xpd = T)
vioplot::vioplot(x[y == "defined"],x[y == "undefined"],drawRect = F,axes = F,las = 1,
                 lwd = 0.5/0.75,
                 ylim = c(0,10),
                 yaxs = "i",xaxs ="i",
                 ylab = "Log2TPM")
set.seed(0)
points(x = runif(n = length(x[y == "defined"]),min = 1-0.25,max = 1+0.25),
       y = x[y == "defined"],pch = 16)
points(x = runif(n = length(x[y == "undefined"]),min = 2-0.25,max = 2+0.25),
       y = x[y == "undefined"],pch = 16)
dev.off()

x <- as.numeric(opc90$log2[opc90$log2$symbol == "Slx1b",9+which(opc90$info$type == "ten-cell")])
pdf(file = "./plots/violin_Slx1b.pdf",width = 2.25/3,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.1,0,0),mgp = c(1.6,0.6,0),xpd = T)
vioplot::vioplot(x[y == "defined"],x[y == "undefined"],drawRect = F,axes = F,las = 1,
                 lwd = 0.5/0.75,
                 ylim = c(0,10),
                 yaxs = "i",xaxs ="i",
                 ylab = "Log2TPM")
set.seed(0)
points(x = runif(n = length(x[y == "defined"]),min = 1-0.25,max = 1+0.25),
       y = x[y == "defined"],pch = 16)
points(x = runif(n = length(x[y == "undefined"]),min = 2-0.25,max = 2+0.25),
       y = x[y == "undefined"],pch = 16)
dev.off()

x <- as.numeric(opc90$log2[opc90$log2$symbol == "Upf3b",9+which(opc90$info$type == "ten-cell")])
pdf(file = "./plots/violin_Upf3b.pdf",width = 2.25/3,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.1,0,0),mgp = c(1.6,0.6,0),xpd = T)
vioplot::vioplot(x[y == "defined"],x[y == "undefined"],drawRect = F,axes = F,las = 1,
                 lwd = 0.5/0.75,
                 ylim = c(0,10),
                 yaxs = "i",xaxs ="i",
                 ylab = "Log2TPM")
set.seed(0)
points(x = runif(n = length(x[y == "defined"]),min = 1-0.25,max = 1+0.25),
       y = x[y == "defined"],pch = 16)
points(x = runif(n = length(x[y == "undefined"]),min = 2-0.25,max = 2+0.25),
       y = x[y == "undefined"],pch = 16)
dev.off()

# Which samples had OPC-dominating?


# Spearman 90 dpi correlations
a <- as.numeric(opc90$log2[opc90$log2$symbol == "Lgr6",9+which(opc90$info$type == "ten-cell")])
b <- as.numeric(opc90$log2[opc90$log2$symbol == "Pim2",9+which(opc90$info$type == "ten-cell")])
cor.test(a,b,method = "spearman")
b <- as.numeric(opc90$log2[opc90$log2$symbol == "Astn1",9+which(opc90$info$type == "ten-cell")])
cor.test(a,b,method = "spearman")
b <- as.numeric(opc90$log2[opc90$log2$symbol == "Smyd2",9+which(opc90$info$type == "ten-cell")])
cor.test(a,b,method = "spearman")
b <- as.numeric(opc90$log2[opc90$log2$symbol == "Mllt11",9+which(opc90$info$type == "ten-cell")])
cor.test(a,b,method = "spearman")

