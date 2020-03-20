library(readxl)
library(pheatmap)

load("./build/annotables_mouse.RData")

removeType <- c("WC","Microvascular","NFO","Pericyte")

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

cell_types <- factor(cell_types)

sig_avg <- cbind(sig_join[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig_join[,1 + which(cell_types == levels(cell_types)[i])]))))
names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
names(sig_avg)[names(sig_avg) == "Mural"] <- "Pericyte"

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

# pheatmap(sig_matrix[,2:ncol(sig_matrix)],cluster_rows = F,cluster_cols = F,scale = "row")
write.table(x = sig_matrix,file = "./Figure S3/signature_matrix_endothelial.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

# Now include pericytes, but remove endothelial
removeType <- c("WC","Microvascular","NFO","Endothelial")

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

sig_matrix_pericyte_noEndothelial <- sig_avg_tpm[match(sig_genes_all,sig_avg_tpm$symbol),]

# pheatmap(sig_matrix_pericyte_noEndothelial[,2:ncol(sig_matrix_pericyte_noEndothelial)],cluster_rows = F,cluster_cols = F,scale = "row")

write.table(x = sig_matrix_pericyte_noEndothelial,file = "./Figure S3/signature_matrix_pericytes.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

# Figure S3 C + E
source("./opcBulk/import_opcBulk.R")

f.cibersort <- "./Figure S3/CIBERSORTx_FigureS3_endothelial.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.
names(cibersort_signature)[names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"

cibersort_bulk150 <- cibersort_signature[c(bulk$info$day == 150 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]

pdf(file = "./plots/figure_s3c.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.6,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk150[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia"))] * 100,lty = 1,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

cibersort_bulk12 <- cibersort_signature[c(bulk$info$day == 12 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]

pdf(file = "./plots/figure_s3e.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.6,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk12[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia"))] * 100,lty = 1,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Endothelial","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

# Figure S3 C + D
source("./opcBulk/import_opcBulk.R")

f.cibersort <- "./Figure S3/CIBERSORTx_FigureS3_pericytes.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.
names(cibersort_signature)[names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"

cibersort_bulk150 <- cibersort_signature[c(bulk$info$day == 150 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]

pdf(file = "./plots/figure_s3d.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.6,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk150[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Pericyte","Microglia"))] * 100,lty = 1,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Pericyte","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()

cibersort_bulk12 <- cibersort_signature[c(bulk$info$day == 12 & bulk$info$genotype == "WT",rep(FALSE,96*2)),]

pdf(file = "./plots/figure_s3f.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.5,0.6,0,0),mgp = c(1.6,0.6,0))
graphics::plot(x = NULL,y = NULL,
               xlim = c(0,80),
               ylim = c(0.5,6.5),
               xlab = "Relative proportion",
               ylab = NA,
               xaxs = "i",
               yaxs = "i",
               axes = F)
graphics::boxplot(x = cibersort_bulk12[,rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Pericyte","Microglia"))] * 100,lty = 1,
                  ylim = c(0,0.8),
                  xaxs = "i",
                  yaxs = "i",
                  horizontal = TRUE,
                  axes = F,
                  lty = 1,
                  lwd = 0.5/0.75,
                  add = T)
axis(side = 1,las = 1,lwd = 0.5/0.75,at = seq(0,80,20),labels = c(0,20,40,60,"80%"))
axis(side = 2,las = 1,lwd = NA,at = 1:6,labels = rev(c("OPC","Astrocyte","Neuron","Myelinating\noligodendrocyte","Pericyte","Microglia")),par(mgp = c(0,0.125,0)))
axis(side = 2,las = 1,lwd = 0.5/0.75,at = 0:6 + 0.5,labels = NA)
dev.off()